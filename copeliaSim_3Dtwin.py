"""
Real-time WEC demo: one plant clock, ESP32 as a monitor, CoppeliaSim as a viewer.

  plant thread : FMU paced to the wall clock (x --speed). Applies the latest Bpto
                 it has received. Single source of truth for waves, buoy, power.
  ESP thread   : sends the LATEST (t, z) to the board as fast as the board answers,
                 reads back "bpto,tp". Slow board = skipped samples, never a queue.
  main thread  : CoppeliaSim free-runs; this thread redraws the plant's current
                 state at ~FPS. All CoppeliaSim calls stay in this thread.

Serial protocol (needs the small firmware change):
  host -> board : "r\n"            board -> host : "ok\n"        (reset, per sea)
  host -> board : "<t>,<z>,<p_mW>\n"   board -> host : "<bpto>,<tp>\n"
  t is plant time in seconds; the firmware computes dt from consecutive t values.

Baseline: --fixed <Bpto> holds the damping constant, no ESP32 needed. Every run
(controller or fixed) appends its mean power per sea to live_results.csv.

Keep the lockstep loop for your evaluation numbers; this one has timing jitter.
"""
import argparse
import os
import threading
import time
from collections import deque
from dataclasses import dataclass

import numpy as np
import pandas as pd
import serial
import yaml
from fmpy import extract, read_model_description
from fmpy.fmi2 import FMU2Slave
from coppeliasim_zmqremoteapi_client import RemoteAPIClient

# ---------------------------------------------------------------- config
FMU = r"G:\Projects\Wave_Energy_Converter\WEC_Heave.fmu"
CAD = r"G:\Projects\Wave_Energy_Converter\wec-twin\cad_out"
# (Tp s, Hs m, duration s): resonance (natural period of WEC_Heave ~0.81 s), mid range, long period
SEAS = [(0.81, 0.04, 10.0), (0.95, 0.04, 20.0), (1.2, 0.04, 20.0)]
H = 0.0005                                            # FMU step (plant time)
VIS = 5.0                                             # visual exaggeration of motion
NX = 25                                               # water columns
FPS = 25.0
G, DEPTH = 9.81, 1.0
B_INIT = 7.0
BPTO_SLEW = 4.0   # max change of the applied damping, per plant second (keeps power from spiking)

# Rough generator model (adjust to your real generator data in params.yaml)
KF = 20.0        # force constant, N/A  (= EMF constant in V per m/s, ideal lossless coupling)
ETA = 1.0        # efficiency applied to the power; 1.0 = the FMU's own P = Bpto*v^2 (same basis as pto_table)
P_WIN = 5.0      # seconds of plant time for the moving-average power
VEL_NAMES = ["zdot", "vz", "v", "heave_vel", "der(z)"]   # FMU velocity output, if it has one

ap = argparse.ArgumentParser()
ap.add_argument("--fmu", default=FMU)
ap.add_argument("--speed", type=float, default=0.5,
                help="plant seconds per wall second (0.25 = slow motion, gives Wokwi more samples per wave)")
ap.add_argument("--port", default="rfc2217://localhost:4000")
ap.add_argument("--fixed", type=float, default=None,
                help="baseline: hold Bpto at this value, no ESP32 needed")
args = ap.parse_args()
B_START = args.fixed if args.fixed is not None else B_INIT


# ---------------------------------------------------------------- waves (same as FMU model)
# The FMU's 20 wave components sit at fixed frequencies wi[1..20] (5.1-8.9 rad/s = Tp 0.71-1.23 s,
# the same band as pto_table.csv). The amplitudes must be computed on exactly that grid.
N, gam = 20, 3.3
_md = read_model_description(args.fmu)
_wi = {v.name: float(v.start) for v in _md.modelVariables if v.name.startswith("wi[")}
w = np.array([_wi[f"wi[{i + 1}]"] for i in range(N)])
dw = w[1] - w[0]
w0, w1 = w[0] - dw / 2, w[-1] + dw / 2                 # band edges, used for the out-of-band warning
ph = np.random.default_rng(1).uniform(0, 2 * np.pi, N)


def amps(Hs, Tp):
    wp = 2 * np.pi / Tp
    sig = np.where(w <= wp, 0.07, 0.09)
    S = w**-5 * np.exp(-1.25 * (wp / w)**4) * gam**np.exp(-(w - wp)**2 / (2 * sig**2 * wp**2))
    a = np.sqrt(2 * S * dw)
    return a * (Hs / 4) / np.sqrt(np.sum(a**2) / 2)


kw = w**2 / G
for _ in range(50):
    kw = w**2 / (G * np.tanh(kw * DEPTH))


# ---------------------------------------------------------------- shared state
@dataclass
class State:
    t: float = 0.0          # plant time of latest state
    z: float = 0.0          # buoy heave at that time
    b: float = B_INIT       # latest damping from the ESP32
    b_ver: int = 0          # bumps when b changes (plant applies on change)
    tp: float = 0.0         # latest period estimate
    v: float = 0.0          # heave velocity
    V: float = 0.0          # instantaneous generator voltage (rough)
    I: float = 0.0          # instantaneous current (rough)
    P: float = 0.0          # instantaneous electrical power
    Vrms: float = 0.0       # moving-window values
    Irms: float = 0.0
    Pavg: float = 0.0
    E: float = 0.0          # cumulative electrical energy, J
    done: bool = False


S = State()
lock = threading.Lock()
stop = threading.Event()


# ---------------------------------------------------------------- threads
def plant_thread(fmu, vr):
    t0 = time.perf_counter()
    t, b_cur = 0.0, B_START
    vname = next((n for n in VEL_NAMES if n in vr), None)
    t_prev, z_prev, v = 0.0, 0.0, 0.0
    E = EV = EI = 0.0                        # integrals of P, V^2, I^2
    hist = deque([(0.0, 0.0, 0.0, 0.0)])     # (t, E, EV, EI) snapshots for the moving window
    next_snap = 0.1
    while not stop.is_set() and t < tsea_cur:
        target = (time.perf_counter() - t0) * args.speed
        with lock:
            b_target = S.b
        # ramp the applied damping toward the ESP32 command instead of jumping
        max_move = BPTO_SLEW * max(target - t, 0.0)
        new_b = b_cur + max(-max_move, min(max_move, b_target - b_cur))
        if abs(new_b - b_cur) > 1e-4:
            b_cur = new_b
            fmu.setReal([vr["Bpto"]], [b_cur])
        b = b_cur
        while t < target and t < tsea_cur:
            fmu.doStep(t, H)
            t += H
        z = fmu.getReal([vr["z"]])[0]
        dtb = t - t_prev
        if dtb > 0:
            v = fmu.getReal([vr[vname]])[0] if vname else (z - z_prev) / dtb
            F = b * v                         # PTO force
            I = F / KF                        # A
            V = KF * v                        # V  (V*I = F*v)
            P = ETA * b * v * v               # W electrical
            E += P * dtb
            EV += V * V * dtb
            EI += I * I * dtb
            t_prev, z_prev = t, z
        else:
            V = I = P = 0.0
        Vrms = Irms = Pavg = 0.0
        if t >= next_snap:
            hist.append((t, E, EV, EI))
            next_snap = t + 0.1
            while len(hist) > 2 and t - hist[1][0] >= P_WIN:
                hist.popleft()
        ta, Ea, EVa, EIa = hist[0]
        if t - ta > 0:
            Pavg = (E - Ea) / (t - ta)
            Vrms = ((EV - EVa) / (t - ta)) ** 0.5
            Irms = ((EI - EIa) / (t - ta)) ** 0.5
        with lock:
            S.t, S.z, S.v, S.V, S.I, S.P = t, z, v, V, I, P
            S.Vrms, S.Irms, S.Pavg, S.E = Vrms, Irms, Pavg, E
        time.sleep(0.001)
    with lock:
        S.done = True


def esp_thread(ser, log):
    last_t = -1.0
    t_wall0 = time.perf_counter()
    while not stop.is_set():
        with lock:
            t, z, p_send = S.t, S.z, S.Pavg
        if t == last_t:                       # nothing new from the plant yet
            time.sleep(0.002)
            continue
        last_t = t
        try:
            # third field = 5 s average power in mW, shown on the OLED only
            ser.write(f"{t:.4f},{z:.6f},{p_send * 1000:.3f}\n".encode())
            deadline = time.perf_counter() + 2.0
            reply = None
            while time.perf_counter() < deadline and not stop.is_set():
                line = ser.readline().decode(errors="ignore").strip()
                try:
                    b, tp = map(float, line.split(","))
                    reply = (b, tp)
                    break
                except ValueError:
                    continue
            if reply is None:
                continue                      # sample lost, move on with the next latest
        except serial.SerialException as e:
            print("serial error:", e)
            break
        b, tp = reply
        with lock:
            t_now = S.t
            if b != S.b:
                S.b, S.b_ver = b, S.b_ver + 1
            S.tp = tp
            Vr, Ir, Pa = S.Vrms, S.Irms, S.Pavg
        log.append((t, z, b, tp, t_now - t, time.perf_counter() - t_wall0, Vr, Ir, Pa))


# ---------------------------------------------------------------- CoppeliaSim scene
client = RemoteAPIClient()
sim = client.require("sim")
sim.setStepping(False)                        # free-running, we only draw


def setpos(h, p):
    try:
        sim.setObjectPosition(h, p, sim.handle_world)       # 4.6 and newer
    except Exception:
        sim.setObjectPosition(h, sim.handle_world, p)       # older versions


def colour(h, rgb):
    sim.setShapeColor(h, None, sim.colorcomponent_ambient_diffuse, rgb)


def prep(h, rgb):
    sim.setObjectInt32Param(h, sim.shapeintparam_static, 1)
    colour(h, rgb)


def shape(size, pos, rgb):
    h = sim.createPrimitiveShape(sim.primitiveshape_cuboid, size, 0)
    prep(h, rgb)
    setpos(h, pos)
    return h


pose = yaml.safe_load(open(os.path.join(CAD, "derived.yaml")))["nominal_pose_bbox_centre_m"]


def part(name, rgb):
    h = sim.importShape(4, os.path.join(CAD, f"{name}_local.stl"), 0, 0.0001, 1.0)
    prep(h, rgb)
    p = pose[name]
    setpos(h, [p["x"], p["y"], p["z"]])
    return h, p


buoy, bp = part("buoy", [0.95, 0.95, 0.95])
part("shaft", [0.35, 0.35, 0.4])
part("base", [0.3, 0.3, 0.3])
shape([2.4, 0.8, 0.02], [0, 0, -1.0], [0.76, 0.68, 0.45])                 # seabed
pto = shape([0.08, 0.08, 0.06], [0, 0, 0.23], [0.0, 1.0, 0.0])            # colour = damping

xs = np.linspace(-1.0, 1.0, NX)
dx = xs[1] - xs[0]
water = [(x, shape([dx, 0.6, 0.5], [x, 0, -0.25], [0.15, 0.5, 0.9])) for x in xs if abs(x) > 0.14]
for _, h in water:
    sim.setShapeColor(h, None, sim.colorcomponent_transparency, [0.55])

seg = [shape([0.06, 0.06, 0.04], [0, -0.55, i * 0.055], [0.2, 0.2, 0.2]) for i in range(8)]
lit = -1


def dcol(b):                                  # green (low) to red (high)
    f = float(np.clip((b - 3.0) / 37.0, 0, 1))
    return [f, 1.0 - f, 0.0]


def meter(b):
    global lit
    n = int(np.clip(round(8 * (b - 3.0) / 37.0), 0, 8))
    if n != lit:
        for i in range(8):
            colour(seg[i], dcol(b) if i < n else [0.2, 0.2, 0.2])
        lit = n


def draw(a_i, t, z, b):
    setpos(buoy, [bp["x"], bp["y"], bp["z"] + VIS * z])
    colour(pto, dcol(b))
    meter(b)
    for x, h in water:
        eta = float(np.sum(a_i * np.cos(w * t - kw * x + ph)))
        setpos(h, [x, 0, VIS * eta - 0.25])


# ---------------------------------------------------------------- run
md = read_model_description(args.fmu)
vr = {v.name: v.valueReference for v in md.modelVariables}
ud = extract(args.fmu)
baseline = args.fixed is not None
ser = None if baseline else serial.serial_for_url(args.port, baudrate=115200, timeout=0.5)
mode = f"fixed {args.fixed:g}" if baseline else "ESP32 controller"

sim.startSimulation()
try:
    for Tp, Hs, tsea_cur in SEAS:
        print(f"--- {mode}: sea Tp={Tp} s, Hs={Hs} m, {tsea_cur:.0f} s, speed x{args.speed} ---")
        if not (w0 <= 2 * np.pi / Tp <= w1):
            print(f"    warning: Tp={Tp} s is outside the wave-component band "
                  f"({2*np.pi/w1:.2f}-{2*np.pi/w0:.2f} s) set by pto_table.csv")
        if ser:
            ser.reset_input_buffer()
            ser.write(b"r\n")
            while ser.readline().decode(errors="ignore").strip() != "ok":
                pass

        a_i = amps(Hs, Tp)
        fmu = FMU2Slave(guid=md.guid, unzipDirectory=ud,
                        modelIdentifier=md.coSimulation.modelIdentifier, instanceName="x")
        fmu.instantiate()
        fmu.setupExperiment(startTime=0.0)
        fmu.setReal([vr[f"ai[{i+1}]"] for i in range(N)], [float(v) for v in a_i])
        fmu.setReal([vr["Bpto"]], [B_START])
        fmu.enterInitializationMode()
        fmu.exitInitializationMode()

        with lock:
            S.v = S.V = S.I = S.P = S.Vrms = S.Irms = S.Pavg = S.E = 0.0
            S.t, S.z, S.b, S.b_ver, S.tp, S.done = 0.0, 0.0, B_START, 0, 0.0, False
        stop.clear()
        log = []
        tp_thread = threading.Thread(target=plant_thread, args=(fmu, vr), daemon=True)
        te_thread = threading.Thread(target=esp_thread, args=(ser, log), daemon=True) if ser else None
        tp_thread.start()
        if te_thread:
            te_thread.start()

        frame = 1.0 / FPS
        nxt = time.perf_counter()
        last_print = -1.0
        while True:
            with lock:
                t, z, b, tp, done = S.t, S.z, S.b, S.tp, S.done
                Vr, Ir, Pa = S.Vrms, S.Irms, S.Pavg
            if done:
                break
            draw(a_i, t, z, b)
            if t - last_print >= args.speed:
                print(f"t={t:5.1f}  Tp_est={tp:.2f}  Bpto={b:5.1f}  "
                      f"V={Vr:5.2f} Vrms  I={Ir:5.3f} Arms  P={Pa*1000:7.2f} mW (5 s avg)")
                last_print = t
            nxt += frame
            sleep = nxt - time.perf_counter()
            if sleep > 0:
                time.sleep(sleep)
            else:
                nxt = time.perf_counter()     # fell behind, don't try to catch up

        stop.set()
        tp_thread.join()
        if te_thread:
            te_thread.join(timeout=3)
        fmu.terminate()
        fmu.freeInstance()
        with lock:
            E_tot, t_end = S.E, S.t
        mean_mw = E_tot / max(t_end, 1e-9) * 1000
        print(f"    mean electrical power over run: {mean_mw:.2f} mW "
              f"({E_tot:.3f} J in {t_end:.1f} s)")
        row = dict(mode=mode, speed=args.speed, Tp=Tp, Hs=Hs,
                   duration_s=round(t_end, 1), mean_mW=round(mean_mw, 2))
        pd.DataFrame([row]).to_csv("live_results.csv", mode="a", index=False,
                                   header=not os.path.exists("live_results.csv"))
        if ser:
            df = pd.DataFrame(log, columns=["t", "z", "bpto", "tp_est", "cmd_latency_s", "wall_s",
                                            "v_rms", "i_rms", "p_avg_w"])
            df.to_csv(f"esp_rt_log_Tp{Tp}.csv", index=False)
            if len(df):
                rate = len(df) / df.wall_s.iloc[-1]
                print(f"    {len(df)} samples answered (~{rate:.1f}/s wall), "
                      f"mean command latency {df.cmd_latency_s.mean():.3f} s plant time")
finally:
    stop.set()
    sim.stopSimulation()