"""
Real-time WEC demo: one plant clock, ESP32 as a monitor, CoppeliaSim as a viewer.

  plant thread : WEC_Heave_Wide FMU (Model Exchange) integrated with CVode, paced to the
                 wall clock (x --speed). Applies the latest Bpto it has received. Single
                 source of truth for waves, buoy, power.
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

Solver (--solver):
  cvode : Model Exchange interface of the FMU, integrated by fmpy's CVode. fmpy calls a hook every
          OUT_DT plant seconds; pacing, Bpto slew and the shared-state update happen there.
  cs    : co-simulation interface, doStep() loop with the FMU's own fixed-step solver (what the
          previous version of this file did). The exported DLL has no SUNDIALS, so no CVode here.
  auto  : (default) try cvode, and if it fails before the first sample fall back to cs, loudly.

Sea state: Python computes the wave amplitudes (JONSWAP, scaled to Hs, same formula as the .mo)
on the FMU's own frequency grid wi[], 2-16 rad/s, so Tp 0.5-2.0 s all fit. The FMU's ai[] is read
back and compared as a sanity check.

Keep the lockstep loop for your evaluation numbers; this one has timing jitter.
"""
import argparse
import os
import threading
import time
import traceback
from collections import deque
from dataclasses import dataclass

import numpy as np
import pandas as pd
import serial
import yaml
import fmpy
from fmpy import read_model_description, simulate_fmu, extract
from fmpy.fmi2 import FMU2Slave
from coppeliasim_zmqremoteapi_client import RemoteAPIClient

# ---------------------------------------------------------------- config
FMU = r"G:\Projects\Wave_Energy_Converter\WEC_Heave_Wide.fmu"
CAD = r"G:\Projects\Wave_Energy_Converter\wec-twin\cad_out"
# (Tp s, Hs m, duration s), in order: low -> resonance (natural period ~0.81 s) -> bigger ->
# high -> higher -> even higher. Hs stays 0.04 m, the basis of pto_table_wide.csv. Longer periods
# get longer runs so the controller sees ~15-20 waves. Table / FMU range is Tp 0.5 .. 2.0 s.
SEAS = [(0.6, 0.04, 12.0), (0.81, 0.04, 12.0), (0.95, 0.04, 20.0),
        (1.2, 0.04, 20.0), (1.6, 0.04, 30.0), (2.0, 0.04, 40.0)]
TP_MIN, TP_MAX = 0.5, 2.0                             # range of pto_table_wide.csv
OUT_DT = 0.005                                        # communication step (plant time)
H = 0.0005                                            # doStep size for --solver cs
RTOL = 1e-6                                           # CVode relative tolerance
VIS = 5.0                                             # visual exaggeration of motion
NX = 41                                               # water columns (Tp 0.5 s waves are only 0.4 m long; use 25 if CoppeliaSim lags)
FPS = 25.0
G, DEPTH = 9.81, 1.0
B_INIT = 7.0     # fallback only; every sea now starts at its own pto_table_wide.csv value (see pto_start)
_HERE = os.path.dirname(os.path.abspath(__file__))
_PTO = None


def _load_pto():
    """Find pto_table_wide.csv: --pto, else next to the script, in results/ or wec-twin/results/ (script folder, then cwd)."""
    cands = ([args.pto] if args.pto else []) + [
        os.path.join(d, sub, "pto_table_wide.csv")
        for d in (_HERE, os.getcwd()) for sub in ("", "results", os.path.join("wec-twin", "results"))]
    for c in cands:
        if os.path.isfile(c):
            print(f"pto table: {c}")
            return pd.read_csv(c).sort_values("Tp")
    raise FileNotFoundError("pto_table_wide.csv not found; pass --pto <path>. Looked in:\n  " + "\n  ".join(cands))


def pto_start(tp):
    """Bpto the table gives for this sea's Tp (linear interpolation, same as pto_lookup on the board)."""
    global _PTO
    if _PTO is None:
        _PTO = _load_pto()
    return float(np.interp(tp, _PTO["Tp"], _PTO["best_Bpto"]))

B_LO, B_HI = 3.0, 100.0                               # LED / colour range, matches main.cpp B_MIN / B_MAX
BPTO_SLEW = 10.0  # max change of the applied damping, per plant second (was 4 for the 3..40 range)

# Rough generator model (adjust to your real generator data in params.yaml)
KF = 20.0        # force constant, N/A  (= EMF constant in V per m/s, ideal lossless coupling)
ETA = 1.0        # efficiency applied to the power; 1.0 = the FMU's own P = Bpto*v^2 (same basis as pto_table_wide.csv)
P_WIN = 5.0      # seconds of plant time for the moving-average power

ap = argparse.ArgumentParser()
ap.add_argument("--fmu", default=FMU)
ap.add_argument("--speed", type=float, default=0.5,
                help="plant seconds per wall second (0.25 = slow motion, gives Wokwi more samples per wave)")
ap.add_argument("--port", default="rfc2217://localhost:4000")
ap.add_argument("--pto", default=None, help="path to pto_table_wide.csv (default: next to this script or in results/)")
ap.add_argument("--solver", choices=["auto", "cvode", "cs"], default="auto")
ap.add_argument("--fixed", type=float, default=None,
                help="baseline: hold Bpto at this value, no ESP32 needed")
args = ap.parse_args()
B_START = args.fixed if args.fixed is not None else B_INIT

md = read_model_description(args.fmu)
vr = {v.name: v.valueReference for v in md.modelVariables}
N = sum(1 for n in vr if n.startswith("wi["))         # wave components (40 in WEC_Heave_Wide)


GAM = 3.3


def amps(w, dw, Hs, Tp):
    """JONSWAP amplitudes on the FMU's frequency grid, scaled so that sum(a^2)/2 = (Hs/4)^2."""
    wp = 2 * np.pi / Tp
    sig = np.where(w <= wp, 0.07, 0.09)
    Sw = w**-5 * np.exp(-1.25 * (wp / w)**4) * GAM**np.exp(-(w - wp)**2 / (2 * sig**2 * wp**2))
    a = np.sqrt(2 * Sw * dw)
    return a * (Hs / 4) / np.sqrt(np.sum(a**2) / 2)


def dispersion(w):
    k = w**2 / G
    for _ in range(50):
        k = w**2 / (G * np.tanh(k * DEPTH))
    return k


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
WAVE = None                 # (w, a, ph, kw) read from the FMU after initialisation, same for plant and visuals
lock = threading.Lock()
stop = threading.Event()


# ---------------------------------------------------------------- threads
class Sampler:
    """Everything that happens once per communication step, whichever solver produced the state."""

    def __init__(self, Tp, Hs):
        self.Tp, self.Hs = Tp, Hs
        self.first = True
        self.started = False
        self.t0 = 0.0
        self.b = B_START
        self.t_prev = 0.0
        self.EV = self.EI = 0.0
        self.next_snap = 0.1
        self.hist = deque([(0.0, 0.0, 0.0, 0.0)])     # (t, E, EV, EI) for the moving window

    def sample(self, t, get, put):
        """get(name) -> float, put(name, value). Returns False to stop the run."""
        global WAVE
        if stop.is_set():
            return False
        if self.first:
            self.first = False
            wi = np.array([get(f"wi[{i + 1}]") for i in range(N)])
            ph = np.array([get(f"ph[{i + 1}]") for i in range(N)])
            ai_fmu = np.array([get(f"ai[{i + 1}]") for i in range(N)])
            dw = (get("wmax") - get("wmin")) / N
            a = amps(wi, dw, self.Hs, self.Tp)                       # sea state, Python side
            dev = float(np.max(np.abs(a - ai_fmu)) / np.max(ai_fmu))
            print(f"    sea state: {N} components {wi.min():.1f}-{wi.max():.1f} rad/s, "
                  f"max deviation from the FMU's ai[] {dev * 100:.2g} %"
                  + ("   <-- check Hs/Tp/gam" if dev > 0.02 else ""))
            with lock:
                WAVE = (wi, a, ph, dispersion(wi))
            self.t0 = time.perf_counter() - t / args.speed
        self.started = True

        # pace to the wall clock; if we are behind, just carry on
        lag = self.t0 + t / args.speed - time.perf_counter()
        if lag > 0:
            time.sleep(lag)

        z, v, E = get("z"), get("v"), get("E") * ETA
        b = self.b                                 # damping that was applied over the last step
        dtb = t - self.t_prev
        if dtb > 0:
            I = b * v / KF                         # A      (F = b*v, F/KF)
            V = KF * v                             # V      (V*I = F*v)
            P = ETA * b * v * v                    # W electrical
            self.EV += V * V * dtb
            self.EI += I * I * dtb
        else:
            V = I = P = 0.0
        self.t_prev = t

        if t >= self.next_snap:
            self.hist.append((t, E, self.EV, self.EI))
            self.next_snap = t + 0.1
            while len(self.hist) > 2 and t - self.hist[1][0] >= P_WIN:
                self.hist.popleft()
        ta, Ea, EVa, EIa = self.hist[0]
        Vrms = Irms = Pavg = 0.0
        if t - ta > 0:
            Pavg = (E - Ea) / (t - ta)
            Vrms = ((self.EV - EVa) / (t - ta)) ** 0.5
            Irms = ((self.EI - EIa) / (t - ta)) ** 0.5

        with lock:
            b_target = S.b
            S.t, S.z, S.v, S.V, S.I, S.P = t, z, v, V, I, P
            S.Vrms, S.Irms, S.Pavg, S.E = Vrms, Irms, Pavg, E

        # ramp the applied damping toward the ESP32 command instead of jumping
        max_move = BPTO_SLEW * max(dtb, 0.0)
        new_b = b + max(-max_move, min(max_move, b_target - b))
        if abs(new_b - b) > 1e-4:
            self.b = new_b
            put("Bpto", new_b)
        return True


def run_cvode(sm, Tp, Hs, tsea):
    """Model Exchange + CVode, driven by fmpy; fmpy calls the hook after every output step."""
    def hook(t, recorder):
        fmu = getattr(recorder, "fmu", None)
        if fmu is None:
            raise RuntimeError("this fmpy version does not expose recorder.fmu to step_finished")
        return sm.sample(t, lambda n: fmu.getReal([vr[n]])[0], lambda n, x: fmu.setReal([vr[n]], [x]))

    simulate_fmu(args.fmu, start_time=0.0, stop_time=tsea, fmi_type="ModelExchange", solver="CVode",
                 relative_tolerance=RTOL, output_interval=OUT_DT,
                 start_values={"Tp": Tp, "Hs": Hs, "Bpto": B_START},
                 output=["z"], record_events=False, validate=False, step_finished=hook)


def run_cs(sm, Tp, Hs, tsea):
    """Co-simulation FMU, doStep loop (fixed internal step), same as the previous version of this file."""
    ud = extract(args.fmu)
    fmu = FMU2Slave(guid=md.guid, unzipDirectory=ud,
                    modelIdentifier=md.coSimulation.modelIdentifier, instanceName="x")
    fmu.instantiate()
    fmu.setupExperiment(startTime=0.0)
    fmu.setReal([vr["Tp"], vr["Hs"], vr["Bpto"]], [Tp, Hs, B_START])
    fmu.enterInitializationMode()
    fmu.exitInitializationMode()
    get = lambda n: fmu.getReal([vr[n]])[0]
    put = lambda n, x: fmu.setReal([vr[n]], [x])
    try:
        t = 0.0
        sm.sample(t, get, put)
        while t < tsea - 1e-9:
            t_next = min(t + OUT_DT, tsea)
            while t < t_next - 1e-9:
                fmu.doStep(t, H)
                t += H
            if not sm.sample(t, get, put):
                break
    finally:
        fmu.terminate()
        fmu.freeInstance()


def plant_thread(Tp, Hs, tsea):
    sm = Sampler(Tp, Hs)
    try:
        if args.solver in ("auto", "cvode"):
            try:
                run_cvode(sm, Tp, Hs, tsea)
            except Exception:
                traceback.print_exc()
                if args.solver == "cvode" or sm.started:
                    print("!!! plant (CVode) failed, see traceback above")
                else:
                    print("!!! CVode path failed before the first sample, falling back to --solver cs "
                          "(co-simulation, FMU's own fixed-step solver)")
                    run_cs(Sampler(Tp, Hs), Tp, Hs, tsea)
        else:
            run_cs(sm, Tp, Hs, tsea)
    except Exception:
        traceback.print_exc()
        print("!!! plant failed, see traceback above")
    finally:
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
    f = float(np.clip((b - B_LO) / (B_HI - B_LO), 0, 1))
    return [f, 1.0 - f, 0.0]


def meter(b):
    global lit
    n = int(np.clip(round(8 * (b - B_LO) / (B_HI - B_LO)), 0, 8))
    if n != lit:
        for i in range(8):
            colour(seg[i], dcol(b) if i < n else [0.2, 0.2, 0.2])
        lit = n


def draw(wave, t, z, b):
    w, a_i, ph, kw = wave
    setpos(buoy, [bp["x"], bp["y"], bp["z"] + VIS * z])
    colour(pto, dcol(b))
    meter(b)
    for x, h in water:
        eta = float(np.sum(a_i * np.cos(w * t - kw * x + ph)))
        setpos(h, [x, 0, VIS * eta - 0.25])


# ---------------------------------------------------------------- run
baseline = args.fixed is not None
ser = None if baseline else serial.serial_for_url(args.port, baudrate=115200, timeout=0.5)
mode = f"fixed {args.fixed:g}" if baseline else "ESP32 controller"
print(f"fmpy {fmpy.__version__}, WEC_Heave_Wide, --solver {args.solver} "
      f"(cvode: Model Exchange, rtol {RTOL:g}; cs: doStep H={H * 1000:g} ms), output step {OUT_DT * 1000:g} ms")

sim.startSimulation()
try:
    for Tp, Hs, tsea_cur in SEAS:
        print(f"--- {mode}: sea Tp={Tp} s, Hs={Hs} m, {tsea_cur:.0f} s, speed x{args.speed} ---")
        if not (TP_MIN <= Tp <= TP_MAX):
            print(f"    warning: Tp={Tp} s is outside pto_table_wide.csv ({TP_MIN}-{TP_MAX} s)")
        B_START = args.fixed if args.fixed is not None else pto_start(Tp)   # start the sea at its table value
        print(f"    start Bpto = {B_START:.1f}")
        if ser:
            ser.reset_input_buffer()
            ser.write(f"r,{Tp:.3f}\n".encode())   # board starts from pto_lookup(Tp)
            while ser.readline().decode(errors="ignore").strip() != "ok":
                pass

        with lock:
            WAVE = None
            S.v = S.V = S.I = S.P = S.Vrms = S.Irms = S.Pavg = S.E = 0.0
            S.t, S.z, S.b, S.b_ver, S.tp, S.done = 0.0, 0.0, B_START, 0, 0.0, False
        stop.clear()
        log = []
        tp_thread = threading.Thread(target=plant_thread, args=(Tp, Hs, tsea_cur), daemon=True)
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
                wave = WAVE
            if done:
                break
            if wave is not None:                 # FMU not initialised yet until the first step
                draw(wave, t, z, b)
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