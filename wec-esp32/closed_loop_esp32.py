import numpy as np, pandas as pd, serial
from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave

FMU = r"G:\Projects\Wave_Energy_Converter\WEC_Heave.fmu"
CSV = r"G:\Projects\Wave_Energy_Converter\wec-twin\results\mat_sweep.csv"
SEAS = [(0.75, 0.04), (1.0, 0.04), (1.2, 0.04)]
T, T0 = 40.0, 15.0            # run length and start of the averaging window (s)
H, NSUB = 0.0005, 40          # FMU step, one sample every 0.02 s

d = pd.read_csv(CSV)
N, gam = 20, 3.3
w0, w1 = d.omega_rad_s.min(), d.omega_rad_s.max()
dw = (w1 - w0) / N
w = w0 + dw * (np.arange(N) + 0.5)

def amps(Hs, Tp):
    wp = 2 * np.pi / Tp
    sig = np.where(w <= wp, 0.07, 0.09)
    S = w**-5 * np.exp(-1.25 * (wp / w)**4) * gam**np.exp(-(w - wp)**2 / (2 * sig**2 * wp**2))
    a = np.sqrt(2 * S * dw)
    return a * (Hs / 4) / np.sqrt(np.sum(a**2) / 2)

md = read_model_description(FMU)
vr = {v.name: v.valueReference for v in md.modelVariables}
ud = extract(FMU)
ser = serial.serial_for_url("rfc2217://localhost:4000", baudrate=115200, timeout=5)

def run_esp(Tp, Hs, bfix=None):
    if bfix is None:
        ser.reset_input_buffer()
        ser.write(b"r\n")
        while ser.readline().decode().strip() != "ok":
            pass                                    # skip boot-log lines
    fmu = FMU2Slave(guid=md.guid, unzipDirectory=ud,
                    modelIdentifier=md.coSimulation.modelIdentifier, instanceName="x")
    fmu.instantiate()
    fmu.setupExperiment(startTime=0.0)
    fmu.setReal([vr[f"ai[{i+1}]"] for i in range(N)], [float(v) for v in amps(Hs, Tp)])
    fmu.enterInitializationMode(); fmu.exitInitializationMode()
    b, tp, t, E0 = float(bfix or 7.0), 0.0, 0.0, 0.0
    log = []
    for k in range(int(T / H)):
        if k % NSUB == 0:                            # one sample every 0.02 s
            if bfix is None:
                z = fmu.getReal([vr["z"]])[0]
                ser.write(f"{z:.6f}\n".encode())
                while True:
                    line = ser.readline().decode().strip()
                    try:
                        b, tp = map(float, line.split(","))
                        break
                    except ValueError:
                        pass                         # ignore stray lines
                log.append((t, z, b, tp))
            fmu.setReal([vr["Bpto"]], [b])
        fmu.doStep(t, H); t += H
        if abs(t - T0) < H / 2:
            E0 = fmu.getReal([vr["E"]])[0]
    E1 = fmu.getReal([vr["E"]])[0]
    fmu.terminate(); fmu.freeInstance()
    if bfix is None:
        pd.DataFrame(log, columns=["t", "z", "bpto", "tp_est"]).to_csv(f"esp_log_Tp{Tp}.csv", index=False)
    return (E1 - E0) / (T - T0), tp, b

for Tp, Hs in SEAS:
    fx = run_esp(Tp, Hs, bfix=12)[0]
    P, tp, b = run_esp(Tp, Hs)
    print(f"Tp={Tp}: fixed 12={fx:.4f} W  ESP32={P:.4f} W | est Tp={tp:.2f}, final Bpto={b:.1f}")