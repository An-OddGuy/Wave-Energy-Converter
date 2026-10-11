import ctypes as C, numpy as np, pandas as pd
from fmpy import read_model_description, extract
from fmpy.fmi2 import FMU2Slave

FMU = r"G:\Projects\Wave_Energy_Converter\WEC_Heave.fmu"
T, T0, H, NSUB = 80.0, 20.0, 0.0005, 20      # controller every 20*H = 0.01 s
SEAS = [(0.75, 0.04), (1.0, 0.04), (1.2, 0.04)]
FIXED = [7, 12]

class Ctrl(C.Structure):
    _fields_ = [("t", C.c_float), ("dc", C.c_float), ("ms", C.c_float),
                ("prev", C.c_float), ("last_cross", C.c_float),
                ("armed", C.c_int), ("survival", C.c_int),
                ("tp_est", C.c_float), ("amp_rms", C.c_float), ("bpto", C.c_float),
                ("tp_est", C.c_float), ("amp_rms", C.c_float), ("bpto", C.c_float),
                ("z1", C.c_float), ("z2", C.c_float), ("fms", C.c_float)]
from pathlib import Path
HERE = Path(__file__).resolve().parent
lib = C.CDLL(str(HERE / "controller.dll"))
lib.ctrl_step.restype = C.c_float
lib.ctrl_step.argtypes = [C.POINTER(Ctrl), C.c_float, C.c_float]
lib.ctrl_tp.restype = C.c_float
lib.ctrl_tp.argtypes = [C.POINTER(Ctrl)]

d = pd.read_csv(HERE.parent / "wec-twin" / "results" / "mat_sweep.csv")
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

def run(Tp, Hs, bfix=None):
    fmu = FMU2Slave(guid=md.guid, unzipDirectory=ud,
                    modelIdentifier=md.coSimulation.modelIdentifier, instanceName="x")
    fmu.instantiate()
    fmu.setupExperiment(startTime=0.0)
    fmu.setReal([vr[f"ai[{i+1}]"] for i in range(N)], [float(v) for v in amps(Hs, Tp)])
    fmu.enterInitializationMode(); fmu.exitInitializationMode()
    c = Ctrl(); lib.ctrl_init(C.byref(c))
    b = float(bfix) if bfix else float(c.bpto)
    t, E0 = 0.0, 0.0
    for k in range(int(T / H)):
        if k % NSUB == 0:
            if bfix is None:
                z = fmu.getReal([vr["z"]])[0]
                b = float(lib.ctrl_step(C.byref(c), z, NSUB * H))
            fmu.setReal([vr["Bpto"]], [b])
        fmu.doStep(t, H); t += H
        if abs(t - T0) < H / 2:
            E0 = fmu.getReal([vr["E"]])[0]
    E1 = fmu.getReal([vr["E"]])[0]
    fmu.terminate(); fmu.freeInstance()
    return (E1 - E0) / (T - T0), lib.ctrl_tp(C.byref(c)), b

rows = []
for Tp, Hs in SEAS:
    res = {f"fixed {b}": run(Tp, Hs, b)[0] for b in FIXED}
    P, tp_est, b_end = run(Tp, Hs)
    res["controller"] = P
    print(f"Tp={Tp}: " + "  ".join(f"{k}={v:.4f} W" for k, v in res.items())
          + f"  | est Tp={tp_est:.2f}, final Bpto={b_end:.1f}")