import numpy as np, pandas as pd
from fmpy import simulate_fmu

FMU = r"G:\Projects\Wave_Energy_Converter\WEC_Heave.fmu"
N, gam, T, T0 = 20, 3.3, 60, 10
Tps, Hss = [round(x, 2) for x in np.arange(0.70, 1.21, 0.05)], [0.04]
Bs = [3, 5, 7, 9, 12, 16, 20, 25, 30, 40, 50, 65]

d = pd.read_csv("results/mat_sweep.csv")
w0, w1 = d.omega_rad_s.min(), d.omega_rad_s.max()
dw = (w1 - w0) / N
w = w0 + dw * (np.arange(N) + 0.5)

def amps(Hs, Tp):
    wp = 2 * np.pi / Tp
    sig = np.where(w <= wp, 0.07, 0.09)
    S = w**-5 * np.exp(-1.25 * (wp / w)**4) * gam**np.exp(-(w - wp)**2 / (2 * sig**2 * wp**2))
    a = np.sqrt(2 * S * dw)
    return a * (Hs / 4) / np.sqrt(np.sum(a**2) / 2)

rows = []
for Tp in Tps:
    for Hs in Hss:
        sv = {f"ai[{i+1}]": float(v) for i, v in enumerate(amps(Hs, Tp))}
        P = []
        for b in Bs:
            r = simulate_fmu(FMU, stop_time=T, start_values={**sv, "Bpto": b},
                             output=["E"], output_interval=0.0005)
            t, E = r["time"], r["E"]
            i0 = np.searchsorted(t, T0)
            P.append((E[-1] - E[i0]) / (T - T0))
        P = np.array(P)
        best = Bs[int(P.argmax())]
        rows.append(dict(Tp=Tp, Hs=Hs, best_Bpto=best, P_best_W=P.max(),
                         **{f"P_B{b}": p for b, p in zip(Bs, P)}))
        print(f"Tp={Tp} Hs={Hs}: best Bpto={best}, P={P.max():.4f} W")

pd.DataFrame(rows).to_csv("results/pto_table.csv", index=False)
print("wrote results/pto_table.csv")
df = pd.DataFrame(rows)
edge = df[(df.best_Bpto == Bs[0]) | (df.best_Bpto == Bs[-1])]
if len(edge):
    print("WARNING: best Bpto at the edge of the search range for Tp =", list(edge.Tp))
print("Tp  =", list(df.Tp))
print("Bpto=", list(df.best_Bpto))