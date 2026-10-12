import numpy as np, csv, os
from fmpy import simulate_fmu
from scipy.optimize import minimize_scalar

FMU     = "WEC_Heave_Wide.fmu"
HS      = 0.04
TP_LIST = np.round(np.arange(0.5, 2.0 + 1e-9, 0.1), 3)
B_GRID  = np.geomspace(1.5, 300.0, 20)       # N s/m, coarse search
T_STOP, T_SKIP, DT = 120.0, 20.0, 0.01       # discard start-up transient
RTOL    = 1e-6

RHO, G = 1000.0, 9.81


def pmax(tp):
    # heave point-absorber upper bound (isotropic capture width = lambda/2pi)
    te  = 0.9 * tp
    lam = G * tp**2 / (2 * np.pi)
    return RHO * G**2 * HS**2 * te / (64 * np.pi) * lam / (2 * np.pi)


def mean_power(tp, bpto):
    r = simulate_fmu(
        FMU,
        fmi_type="ModelExchange",
        solver="CVode",
        relative_tolerance=RTOL,
        start_values={"Tp": float(tp), "Hs": HS, "Bpto": float(bpto)},
        stop_time=T_STOP,
        output_interval=DT,
        output=["E"],
    )
    i = np.searchsorted(r["time"], T_SKIP)
    return (r["E"][-1] - r["E"][i]) / (r["time"][-1] - r["time"][i])


def optimise(tp):
    coarse = np.array([mean_power(tp, b) for b in B_GRID])
    k  = int(np.argmax(coarse))
    lo = B_GRID[max(k - 1, 0)]
    hi = B_GRID[min(k + 1, len(B_GRID) - 1)]
    res = minimize_scalar(lambda b: -mean_power(tp, b),
                          bounds=(lo, hi), method="bounded",
                          options={"xatol": 0.05})
    b_best, p_best = float(res.x), float(-res.fun)

    at_edge = (b_best < B_GRID[0] * 1.05) or (b_best > B_GRID[-1] * 0.95)
    too_big = p_best > pmax(tp)
    valid   = (not at_edge) and (not too_big) and np.isfinite(p_best)
    return (float(tp), b_best, p_best, float(pmax(tp)), at_edge, too_big, valid)


if __name__ == "__main__":
    rows = []
    for tp in TP_LIST:
        row = optimise(tp)
        rows.append(row)
        tp_, b, p, pm, edge, big, ok = row
        flag = "" if ok else ("  <-- EDGE" if edge else "") + ("  <-- > Pmax" if big else "")
        print(f"Tp={tp_:.2f}s  Bpto*={b:7.2f}  Pmean={p:.4f} W  (Pmax={pm:.3f}){flag}", flush=True)

    os.makedirs("results", exist_ok=True)
    with open("results/pto_table.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Tp_s", "Bpto_best_Ns_per_m", "Pmean_W", "Pmax_W", "at_grid_edge", "above_Pmax", "valid"])
        w.writerows(rows)

    good = [r for r in rows if r[6]]
    if len(good) < 2:
        raise SystemExit("Fewer than 2 valid rows, not writing the header. Check results/pto_table.csv")

    tps  = ", ".join(f"{r[0]:.3f}f" for r in good)
    bpto = ", ".join(f"{r[1]:.3f}f" for r in good)
    with open("results/pto_table.h", "w") as f:
        f.write(f"""#pragma once
#define PTO_N {len(good)}
static const float PTO_TP[PTO_N]   = {{{tps}}};
static const float PTO_BPTO[PTO_N] = {{{bpto}}};

static inline float pto_lookup(float tp) {{
    if (tp <= PTO_TP[0]) return PTO_BPTO[0];
    if (tp >= PTO_TP[PTO_N-1]) return PTO_BPTO[PTO_N-1];
    int i = 0;
    while (PTO_TP[i+1] < tp) i++;
    float t = (tp - PTO_TP[i]) / (PTO_TP[i+1] - PTO_TP[i]);
    return PTO_BPTO[i] + t * (PTO_BPTO[i+1] - PTO_BPTO[i]);
}}
""")
    print(f"wrote results/pto_table.csv ({len(rows)} rows) and results/pto_table.h ({len(good)} valid rows)")