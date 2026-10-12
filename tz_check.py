import numpy as np, csv, os
from fmpy import simulate_fmu

FMU  = "WEC_Heave_Wide.fmu"
HS   = 0.04
RTOL = 1e-6
T_STOP, T_SKIP, DT = 300.0, 20.0, 0.01     # longer run = more waves for statistics
HYS, MEAN_TAU = 0.0003, 10.0               # same as the firmware
T_MIN, T_MAX = 0.30, 3.00


def crossing_periods(t, z):
    # firmware-style detector: EMA mean removal, hysteresis, interpolated upward crossings
    zmean, armed = z[0], False
    prev_zc, prev_t, last = 0.0, t[0], None
    out = []
    for k in range(1, len(t)):
        dt = t[k] - prev_t
        zmean += dt / (MEAN_TAU + dt) * (z[k] - zmean)
        zc = z[k] - zmean
        if zc < -HYS:
            armed = True
        if armed and prev_zc < 0.0 <= zc:
            tc = prev_t + (0.0 - prev_zc) / (zc - prev_zc) * dt
            if last is not None and t[k] > T_SKIP:
                T = tc - last
                if T_MIN <= T <= T_MAX:
                    out.append(T)
            last = tc
            armed = False
        prev_t, prev_zc = t[k], zc
    return np.array(out)


if __name__ == "__main__":
    tp_list, b_opt = [], []
    with open("results/pto_table.csv") as f:
        for row in csv.DictReader(f):
            if row["valid"] == "True":
                tp_list.append(float(row["Tp_s"]))
                b_opt.append(float(row["Bpto_best_Ns_per_m"]))

    print(f"{'Tp':>5} {'B*':>7} {'Tz_mean':>8} {'Tz_std':>7} {'Tp/Tz':>7} {'waves':>6}")
    rows = []
    for tp, b in zip(tp_list, b_opt):
        r = simulate_fmu(
            FMU, fmi_type="ModelExchange", solver="CVode", relative_tolerance=RTOL,
            start_values={"Tp": tp, "Hs": HS, "Bpto": b},
            stop_time=T_STOP, output_interval=DT, output=["z"],
        )
        T = crossing_periods(r["time"], r["z"])
        tz, sd = T.mean(), T.std()
        rows.append((tp, b, tz, sd, tp / tz, len(T)))
        print(f"{tp:5.2f} {b:7.2f} {tz:8.3f} {sd:7.3f} {tp/tz:7.3f} {len(T):6d}", flush=True)

    os.makedirs("results", exist_ok=True)
    with open("results/tz_check.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Tp_s", "Bpto", "Tz_mean_s", "Tz_std_s", "Tp_over_Tz", "n_waves"])
        w.writerows(rows)

    ratios = np.array([r[4] for r in rows])
    print(f"\nTp/Tz: min {ratios.min():.3f}, max {ratios.max():.3f}, mean {ratios.mean():.3f}")
    print("Firmware uses TZ_TO_TP = 1.03")