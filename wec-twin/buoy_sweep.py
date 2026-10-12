"""Buoy-only Capytaine heave sweep over a wide band -> results/buoy_wide.csv
Run from wec-twin with the venv active:
    python buoy_sweep.py --omega-min 2.2 --omega-max 13.5 --n 40
Keep omega-max <= 13.5: above ~13.9 rad/s this buoy hits irregular frequencies (~15 s solve total).
"""
import argparse
import math
from pathlib import Path

import numpy as np
import yaml

from mat import build_buoy, solve_case          # reuse your own functions

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("--omega-min", type=float, default=2.2)
ap.add_argument("--omega-max", type=float, default=13.5)
ap.add_argument("--n", type=int, default=40)
ap.add_argument("--faces-max-radius", type=float, default=0.015)
a = ap.parse_args()

p = yaml.safe_load(open(HERE / "params.yaml", encoding="utf-8"))
env, b = p["environment"], p["buoy"]
der = yaml.safe_load(open(HERE / p["export"]["out_dir"] / "derived.yaml", encoding="utf-8"))
mass, cog_z = der["buoy"]["mass_kg"], der["buoy"]["cog_z_m"]
K = env["rho_water"] * env["gravity"] * math.pi * b["radius"] ** 2
omega = np.linspace(a.omega_min, a.omega_max, a.n)

buoy = build_buoy(b["radius"], b["height"], b["draft"], cog_z, a.faces_max_radius)
c = solve_case("no mat", [buoy], omega, env["water_depth"], env["rho_water"],
               env["gravity"], mass, K)

out = HERE / "results"
out.mkdir(exist_ok=True)
np.savetxt(out / "buoy_wide.csv",
           np.column_stack([omega, 2 * np.pi / omega, c["A33"], c["B33"], c["Fe"]]),
           delimiter=",", header="omega_rad_s,period_s,A33_no_mat,B33_no_mat,Fe_no_mat",
           comments="")
print(f"wrote {out / 'buoy_wide.csv'}  (mass {mass:.3f} kg, natural freq {c['w_n']:.2f} rad/s)")
