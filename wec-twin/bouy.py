#!/usr/bin/env python3
"""stage1_buoy.py - Stage 1a: heave hydrodynamics of the buoy alone (Capytaine 3.x).

Reads params.yaml (and cad_out/derived.yaml if present), builds the buoy from
Capytaine's own cylinder mesh (cleaner panels than the STL tessellation), runs a
frequency sweep for heave, and writes results/buoy_only.csv + results/buoy_only.png.

Usage (project folder, venv active):
    python stage1_buoy.py
    python stage1_buoy.py --faces-max-radius 0.012 --omega-min 2 --omega-max 12 --n 21

Self-checks printed at the start (read them before trusting the sweep):
  - displaced volume of the immersed mesh vs pi r^2 draft
  - heave hydrostatic stiffness vs rho g pi r^2
If either is off, the mesh is wrong (missing caps or inward normals) and results are junk.

Quantities (regular waves, per metre of wave amplitude):
  A33, B33    heave added mass (kg) and radiation damping (N s/m)
  |Fe|        excitation force amplitude (N per m of wave amplitude)
  B_opt       optimal passive PTO damping = |intrinsic impedance|
              = sqrt(B33^2 + (w (m + A33) - K / w)^2)
  P_max       maximum absorbable power = |Fe|^2 / (8 B33), W per m^2 of amplitude^2
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np
import xarray as xr
import yaml
import capytaine as cpt
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent   # so the script works from any working directory


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def as_vector(da, n, what):
    """Select-and-squeeze helper: return a 1-D numpy array of length n."""
    arr = np.squeeze(np.asarray(da.values))
    if arr.shape != (n,):
        die(f"{what}: expected shape ({n},) after selection, got {arr.shape}; "
            f"remaining dims: {da.dims}")
    return arr


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=str(HERE / "params.yaml"))
    ap.add_argument("--faces-max-radius", type=float, default=0.015,
                    help="max panel radius in m (smaller = finer mesh)")
    ap.add_argument("--omega-min", type=float, default=2.0)
    ap.add_argument("--omega-max", type=float, default=12.0,
                    help="keep below ~14 rad/s: first irregular frequency of this hull")
    ap.add_argument("--n", type=int, default=21, help="number of frequencies")
    ap.add_argument("--out", default=str(HERE / "results"))
    ap.add_argument("--show", action="store_true", help="open the plot window")
    args = ap.parse_args()

    with open(args.params, "r", encoding="utf-8") as f:
        p = yaml.safe_load(f)
    env, b = p["environment"], p["buoy"]
    r, h, d = b["radius"], b["height"], b["draft"]
    rho, g, depth = env["rho_water"], env["gravity"], env["water_depth"]

    derived_path = HERE / p["export"]["out_dir"] / "derived.yaml"
    if derived_path.exists():
        with open(derived_path, "r", encoding="utf-8") as f:
            der = yaml.safe_load(f)["buoy"]
        mass, cog_z = der["mass_kg"], der["cog_z_m"]
    else:
        print("note: derived.yaml not found; computing mass and CoG from params.yaml")
        mass = b["mass"] if b["mass"] is not None else rho * math.pi * r * r * d
        cog_z = -d + h / 2.0

    # ---- body: full cylinder, then clip at the free surface (z = 0)
    z_center = -d + h / 2.0
    try:
        mesh = cpt.mesh_vertical_cylinder(length=h, radius=r, center=(0, 0, z_center),
                                          faces_max_radius=args.faces_max_radius,
                                          name="buoy")
    except TypeError as e:
        die(f"mesh_vertical_cylinder rejected the arguments ({e}). "
            "Run help(cpt.mesh_vertical_cylinder) and send me the signature.")

    body = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, cog_z)),
        center_of_mass=(0, 0, cog_z),
        name="buoy",
    )
    body = body.immersed_part()   # 3.x: returns a new body (no in-place keep_immersed_part)

    # ---- self-checks
    try:
        hs = body.compute_hydrostatics(rho=rho, g=g)
    except TypeError:
        hs = body.compute_hydrostatics(rho=rho)
    vol_expected = math.pi * r * r * d
    vol_got = float(np.squeeze(np.asarray(hs["disp_volume"])))
    k_expected = rho * g * math.pi * r * r
    print(f"displaced volume : {vol_got:.6f} m^3   expected {vol_expected:.6f} "
          f"({abs(vol_got - vol_expected) / vol_expected * 100:.2f}% off)")
    try:
        K_da = hs["hydrostatic_stiffness"]
        k_got = float(K_da.sel(influenced_dof="Heave", radiating_dof="Heave"))
    except Exception:
        try:
            k_got = float(np.asarray(hs["hydrostatic_stiffness"].values)[2, 2])
        except Exception:
            k_got = float("nan")
    print(f"heave stiffness  : {k_got:.2f} N/m   expected {k_expected:.2f}")
    bad = (abs(vol_got - vol_expected) / vol_expected > 0.03) or \
          (not math.isnan(k_got) and abs(k_got - k_expected) / k_expected > 0.03)
    if bad:
        print("WARNING: hydrostatics disagree with the analytic cylinder by >3%. "
              "Check the mesh (caps, normal direction, resolution) before trusting results.")

    # ---- BEM sweep (heave only, wave direction 0; the buoy is axisymmetric)
    omega = np.linspace(args.omega_min, args.omega_max, args.n)
    coords = {
        "omega": omega,
        "wave_direction": [0.0],
        "radiating_dof": ["Heave"],
        "water_depth": [depth],
    }
    if abs(rho - 1000.0) > 1e-9:
        coords["rho"] = [rho]
    if abs(g - 9.81) > 1e-9:
        coords["g"] = [g]
    test_matrix = xr.Dataset(coords=coords)

    print(f"running BEM: {args.n} frequencies, {args.omega_min}-{args.omega_max} rad/s ...")
    ds = cpt.BEMSolver().fill_dataset(test_matrix, body)
    print("dataset variables:", list(ds.data_vars))

    heave = dict(radiating_dof="Heave", influenced_dof="Heave")
    A33 = as_vector(ds["added_mass"].sel(**heave), args.n, "added_mass")
    B33 = as_vector(ds["radiation_damping"].sel(**heave), args.n, "radiation_damping")

    if "excitation_force" in ds:
        F = ds["excitation_force"]
    elif "diffraction_force" in ds and "Froude_Krylov_force" in ds:
        F = ds["diffraction_force"] + ds["Froude_Krylov_force"]
    else:
        die("no excitation force found in the dataset; see 'dataset variables' above")
    Fe = np.abs(as_vector(F.sel(wave_direction=0.0, influenced_dof="Heave"),
                          args.n, "excitation force"))

    # ---- derived quantities
    K = k_expected
    X = omega * (mass + A33) - K / omega               # reactance of the intrinsic impedance
    B_opt = np.sqrt(B33 ** 2 + X ** 2)                 # optimal passive PTO damping
    P_max = np.where(B33 > 1e-12, Fe ** 2 / (8.0 * np.maximum(B33, 1e-12)), np.nan)

    f = omega ** 2 * (mass + A33) - K                  # resonance where f crosses zero
    w_n = None
    idx = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
    if len(idx):
        i = idx[0]
        w_n = omega[i] - f[i] * (omega[i + 1] - omega[i]) / (f[i + 1] - f[i])

    # ---- outputs
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    table = np.column_stack([omega, 2 * np.pi / omega, A33, B33, Fe, B_opt, P_max])
    np.savetxt(out / "buoy_only.csv", table, delimiter=",",
               header="omega_rad_s,period_s,A33_kg,B33_Ns_per_m,Fe_N_per_m,B_opt_Ns_per_m,Pmax_W_per_m2",
               comments="")

    fig, ax = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True)
    ax[0, 0].plot(omega, A33, label="A33 (kg)")
    ax[0, 0].plot(omega, B33, label="B33 (N s/m)")
    ax[0, 0].set_title("Heave added mass and radiation damping")
    ax[0, 1].plot(omega, Fe)
    ax[0, 1].set_title("Excitation force amplitude (N per m of wave)")
    ax[1, 0].plot(omega, B_opt)
    ax[1, 0].set_title("Optimal passive PTO damping (N s/m)")
    ax[1, 1].semilogy(omega, P_max)
    ax[1, 1].set_title("Max absorbable power (W per m^2 of amplitude^2)")
    for a in ax.ravel():
        a.set_xlabel("omega (rad/s)")
        a.grid(alpha=0.3)
        if w_n is not None:
            a.axvline(w_n, color="k", ls=":", lw=1)
    ax[0, 0].legend()
    fig.savefig(out / "buoy_only.png", dpi=150)
    if args.show:
        plt.show()

    # ---- summary
    print()
    print(f"buoy mass                 : {mass:.3f} kg")
    print(f"A33 range                 : {A33.min():.3f} .. {A33.max():.3f} kg")
    if w_n is not None:
        print(f"heave resonance (with A33): {w_n:.2f} rad/s  (T = {2 * np.pi / w_n:.3f} s)")
        print(f"  without added mass      : {math.sqrt(K / mass):.2f} rad/s "
              f"(T = {2 * math.pi * math.sqrt(mass / K):.3f} s)")
    else:
        print("no heave resonance inside the sweep range; widen --omega-min/--omega-max")
    print(f"wrote {out / 'buoy_only.csv'} and {out / 'buoy_only.png'}")


if __name__ == "__main__":
    main()