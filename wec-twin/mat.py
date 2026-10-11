#!/usr/bin/env python3
"""stage1b_mat.py - Stage 1b: does a fixed submerged ring mat change the buoy's heave hydrodynamics?

Runs the same heave sweep as stage1_buoy.py for:
  - baseline: buoy alone
  - buoy + mat, at each requested mat depth (top surface depth below SWL)
and compares excitation force, added mass, optimal PTO damping and maximum absorbable power.

The mat mesh comes from cad_out/mat_local.stl (your CadQuery geometry), moved to each
depth by translation (no CAD re-run needed), refined to a uniform panel size, and added as a
second body. The mat is fixed: only the buoy's heave is a radiating dof.

Usage (project folder, venv active):
    python stage1b_mat.py
    python stage1b_mat.py --depths 0.05 0.15 0.30 --n 15
    python stage1b_mat.py --mat-max-edge 0.03        # finer mat mesh (slower)

Notes:
  * BEM treats a thin plate poorly. If mat.thickness in params.yaml is below ~0.03 m, set it
    to 0.03, re-run make_cad.py and check_mesh.py, then run this script.
  * Regular-wave, linear potential-flow result. It is a go/no-go signal, not the final answer;
    the real comparison comes from irregular seas in Stage 3.
"""
import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np
import trimesh
import xarray as xr
import yaml
import capytaine as cpt
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent   # so the script works from any working directory


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def as_vector(da, n, what):
    arr = np.squeeze(np.asarray(da.values))
    if arr.shape != (n,):
        die(f"{what}: expected shape ({n},) after selection, got {arr.shape}; dims {da.dims}")
    return arr


def build_buoy(r, h, d, cog_z, faces_max_radius):
    z_center = -d + h / 2.0
    try:
        mesh = cpt.mesh_vertical_cylinder(length=h, radius=r, center=(0, 0, z_center),
                                          faces_max_radius=faces_max_radius, name="buoy")
    except TypeError as e:
        die(f"mesh_vertical_cylinder rejected the arguments ({e})")
    body = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=(0, 0, cog_z)),
        center_of_mass=(0, 0, cog_z),
        name="buoy",
    )
    return body.immersed_part()


def build_mat(cad_dir, pose, depth, max_edge, tmp_dir, m_par):
    """Structured ring mat (few panels), top surface at 'depth' below SWL."""
    ri, ro, t = m_par["ring_inner_radius"], m_par["ring_outer_radius"], m_par["thickness"]
    nr = max(2, math.ceil((ro - ri) / max_edge))
    na = max(12, math.ceil(2 * math.pi * ro / max_edge))
    zt = -depth
    zb = -depth - t
    z_center = (zt + zb) / 2.0

    radii = np.linspace(ri, ro, nr + 1)
    ang = np.linspace(0, 2 * np.pi, na, endpoint=False)
    verts = []
    for z in (zt, zb):
        for rr in radii:
            for a in ang:
                verts.append([pose["x"] + rr * np.cos(a), pose["y"] + rr * np.sin(a), z])
    verts = np.array(verts)

    def vid(layer, j, k):
        return layer * (nr + 1) * na + j * na + (k % na)

    faces = []
    for layer in (0, 1):                      # top and bottom faces
        for j in range(nr):
            for k in range(na):
                a, b = vid(layer, j, k), vid(layer, j, k + 1)
                c, d = vid(layer, j + 1, k + 1), vid(layer, j + 1, k)
                faces += [[a, b, c], [a, c, d]]
    for j in (0, nr):                         # inner and outer walls
        for k in range(na):
            a, b = vid(0, j, k), vid(0, j, k + 1)
            c, d = vid(1, j, k + 1), vid(1, j, k)
            faces += [[a, b, c], [a, c, d]]

    m = trimesh.Trimesh(vertices=verts, faces=np.array(faces), process=True)
    trimesh.repair.fix_normals(m)
    tmp_dir.mkdir(parents=True, exist_ok=True)
    path = tmp_dir / f"mat_depth_{depth:.3f}.stl"
    m.export(str(path))

    mesh = cpt.load_mesh(str(path), file_format="stl")
    body = cpt.FloatingBody(
        mesh=mesh,
        dofs=cpt.rigid_body_dofs(rotation_center=(pose["x"], pose["y"], z_center)),
        center_of_mass=(pose["x"], pose["y"], z_center),
        name="mat",
    )
    return body, t, len(m.faces)


def find_heave_dof(system, single):
    keys = list(system.dofs)
    if single:
        if "Heave" not in keys:
            die(f"no 'Heave' dof; dofs are {keys}")
        return "Heave"
    cand = [k for k in keys if k.endswith("Heave") and "buoy" in k.lower()]
    if not cand:
        die(f"could not find the buoy heave dof among {keys}")
    return cand[0]


def solve_case(label, bodies, omega, depth, rho, g, mass, K):
    single = len(bodies) == 1
    system = bodies[0] if single else cpt.Multibody(bodies)
    dof = find_heave_dof(system, single)

    coords = {"omega": omega, "wave_direction": [0.0],
              "radiating_dof": [dof], "water_depth": [depth]}
    if abs(rho - 1000.0) > 1e-9:
        coords["rho"] = [rho]
    if abs(g - 9.81) > 1e-9:
        coords["g"] = [g]
    test_matrix = xr.Dataset(coords=coords)

    solver = cpt.BEMSolver()
    t0 = time.time()
    try:
        ds = solver.fill_dataset(test_matrix, system, hydrostatics=False)
    except TypeError:
        ds = solver.fill_dataset(test_matrix, system)
    print(f"  [{label}] solved in {time.time() - t0:.0f} s")

    n = len(omega)
    sel = dict(radiating_dof=dof, influenced_dof=dof)
    A33 = as_vector(ds["added_mass"].sel(**sel), n, "added_mass")
    B33 = as_vector(ds["radiation_damping"].sel(**sel), n, "radiation_damping")
    if "excitation_force" in ds:
        F = ds["excitation_force"]
    elif "diffraction_force" in ds and "Froude_Krylov_force" in ds:
        F = ds["diffraction_force"] + ds["Froude_Krylov_force"]
    else:
        die(f"no excitation force in dataset; variables: {list(ds.data_vars)}")
    Fe = np.abs(as_vector(F.sel(wave_direction=0.0, influenced_dof=dof), n, "excitation force"))

    X = omega * (mass + A33) - K / omega
    B_opt = np.sqrt(B33 ** 2 + X ** 2)
    P_max = np.where(B33 > 1e-12, Fe ** 2 / (8.0 * np.maximum(B33, 1e-12)), np.nan)

    fz = omega ** 2 * (mass + A33) - K
    w_n = np.nan
    idx = np.where(np.sign(fz[:-1]) != np.sign(fz[1:]))[0]
    if len(idx):
        i = idx[0]
        w_n = omega[i] - fz[i] * (omega[i + 1] - omega[i]) / (fz[i + 1] - fz[i])
    return dict(label=label, A33=A33, B33=B33, Fe=Fe, B_opt=B_opt, P_max=P_max, w_n=w_n)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default=str(HERE / "params.yaml"))
    ap.add_argument("--depths", type=float, nargs="+", default=[0.05, 0.15, 0.30],
                    help="mat top-surface depths below SWL (m)")
    ap.add_argument("--faces-max-radius", type=float, default=0.015, help="buoy mesh")
    ap.add_argument("--mat-max-edge", type=float, default=0.04, help="mat panel edge (m)")
    ap.add_argument("--omega-min", type=float, default=2.0)
    ap.add_argument("--omega-max", type=float, default=12.0)
    ap.add_argument("--n", type=int, default=15, help="number of frequencies")
    ap.add_argument("--out", default=str(HERE / "results"))
    ap.add_argument("--show", action="store_true")
    args = ap.parse_args()

    with open(args.params, "r", encoding="utf-8") as f:
        p = yaml.safe_load(f)
    env, b, m_par = p["environment"], p["buoy"], p["mat"]
    r, h, d = b["radius"], b["height"], b["draft"]
    rho, g, depth_w = env["rho_water"], env["gravity"], env["water_depth"]

    cad_dir = HERE / p["export"]["out_dir"]
    derived_path = cad_dir / "derived.yaml"
    if not derived_path.exists():
        die(f"{derived_path} not found. Run make_cad.py first.")
    with open(derived_path, "r", encoding="utf-8") as f:
        der = yaml.safe_load(f)
    mass, cog_z = der["buoy"]["mass_kg"], der["buoy"]["cog_z_m"]
    pose = der["nominal_pose_bbox_centre_m"].get("mat")
    if pose is None:
        die("no mat in derived.yaml; set mat.enabled: true in params.yaml and re-run make_cad.py")

    K = rho * g * math.pi * r * r
    omega = np.linspace(args.omega_min, args.omega_max, args.n)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    # clearance check for requested depths
    seabed_top = -depth_w + p["base"]["thickness"]
    for dep in args.depths:
        if dep < 0 or -dep - m_par["thickness"] < seabed_top:
            die(f"mat depth {dep} m is outside the water column")

    print(f"baseline + {len(args.depths)} mat depths, {args.n} frequencies each")
    buoy = build_buoy(r, h, d, cog_z, args.faces_max_radius)
    cases = [solve_case("no mat", [buoy], omega, depth_w, rho, g, mass, K)]

    for dep in args.depths:
        mat, thickness, nfaces = build_mat(cad_dir, pose, dep, args.mat_max_edge,
                                           out / "tmp_meshes", m_par)
        if thickness < 0.029:
            print(f"  WARNING: mat thickness is {thickness * 1000:.0f} mm; thin plates are "
                  "poorly handled by BEM. Consider mat.thickness: 0.03 in params.yaml.")
        print(f"  mat depth {dep:.3f} m: mat mesh has {nfaces} panels")
        buoy = build_buoy(r, h, d, cog_z, args.faces_max_radius)   # fresh body each case
        cases.append(solve_case(f"mat {dep:.2f} m", [buoy, mat], omega, depth_w,
                                rho, g, mass, K))

    base = cases[0]

    # ---- CSV
    cols = [omega, 2 * np.pi / omega]
    head = ["omega_rad_s", "period_s"]
    for c in cases:
        tag = c["label"].replace(" ", "_")
        for key in ("A33", "B33", "Fe", "B_opt", "P_max"):
            cols.append(c[key])
            head.append(f"{key}_{tag}")
    np.savetxt(out / "mat_sweep.csv", np.column_stack(cols), delimiter=",",
               header=",".join(head), comments="")

    # ---- plot
    fig, ax = plt.subplots(2, 3, figsize=(14, 7.5), constrained_layout=True)
    panels = [("A33", "Added mass A33 (kg)", False), ("B33", "Radiation damping B33 (N s/m)", False),
              ("Fe", "Excitation force (N per m of wave)", False),
              ("B_opt", "Optimal PTO damping (N s/m)", False),
              ("P_max", "Max absorbable power (W per m^2 amp^2)", True)]
    for a, (key, title, logy) in zip(ax.ravel(), panels):
        for c in cases:
            (a.semilogy if logy else a.plot)(omega, c[key], label=c["label"])
        a.set_title(title)
    for c in cases[1:]:
        ax[1, 2].plot(omega, c["P_max"] / base["P_max"], label=c["label"])
    ax[1, 2].axhline(1.0, color="k", lw=0.8, ls=":")
    ax[1, 2].set_title("Pmax relative to no mat")
    for a in ax.ravel():
        a.set_xlabel("omega (rad/s)")
        a.grid(alpha=0.3)
    ax[0, 0].legend(fontsize=8)
    fig.savefig(out / "mat_sweep.png", dpi=150)
    if args.show:
        plt.show()

    # ---- summary
    i0 = int(np.nanargmin(np.abs(omega - base["w_n"]))) if not np.isnan(base["w_n"]) else 0
    print()
    print(f"{'case':<14}{'w_n rad/s':<11}{'A33 @w_n kg':<13}{'Pmax gain @base w_n':<21}"
          f"{'max gain':<10}{'mean gain (unweighted)'}")
    for c in cases:
        gain = c["P_max"] / base["P_max"]
        print(f"{c['label']:<14}{c['w_n']:<11.2f}{c['A33'][i0]:<13.3f}{gain[i0]:<21.3f}"
              f"{np.nanmax(gain):<10.3f}{np.nanmean(gain):.3f}")
    print()
    print("Reading it: gain > 1 means more absorbable power than the bare buoy at that frequency.")
    print("The mean is unweighted by any wave spectrum, so treat it as a first signal only.")
    print(f"wrote {out / 'mat_sweep.csv'} and {out / 'mat_sweep.png'}")


if __name__ == "__main__":
    main()