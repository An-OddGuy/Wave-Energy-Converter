#!/usr/bin/env python3
"""check_mesh.py - sanity-check the exported STL meshes (Stage 0).

A bad mesh gives wrong hydrodynamic coefficients without an obvious error, so run
this after every make_cad.py run.

Usage:
    python check_mesh.py                 # checks cad_out/
    python check_mesh.py --dir cad_out --tol 0.01

Checks per mesh: watertight, consistent winding, positive volume (outward normals),
and volume agreement with the CAD volume in derived.yaml (tessellation makes the mesh
volume slightly smaller than the exact one, so a small tolerance is normal).
Exit code is 0 if everything passes, 1 otherwise.
"""
import argparse
import sys
from pathlib import Path

import trimesh
import yaml

HERE = Path(__file__).resolve().parent   # so the script works from any working directory


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=str(HERE / "cad_out"))
    ap.add_argument("--tol", type=float, default=0.02,
                    help="allowed relative volume error vs CAD (default 2%%)")
    args = ap.parse_args()

    d = Path(args.dir)
    derived_path = d / "derived.yaml"
    cad_vol = {}
    if derived_path.exists():
        with open(derived_path, "r", encoding="utf-8") as f:
            cad_vol = (yaml.safe_load(f) or {}).get("cad_volumes_m3", {})
    else:
        print(f"WARNING: {derived_path} not found; skipping volume comparison")

    stls = sorted(d.glob("*.stl"))
    if not stls:
        print(f"No STL files found in {d}. Run make_cad.py first.")
        sys.exit(1)

    ok_all = True
    print(f"{'file':<18}{'watertight':<12}{'winding':<9}{'volume m^3':<14}"
          f"{'vs CAD':<10}{'faces':<8}status")
    for path in stls:
        m = trimesh.load_mesh(str(path), force="mesh")
        wt = bool(m.is_watertight)
        wind = bool(m.is_winding_consistent)
        vol = float(m.volume) if wt else float("nan")
        status = []
        if not wt:
            status.append("NOT watertight")
        if not wind:
            status.append("inconsistent winding")
        if wt and vol <= 0:
            status.append("negative volume (normals point inward)")

        rel = ""
        base = path.stem.replace("_local", "")
        if wt and base in cad_vol and cad_vol[base] > 0:
            err = abs(vol - cad_vol[base]) / cad_vol[base]
            rel = f"{err * 100:.2f}%"
            if err > args.tol:
                status.append(f"volume off by {err * 100:.1f}%")

        ok = not status
        ok_all &= ok
        print(f"{path.name:<18}{str(wt):<12}{str(wind):<9}{vol:<14.6f}"
              f"{rel:<10}{len(m.faces):<8}{'OK' if ok else '; '.join(status)}")

    print()
    if ok_all:
        print("All meshes passed.")
        sys.exit(0)
    print("Some meshes failed. Fix the CAD parameters or tessellation tolerance "
          "(export.linear_tolerance in params.yaml) before running Capytaine.")
    sys.exit(1)


if __name__ == "__main__":
    main()
