#!/usr/bin/env python3
"""make_cad.py - generate the WEC digital-twin CAD from params.yaml (Stage 0).

Usage:
    python make_cad.py                       # uses params.yaml
    python make_cad.py --mat-depth 0.25      # override mat depth (m), e.g. for sweeps
    python make_cad.py --params other.yaml --out my_out

Frame: z up, z = 0 at still water level (SWL), units metres.

Outputs (in export.out_dir):
    <part>.step            neutral CAD (open in FreeCAD / SolidWorks)
    <part>.stl             world-frame mesh (use for Capytaine / hydrodynamics)
    <part>_local.stl       mesh centred on its bounding-box centre (use for CoppeliaSim)
    assembly.step          all parts in one file
    derived.yaml           mass, inertia, areas, nominal poses, sanity numbers
"""
import argparse
import math
import sys
from pathlib import Path

import yaml
import cadquery as cq


# ----------------------------------------------------------------- primitives
def cylinder(radius, height, z0):
    return cq.Workplane("XY").workplane(offset=z0).circle(radius).extrude(height)


def ring(r_in, r_out, thickness, z0):
    return (cq.Workplane("XY").workplane(offset=z0)
            .circle(r_out).circle(r_in).extrude(thickness))


def plate(length, width, thickness, x_off, z0):
    return (cq.Workplane("XY").workplane(offset=z0)
            .center(x_off, 0).rect(length, width).extrude(thickness))


# -------------------------------------------------------------------- helpers
def fail(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def bbox_center(wp):
    bb = wp.val().BoundingBox()
    return ((bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2, (bb.zmin + bb.zmax) / 2)


def validate(p):
    env, b, sh, ba, m = p["environment"], p["buoy"], p["shaft"], p["base"], p["mat"]
    if b["draft"] <= 0 or b["draft"] >= b["height"]:
        fail("buoy.draft must be > 0 and < buoy.height")
    if b["radius"] <= sh["radius"]:
        fail("buoy.radius must be larger than shaft.radius")
    if ba["radius"] <= sh["radius"]:
        fail("base.radius must be larger than shaft.radius")
    if b["draft"] + 0.05 > env["water_depth"]:
        fail("water_depth is too shallow for the buoy draft")
    if m["enabled"]:
        seabed_top = -env["water_depth"] + ba["thickness"]
        mat_bottom = -m["depth"] - m["thickness"]
        if m["depth"] < 0:
            fail("mat.depth must be >= 0 (top of mat below SWL)")
        if mat_bottom < seabed_top:
            fail("mat would hit the seabed base; reduce mat.depth or thickness")
        if not (m["depth_min"] <= m["depth"] <= m["depth_max"]):
            print("WARNING: mat.depth is outside [depth_min, depth_max]")
        if m["shape"] == "ring":
            if m["ring_inner_radius"] <= b["radius"]:
                fail("ring_inner_radius must exceed buoy.radius so the buoy can heave freely")
            if m["ring_outer_radius"] <= m["ring_inner_radius"]:
                fail("ring_outer_radius must exceed ring_inner_radius")
        elif m["shape"] == "plate":
            near_edge = m["plate_x_offset"] + m["plate_length"] / 2
            if near_edge > -b["radius"]:
                print("WARNING: plate overlaps the buoy footprint in x; "
                      "check clearance as the buoy heaves")
        else:
            fail("mat.shape must be 'ring' or 'plate'")


# ----------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--params", default="params.yaml")
    ap.add_argument("--mat-depth", type=float, default=None,
                    help="override mat.depth (m)")
    ap.add_argument("--out", default=None, help="override export.out_dir")
    args = ap.parse_args()

    with open(args.params, "r", encoding="utf-8") as f:
        p = yaml.safe_load(f)
    if args.mat_depth is not None:
        p["mat"]["depth"] = args.mat_depth
    validate(p)

    env, b, sh, ba, m, ex = (p["environment"], p["buoy"], p["shaft"],
                             p["base"], p["mat"], p["export"])
    out = Path(args.out or ex["out_dir"])
    out.mkdir(parents=True, exist_ok=True)

    # ---- parts (world frame, z = 0 at SWL)
    z_buoy_bottom = -b["draft"]
    z_buoy_top = b["height"] - b["draft"]
    z_base_bottom = -env["water_depth"]
    z_base_top = z_base_bottom + ba["thickness"]
    z_shaft_top = z_buoy_top + sh["top_above_buoy"]

    parts = {
        "buoy": cylinder(b["radius"], b["height"], z_buoy_bottom),
        "shaft": cylinder(sh["radius"], z_shaft_top - z_base_top, z_base_top),
        "base": cylinder(ba["radius"], ba["thickness"], z_base_bottom),
    }
    if m["enabled"]:
        z_mat_bottom = -m["depth"] - m["thickness"]
        if m["shape"] == "ring":
            parts["mat"] = ring(m["ring_inner_radius"], m["ring_outer_radius"],
                                m["thickness"], z_mat_bottom)
        else:
            parts["mat"] = plate(m["plate_length"], m["plate_width"],
                                 m["thickness"], m["plate_x_offset"], z_mat_bottom)

    # ---- export
    tol, atol = ex["linear_tolerance"], ex["angular_tolerance"]
    nominal_pose = {}
    for name, wp in parts.items():
        cq.exporters.export(wp, str(out / f"{name}.step"))
        cq.exporters.export(wp, str(out / f"{name}.stl"),
                            tolerance=tol, angularTolerance=atol)
        cx, cy, cz = bbox_center(wp)
        nominal_pose[name] = {"x": round(cx, 6), "y": round(cy, 6), "z": round(cz, 6)}
        local = wp.translate((-cx, -cy, -cz))
        cq.exporters.export(local, str(out / f"{name}_local.stl"),
                            tolerance=tol, angularTolerance=atol)

    compound = cq.Compound.makeCompound([wp.val() for wp in parts.values()])
    cq.exporters.export(compound, str(out / "assembly.step"))

    # ---- derived quantities (analytic, uniform-density buoy)
    rho, g = env["rho_water"], env["gravity"]
    r, h, d = b["radius"], b["height"], b["draft"]
    area_wp = math.pi * r * r
    mass = b["mass"] if b["mass"] is not None else rho * area_wp * d
    displaced_at_draft = area_wp * d
    k_heave = rho * g * area_wp                      # hydrostatic stiffness, N/m
    t_heave_dry = 2 * math.pi * math.sqrt(mass / k_heave)   # no added mass
    ixx = mass * (3 * r * r + h * h) / 12.0
    izz = mass * r * r / 2.0
    cog_z = z_buoy_bottom + h / 2.0                  # uniform-density assumption

    derived = {
        "note": "Analytic values assume a uniform-density buoy. Real ballast moves the CoG.",
        "buoy": {
            "mass_kg": round(mass, 6),
            "cog_z_m": round(cog_z, 6),
            "inertia_Ixx_Iyy_kgm2": round(ixx, 8),
            "inertia_Izz_kgm2": round(izz, 8),
            "waterplane_area_m2": round(area_wp, 8),
            "displaced_volume_at_design_draft_m3": round(displaced_at_draft, 8),
            "cad_volume_m3": round(parts["buoy"].val().Volume(), 8),
            "analytic_volume_m3": round(area_wp * h, 8),
            "heave_stiffness_N_per_m": round(k_heave, 4),
            "heave_natural_period_no_added_mass_s": round(t_heave_dry, 4),
            "equilibrium_check": ("neutrally buoyant at design draft"
                                  if abs(mass - rho * displaced_at_draft) < 1e-9
                                  else "mass differs from displaced mass; draft will shift"),
        },
        "mat": ({
            "enabled": True, "shape": m["shape"], "depth_m": m["depth"],
            "depth_min_m": m["depth_min"], "depth_max_m": m["depth_max"],
            "cad_volume_m3": round(parts["mat"].val().Volume(), 8),
        } if m["enabled"] else {"enabled": False}),
        "cad_volumes_m3": {k: round(v.val().Volume(), 8) for k, v in parts.items()},
        "nominal_pose_bbox_centre_m": nominal_pose,
        "coppelia_note": ("Import <part>_local.stl (metres) and place each object at its "
                          "nominal_pose; then drive: buoy z += heave, mat z = nominal z - "
                          "(depth - depth_nominal)."),
    }
    with open(out / "derived.yaml", "w", encoding="utf-8") as f:
        yaml.safe_dump(derived, f, sort_keys=False)

    # ---- summary
    print(f"Wrote {len(parts)} parts + assembly to {out.resolve()}")
    print(f"  buoy mass            : {mass:.3f} kg")
    print(f"  heave stiffness      : {k_heave:.1f} N/m")
    print(f"  heave period (no A33): {t_heave_dry:.3f} s  (added mass will lengthen it)")
    if m["enabled"]:
        print(f"  mat ({m['shape']}) top depth : {m['depth']:.3f} m below SWL")
    print("Next: python check_mesh.py")


if __name__ == "__main__":
    main()
