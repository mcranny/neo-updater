#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, sys, math
from pathlib import Path
from typing import Any, Dict, List

REQUIRED_KEYS = ["schema","source","objects"]

def au_km(au: float) -> float: return au * 149_597_870.7

def main():
    ap = argparse.ArgumentParser(description="Sanity-check a CAD+SBDB feed.")
    ap.add_argument("path", help="Path to latest_intercepts.json")
    args = ap.parse_args()

    p = Path(args.path)
    if not p.exists():
        print(f"ERROR: file not found: {p}", file=sys.stderr)
        return 2

    j = json.loads(p.read_text())
    for k in REQUIRED_KEYS:
        if k not in j:
            print(f"ERROR: missing top-level key: {k}", file=sys.stderr)
            return 2

    objs: List[Dict[str, Any]] = j.get("objects", [])
    if not objs:
        print("ERROR: no objects in feed", file=sys.stderr)
        return 2

    n = len(objs)
    n_hud = 0
    n_intercepts = 0
    n_warnings = 0

    for o in objs:
        des = o.get("des","<unknown>")
        ca  = o.get("next_ca") or {}
        el  = o.get("elements") or {}
        hud = o.get("hud")
        itc = o.get("intercept")

        # Basic fields present
        missing = [k for k in ("des","next_ca","elements") if o.get(k) is None]
        if missing:
            print(f"WARNING [{des}]: missing keys {missing}")
            n_warnings += 1

        # HUD ≈ CAD distance check (tolerant; Earth circular model)
        if hud and ca and ca.get("dist_au") is not None:
            n_hud += 1
            d_hud = float(hud.get("earth_distance_au") or 0.0)
            d_cad = float(ca.get("dist_au") or 0.0)
            if d_hud > 0 and d_cad > 0:
                abs_err = abs(d_hud - d_cad)
                tol = max(0.01, 0.05 * d_cad)   # 0.01 AU absolute or 5% relative
                if abs_err > tol:
                    print(f"WARNING [{des}]: HUD vs CAD Earth distance mismatch "
                          f"(hud={d_hud:.6f} au, cad={d_cad:.6f} au, err={abs_err:.6f} au > tol={tol:.6f} au)")
                    n_warnings += 1

        # Intercept quick plausibility
        if itc:
            n_intercepts += 1
            dv_dep = float(itc.get("dv_depart_kms") or 0.0)
            dv_arr = float(itc.get("dv_arrive_kms") or 0.0)
            if dv_dep <= 0 or dv_arr < 0 or dv_dep > 50 or dv_arr > 50:
                print(f"WARNING [{des}]: suspect Δv (dep={dv_dep:.2f} km/s, arr={dv_arr:.2f} km/s)")
                n_warnings += 1
            poly = itc.get("lambert_polyline_xy_au") or []
            if len(poly) < 10:
                print(f"WARNING [{des}]: short transfer polyline (len={len(poly)})")
                n_warnings += 1

    print(f"OK: {n} objects | HUD on {n_hud} | intercepts on {n_intercepts} | warnings: {n_warnings}")
    # return nonzero only on zero objects or fatal structure issues (handled earlier)
    return 0

if __name__ == "__main__":
    sys.exit(main())
