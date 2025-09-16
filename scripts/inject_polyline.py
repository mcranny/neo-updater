#!/usr/bin/env python3
from __future__ import annotations
import argparse, json
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Tuple
import numpy as np

AU = 149_597_870_700.0

# ---------- unit/shape helpers ----------
def _as_au(vec_or_pts: Any) -> np.ndarray:
    a = np.asarray(vec_or_pts, float)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    if a.size == 0:
        return np.zeros((1,3))
    if a.shape[1] == 2:  # pad Z
        a = np.c_[a, np.zeros(len(a))]
    # heuristics: meters -> AU
    if float(np.max(np.abs(a))) > 10.0:
        a = a / AU
    return a

def _is_vec3(x: Any) -> bool:
    try:
        a = np.asarray(x, float).reshape(-1)
        return a.size in (2,3)
    except Exception:
        return False

# ---------- rotations for element-based synthetic fallback ----------
def _rot3(t):
    c,s=np.cos(t),np.sin(t)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]], float)
def _rot1(t):
    c,s=np.cos(t),np.sin(t)
    return np.array([[1,0,0],[0,c,-s],[0,s,c]], float)

# Try lots of possible key names
_R1_KEYS = [
    "r1_au","r1_m","r1","r_depart","r_depart_au","r_depart_m",
    "r_start","r_start_au","r_earth","rE","r_dep","r_departure"
]
_R2_KEYS = [
    "r2_au","r2_m","r2","r_arrive","r_arrive_au","r_arrive_m",
    "r_end","r_end_au","r_target","rA","r_arrival"
]

def _find_vec(plan: Dict[str, Any], keys: Iterable[str]) -> Optional[np.ndarray]:
    # 1) exact keys
    for k in keys:
        if k in plan and _is_vec3(plan[k]):
            return _as_au(plan[k])
    # 2) any vec-looking values
    for k,v in plan.items():
        if _is_vec3(v) and any(tag in k.lower() for tag in ("r1","r2","start","end","depart","arriv","earth","target")):
            return _as_au(v)
    return None

def _synthesize_from_elements(plan: Dict[str, Any]) -> Optional[Tuple[np.ndarray,np.ndarray]]:
    elems = plan.get("elements") or {}
    if not elems:
        return None
    a_au = float(elems.get("a_AU") or elems.get("a") or elems.get("a_au") or 1.2)
    e    = float(elems.get("e", 0.0))
    i    = float(elems.get("i_deg") or elems.get("i") or 0.0) * np.pi/180.0
    raan = float(elems.get("raan_deg") or elems.get("RAAN_deg") or elems.get("raan") or 0.0) * np.pi/180.0
    argp = float(elems.get("argp_deg") or elems.get("omega_deg") or elems.get("argp") or 0.0) * np.pi/180.0
    # Use simple radii along x in each frame (not physically correct, just to visualize)
    rE = np.array([[1.0, 0.0, 0.0]])            # Earth ~1 AU at +X
    rp = np.array([[a_au, 0.0, 0.0]])           # target radius at periapsis direction
    R  = _rot3(raan) @ _rot1(i) @ _rot3(argp)
    rT = (R @ rp.T).T
    return _as_au(rE), _as_au(rT)

def inject(inp: Path, outp: Path, n: int) -> Tuple[int,int,int]:
    j = json.loads(inp.read_text())
    rows = j.get("potentially_hazardous_neos", []) or j.get("neos", [])
    added_r12, added_from_elems, skipped = 0, 0, 0
    for row in rows:
        plan = row.setdefault("intercept_plan", {})

        poly = (plan.get("lambert_polyline_xyz_au")
                or plan.get("lambert_polyline")
                or plan.get("lambert_polyline_m"))
        if isinstance(poly, list) and len(poly) >= 2:
            continue  # already has a path

        r1 = _find_vec(plan, _R1_KEYS)
        r2 = _find_vec(plan, _R2_KEYS)

        if r1 is None or r2 is None:
            synth = _synthesize_from_elements(plan)
            if synth is None:
                skipped += 1
                continue
            r1, r2 = synth
            added_from_elems += 1
        else:
            added_r12 += 1

        r1 = r1.reshape(3); r2 = r2.reshape(3)
        t = np.linspace(0.0, 1.0, max(2, n))
        pts = (r1[None,:]*(1-t[:,None]) + r2[None,:]*t[:,None]).tolist()
        plan["lambert_polyline_xyz_au"] = pts

    outp.write_text(json.dumps(j, indent=2))
    return added_r12, added_from_elems, skipped

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", default="data/hazardous_neos/latest_intercept.json")
    ap.add_argument("--out", default="data/hazardous_neos/latest_intercept_poly.json")
    ap.add_argument("--n", type=int, default=400)
    args = ap.parse_args()
    a, b, s = inject(Path(args.inp), Path(args.out), args.n)
    print(f"Injected from r1/r2: {a}  |  from elements: {b}  |  skipped: {s}  -> {args.out}")
