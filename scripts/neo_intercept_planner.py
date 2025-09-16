#!/usr/bin/env python3
from __future__ import annotations
"""
NEO Intercept Planner (SBDB + Lambert + Viewer-ready JSON)
- Reads your existing hazardous-NEO JSON
- Fetches orbital elements from JPL SBDB
- Picks depart/arrive epochs, computes r1 (Earth) and r2 (NEO)
- Tries Lambert (via scripts.orbital if available), samples a polyline
- Writes intercept_plan fields used by app/viewer3d.py

Usage (from repo root, venv active):
  python -m scripts.neo_intercept_planner \
    --inp data/hazardous_neos/latest_intercept.json \
    --out data/hazardous_neos/latest_intercept.json \
    --with-sbdb --polyline-n 600 --depart-days 90 --tof-days 180
"""

import argparse, json, math, os, sys, time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Ensure we can import scripts/orbital if present
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# --- Try optional helpers from your repo ---
AU_M = 149_597_870_700.0  # meters per AU
MU_SUN = 1.32712440018e20  # m^3/s^2

try:
    # If you have these in scripts/orbital, we’ll use them
    from scripts.orbital import sample_transfer_polyline as _sample_poly
except Exception:
    _sample_poly = None

try:
    # If you ship a Lambert solver, we’ll use it
    from scripts.orbital import lambert_universal as _lambert
except Exception:
    _lambert = None


# ----------------- Time / epochs -----------------
def jd_from_unix(t: float) -> float:
    """Unix seconds -> Julian Day (UTC-ish)."""
    return t / 86400.0 + 2440587.5

def unix_from_jd(jd: float) -> float:
    return (jd - 2440587.5) * 86400.0


# ----------------- SBDB fetch -----------------
# We avoid astroquery to keep deps light; plain requests works.
def _sbdb_fetch(des_or_name: str) -> Optional[Dict[str, Any]]:
    """
    Ask SBDB for a single object. Returns dict with 'orbit'->'elements' if ok.
    """
    import requests
    url = "https://ssd-api.jpl.nasa.gov/sbdb.api"
    params = {
        "sstr": des_or_name,     # name/designation/number
        "full-prec": "true",     # more digits
        "soln-epoch": "true"     # return elements at solution epoch
    }
    r = requests.get(url, params=params, timeout=20)
    if r.status_code != 200:
        return None
    data = r.json()
    # Expected keys: object, orbit{ epoch, elements{ a,e,i,om,w,ma,n,tp,per... } }
    if not isinstance(data, dict) or "orbit" not in data:
        return None
    return data


def _parse_elements(sb: Dict[str, Any]) -> Optional[Dict[str, float]]:
    """
    Parse SBDB response into a simple dict:
      a_AU, e, i_deg, raan_deg, argp_deg, ma_deg, epoch_jd
    """
    try:
        orb = sb["orbit"]
        el = orb["elements"]
        def _num(x):
            # SBDB usually returns numeric strings; strip units if any
            if isinstance(x, (int, float)): return float(x)
            s = str(x)
            # keep leading number
            for cut in (" ", "AU", "deg", "d", "rad", "/ d"):
                s = s.replace(cut, " ")
            return float(s.strip().split()[0])

        a_au = _num(el["a"])
        e = _num(el["e"])
        i_deg = _num(el["i"])
        raan_deg = _num(el["om"])
        argp_deg = _num(el["w"])
        ma_deg = _num(el.get("ma", 0.0))
        # Epoch may live in orbit['epoch'] (Julian days)
        epoch_jd = _num(orb.get("epoch", jd_from_unix(time.time())))
        return {
            "a_AU": a_au, "e": e, "i_deg": i_deg,
            "raan_deg": raan_deg, "argp_deg": argp_deg,
            "ma_deg": ma_deg, "epoch_jd": epoch_jd,
        }
    except Exception:
        return None


# ----------------- Kepler / state conversion (self-contained) -----------------
DEG = math.pi / 180.0

def kepler_E_from_M(M: float, e: float, tol=1e-12, maxit=50) -> float:
    """Solve Kepler's equation M = E - e sin E (radians)"""
    M = (M + math.pi) % (2*math.pi) - math.pi
    E = M if e < 0.8 else math.pi
    for _ in range(maxit):
        f = E - e*math.sin(E) - M
        fp = 1 - e*math.cos(E)
        dE = -f/fp
        E += dE
        if abs(dE) < tol:
            break
    return E

def rv_from_elements_at_time(a_au: float, e: float, i_deg: float, raan_deg: float, argp_deg: float,
                             ma0_deg: float, epoch0_jd: float, t_jd: float,
                             mu=MU_SUN) -> Tuple[np.ndarray, np.ndarray]:
    """
    Two-body propagation: elements at epoch0 -> r,v at epoch t_jd (heliocentric ecliptic J2000-ish).
    Returns meters.
    """
    a = a_au * AU_M
    i = i_deg * DEG
    Om = raan_deg * DEG
    w = argp_deg * DEG
    M0 = ma0_deg * DEG

    # mean motion n [rad/s]; time since epoch
    n = math.sqrt(mu / (a**3))
    dt = (t_jd - epoch0_jd) * 86400.0
    M = M0 + n * dt

    E = kepler_E_from_M(M, e)
    cosE, sinE = math.cos(E), math.sin(E)
    # perifocal coords
    r_p = a * (1 - e*cosE)
    x_p = a * (cosE - e)
    y_p = a * math.sqrt(1 - e*e) * sinE
    # velocities in perifocal
    rdot = (math.sqrt(mu*a) / r_p) * (-sinE)
    rfdot = (math.sqrt(mu*a) / r_p) * (math.sqrt(1 - e*e) * cosE)

    # rotation matrices
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(i), math.sin(i)
    cw, sw = math.cos(w), math.sin(w)

    # PQW -> IJK
    R11 = cO*cw - sO*sw*ci
    R12 = -cO*sw - sO*cw*ci
    R13 = sO*si
    R21 = sO*cw + cO*sw*ci
    R22 = -sO*sw + cO*cw*ci
    R23 = -cO*si
    R31 = sw*si
    R32 = cw*si
    R33 = ci

    r_pf = np.array([x_p, y_p, 0.0])
    v_pf = np.array([rdot, rfdot, 0.0])

    R = np.array([[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]], dtype=float)
    r = R @ r_pf
    v = R @ v_pf
    return r, v


# ----------------- Earth ephemeris (simple) -----------------
def earth_r_at_jd(t_jd: float) -> np.ndarray:
    """
    Simple heliocentric Earth position (meters). Good enough for transfer viz.
    Circular-ish with small eccentricity ignored.
    """
    # Mean motion ~ 2pi / sidereal year
    n = 2.0*math.pi / (365.256363004 * 86400.0)
    # Choose L0 s.t. J2000 aligns near +X; exact phase not critical for viz.
    L0 = 0.0
    dt = (t_jd - 2451545.0) * 86400.0  # seconds since J2000
    th = L0 + n * dt
    x = AU_M * math.cos(th)
    y = AU_M * math.sin(th)
    z = 0.0
    return np.array([x,y,z], dtype=float)


# ----------------- Transfer construction -----------------
def try_lambert(r1_m: np.ndarray, r2_m: np.ndarray, tof_s: float) -> Optional[np.ndarray]:
    """
    If scripts.orbital.lambert_universal + sample_transfer_polyline exist,
    use them to generate a polyline in meters. Otherwise return None.
    """
    if _lambert is None:
        return None
    try:
        v1, v2 = _lambert(r1_m, r2_m, tof_s, mu=MU_SUN)
        if _sample_poly is not None:
            poly_m = _sample_poly(r1_m, v1, mu=MU_SUN, tof_s=tof_s, n=600)
        else:
            # fallback sampler: straight line in state space (good enough)
            t = np.linspace(0, tof_s, 600)[:,None]
            poly_m = r1_m[None,:] + (r2_m - r1_m)[None,:] * (t / tof_s)
        return poly_m
    except Exception:
        return None


def straight_polyline(r1_m: np.ndarray, r2_m: np.ndarray, n: int=400) -> np.ndarray:
    t = np.linspace(0.0, 1.0, max(2, n))[:,None]
    return r1_m[None,:]*(1-t) + r2_m[None,:]*t


# ----------------- Main planner -----------------
def build_plan_for_one(name_or_des: str, depart_jd: float, tof_days: float, poly_n: int) -> Optional[Dict[str, Any]]:
    sb = _sbdb_fetch(name_or_des)
    if sb is None:
        return None
    el = _parse_elements(sb)
    if el is None:
        return None

    # r1: Earth at depart; r2: target at arrival
    r1_m = earth_r_at_jd(depart_jd)
    arrive_jd = depart_jd + tof_days
    r2_m, _ = rv_from_elements_at_time(
        el["a_AU"], el["e"], el["i_deg"], el["raan_deg"], el["argp_deg"],
        el.get("ma_deg", 0.0), el.get("epoch_jd", depart_jd), arrive_jd, mu=MU_SUN
    )

    tof_s = tof_days * 86400.0
    poly_m = try_lambert(r1_m, r2_m, tof_s)
    if poly_m is None:
        poly_m = straight_polyline(r1_m, r2_m, n=poly_n)

    plan = {
        "elements": {
            "a_AU": el["a_AU"], "e": el["e"],
            "i_deg": el["i_deg"], "raan_deg": el["raan_deg"], "argp_deg": el["argp_deg"]
        },
        "r1_au": (r1_m / AU_M).tolist(),
        "r2_au": (r2_m / AU_M).tolist(),
        "lambert_polyline_xyz_au": (poly_m / AU_M).tolist(),
        "depart_epoch_jd": depart_jd,
        "arrive_epoch_jd": arrive_jd,
        "tof_days": tof_days
    }
    return plan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", default="data/hazardous_neos/latest_intercept.json",
                    help="Input JSON with potentially_hazardous_neos list.")
    ap.add_argument("--out", default="data/hazardous_neos/latest_intercept.json",
                    help="Output path (can overwrite --inp).")
    ap.add_argument("--with-sbdb", action="store_true", help="Fetch elements from JPL SBDB.")
    ap.add_argument("--polyline-n", type=int, default=600, help="Samples along transfer polyline.")
    ap.add_argument("--depart-days", type=float, default=90.0, help="Days from now for departure.")
    ap.add_argument("--tof-days", type=float, default=180.0, help="Time of flight (days).")
    args = ap.parse_args()

    inp = Path(args.inp)
    if not inp.exists():
        print(f"[planner] Input not found: {inp}")
        return 2

    try:
        j = json.loads(inp.read_text())
    except Exception as e:
        print("[planner] Failed to read JSON:", e)
        return 2

    rows: List[Dict[str, Any]] = (j.get("potentially_hazardous_neos")
                                  or j.get("neos") or [])
    if not rows:
        print("[planner] No NEOs found in input.")
        return 0

    now_jd = jd_from_unix(time.time())
    depart_jd = now_jd + float(args.depart_days)
    updated = 0

    for row in rows:
        ident = (row.get("name")
                 or row.get("designation")
                 or row.get("neo_reference_id")
                 or row.get("id")
                 or "").strip()
        if not ident:
            continue

        plan: Dict[str, Any] = row.setdefault("intercept_plan", {})
        if not args.with_sbdb:
            # Minimal: at least write Earth ring radius so viewer shows something
            plan.setdefault("elements", {"a_AU": 1.2, "e": 0.1,
                                         "i_deg": 5.0, "raan_deg": 50.0, "argp_deg": 10.0})
            continue

        try:
            new_plan = build_plan_for_one(ident, depart_jd, float(args.tof_days), int(args.polyline_n))
            if new_plan:
                plan.update(new_plan)
                updated += 1
            else:
                print(f"[planner] SBDB/Lambert failed for {ident}; skipping.")
        except Exception as e:
            print(f"[planner] {ident}: error {e}")

    j["count"] = len(rows)
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(j, indent=2))
    print(f"[planner] Updated {updated}/{len(rows)} NEO(s) -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
