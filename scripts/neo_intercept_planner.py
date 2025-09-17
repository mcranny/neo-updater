#!/usr/bin/env python3
from __future__ import annotations
"""
NEO Intercept Planner (SBDB + Lambert + Viewer-ready JSON)

- Reads hazardous-NEO JSON
- Fetches orbital elements from JPL SBDB
- Picks depart/arrive epochs, computes r1 (Earth) and r2 (NEO)
- Tries Lambert (via scripts.orbital if available), samples a 3D polyline
- Writes viewer-ready fields used by app/viewer3d_pg.py
"""

import argparse, json, math, sys, time, re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Ensure we can import scripts/orbital if present
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# --- Constants ---
AU_M   = 149_597_870_700.0  # meters per AU
MU_SUN = 1.32712440018e20   # m^3/s^2
DEG    = math.pi / 180.0

# --- Optional helpers from your repo (if present) ---
try:
    # lambert_universal(r1_m, r2_m, tof_s, mu=...) -> (v1, v2)
    from scripts.orbital import lambert_universal as _lambert   # type: ignore
except Exception:
    _lambert = None

try:
    # sample_transfer_polyline(r1_m, v1, mu, tof_s, n) -> (N,3) meters
    from scripts.orbital import sample_transfer_polyline as _sample_poly  # type: ignore
except Exception:
    _sample_poly = None


# ----------------- Time helpers -----------------
def jd_from_unix(t: float) -> float:
    return t / 86400.0 + 2440587.5

def unix_from_jd(jd: float) -> float:
    return (jd - 2440587.5) * 86400.0


# ----------------- SBDB fetch -----------------
def _sbdb_fetch(des_or_name: str, *, debug: bool=False) -> Optional[Dict[str, Any]]:
    """Ask SBDB for a single object. Returns dict with 'orbit' if ok.
       Adds retries and a 7-day on-disk cache to survive timeouts."""
    import requests, json, time, hashlib
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry
    cache_dir = Path("data/sbdb_cache"); cache_dir.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha1(des_or_name.encode("utf-8")).hexdigest()[:16]
    cpath = cache_dir / f"{key}.json"

    # cache hit (7 days)
    try:
        if cpath.exists() and (time.time() - cpath.stat().st_mtime) < 7*86400:
            if debug: print(f"[sbdb] cache hit for {des_or_name!r}")
            return json.loads(cpath.read_text(encoding="utf-8"))
    except Exception:
        pass

    url = "https://ssd-api.jpl.nasa.gov/sbdb.api"
    params = {"sstr": des_or_name, "full-prec": "true", "soln-epoch": "true"}

    sess = requests.Session()
    retry = Retry(
        total=4, backoff_factor=0.6,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
        raise_on_status=False,
    )
    sess.mount("https://", HTTPAdapter(max_retries=retry))
    try:
        r = sess.get(url, params=params, timeout=20)
    except Exception as e:
        if debug: print(f"[sbdb] request error for {des_or_name!r}: {e}")
        return None
    if r.status_code != 200:
        if debug: print(f"[sbdb] HTTP {r.status_code} for {des_or_name!r}")
        return None
    data = r.json()
    if not isinstance(data, dict) or "orbit" not in data:
        if debug: print(f"[sbdb] no 'orbit' in payload for {des_or_name!r}")
        return None

    try:
        cpath.write_text(json.dumps(data), encoding="utf-8")
    except Exception:
        pass
    return data

def _make_candidates(row: Dict[str, Any]) -> List[str]:
    cands: list[str] = []
    name = (row.get("name") or "").strip()
    if name:
        clean = name.replace("(", "").replace(")", "").strip()
        if clean: cands.append(clean)         # "2001 SY269"
        ns = clean.replace(" ", "")
        if ns and ns != clean: cands.append(ns)  # "2001SY269"
        cands.append(name)                        # "(2001 SY269)"
    des = (row.get("designation") or "").strip()
    if des:
        if des not in cands: cands.append(des)
        ns = des.replace(" ", "")
        if ns not in cands: cands.append(ns)
    for key in ("neo_reference_id", "id"):
        v = (row.get(key) or "").strip()
        if v and v not in cands:
            cands.append(v)
    # de-dup preserving order
    seen=set(); out=[]
    for s in cands:
        if s and s not in seen:
            seen.add(s); out.append(s)
    return out

def _sbdb_fetch_best(row: Dict[str, Any], *, debug: bool=False) -> Optional[Dict[str, Any]]:
    cands = _make_candidates(row)
    if debug: print(f"[sbdb] candidates: {cands}")
    for s in cands:
        sb = _sbdb_fetch(s, debug=debug)
        if sb is not None:
            if debug: print(f"[sbdb] success: {s!r}")
            return sb
    if debug: print("[sbdb] all candidates failed")
    return None


# ----------------- Parsing helpers -----------------
_num_re = re.compile(r"^[\s]*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")

def _num(x: Any) -> float:
    """Parse numbers from SBDB fields that may include units/strings."""
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x)
    m = _num_re.match(s)
    if m:
        return float(m.group(1))
    # space-delimited fallback
    return float(s.strip().split()[0])

def _as_elem_map(el: Any) -> Dict[str, Any]:
    """
    Normalize SBDB 'elements' into a flat mapping.
    - If dict: return as-is.
    - If list of {name,value,...}: build {name: value}.
    """
    if isinstance(el, dict):
        return dict(el)
    if isinstance(el, list):
        out: Dict[str, Any] = {}
        for item in el:
            if not isinstance(item, dict):
                continue
            key = (item.get("name") or item.get("label") or item.get("key") or "").strip()
            if not key:
                continue
            val = item.get("value")
            if val is None:
                val = item.get("val") or item.get("v")
            out[key] = val
        return out
    return {}

def _pick_map(m: Dict[str, Any], *keys: str) -> Optional[float]:
    """Try multiple key aliases (case-insensitive)."""
    for k in keys:
        if k in m:
            return _num(m[k])
        kl = k.lower()
        ku = k.upper()
        if kl in m: return _num(m[kl])
        if ku in m: return _num(m[ku])
    return None

def _parse_elements(sb: Dict[str, Any], *, debug: bool=False) -> Optional[Dict[str, float]]:
    """
    Parse SBDB 'orbit' -> 'elements' into:
      a_AU, e, i_deg, raan_deg, argp_deg, ma_deg, epoch_jd

    Accepts dict OR list-of-{name,value}. 
    Falls back to: 
      - a from q,e (a = q/(1-e))
      - ma from tp and derived mean motion n if 'ma' not given.
    Returns Dict[str, float] only (no None).
    """
    def _req(x: Optional[float], label: str) -> float:
        if x is None:
            raise ValueError(f"{label}=None")
        xf = float(x)
        if math.isnan(xf):
            raise ValueError(f"{label}=NaN")
        return xf

    try:
        orb = sb["orbit"]
        el_map = _as_elem_map(orb.get("elements"))

        # Primary fields (optionals at first)
        a_au  = _pick_map(el_map, "a")
        e     = _pick_map(el_map, "e")
        i_deg = _pick_map(el_map, "i", "incl", "inclination")
        raan  = _pick_map(el_map, "om", "Omega", "node", "RAAN", "LongNode")
        argp  = _pick_map(el_map, "w", "peri", "arg", "argPeri", "ArgP")
        ma    = _pick_map(el_map, "ma", "M", "mean_anom", "meanAnomaly")

        # Epoch (prefer orbit.epoch; else elements epoch/t0; else 'now')
        epoch_raw = orb.get("epoch", el_map.get("epoch", el_map.get("t0", el_map.get("Epoch"))))
        epoch_j   = _num(epoch_raw if epoch_raw is not None else jd_from_unix(time.time()))

        # If 'a' missing but 'q' and 'e' present, derive a = q/(1-e)
        if (a_au is None or (isinstance(a_au, float) and math.isnan(a_au))) and e is not None:
            q = _pick_map(el_map, "q", "q_au", "Perihelion", "q(AU)")
            if q is not None and e < 1.0:
                a_au = q / (1.0 - e)

        # If MA missing, try from mean motion and periapse passage
        if ma is None or (isinstance(ma, float) and math.isnan(ma)):
            tp = _pick_map(el_map, "tp", "T_p", "t_peri")
            n_deg_day: Optional[float] = _pick_map(el_map, "n", "meanMotion")
            if n_deg_day is None and a_au is not None:
                # derive n from a (Kepler 3rd), convert to deg/day
                a_m = a_au * AU_M
                n_rad_s = math.sqrt(MU_SUN / (a_m**3))
                n_deg_day = n_rad_s * 86400.0 * (180.0 / math.pi)
            if tp is not None and n_deg_day is not None:
                ma = (n_deg_day * (epoch_j - tp)) % 360.0
            else:
                ma = 0.0

        # Narrow everything to float (raises if None/NaN)
        a_f     = _req(a_au,  "a_AU")
        e_f     = _req(e,     "e")
        i_f     = _req(i_deg, "i_deg")
        raan_f  = _req(raan,  "raan_deg")
        argp_f  = _req(argp,  "argp_deg")
        ma_f    = _req(ma,    "ma_deg")
        epoch_f = _req(epoch_j, "epoch_jd")

        return {
            "a_AU": a_f, "e": e_f, "i_deg": i_f,
            "raan_deg": raan_f, "argp_deg": argp_f,
            "ma_deg": ma_f, "epoch_jd": epoch_f,
        }

    except Exception as exc:
        if debug:
            try:
                keys_seen = sorted(_as_elem_map(sb.get("orbit", {}).get("elements")).keys())
            except Exception:
                keys_seen = []
            print(f"[parse] exception: {exc}  keys_seen={keys_seen}")
        return None

# ----------------- Kepler / state conversion -----------------
def kepler_E_from_M(M: float, e: float, tol=1e-12, maxit=50) -> float:
    """Solve Kepler's equation M = E - e sin E (radians)"""
    M = (M + math.pi) % (2*math.pi) - math.pi
    E = M if e < 0.8 else math.pi
    for _ in range(maxit):
        f  = E - e*math.sin(E) - M
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
    Two-body propagation: elements at epoch0 -> r,v at epoch t_jd (meters, heliocentric ecliptic-ish).
    """
    a   = a_au * AU_M
    i   = i_deg * DEG
    Om  = raan_deg * DEG
    w   = argp_deg * DEG
    M0  = ma0_deg * DEG

    n   = math.sqrt(mu / (a**3))
    dt  = (t_jd - epoch0_jd) * 86400.0
    M   = M0 + n * dt
    E   = kepler_E_from_M(M, e)

    cosE, sinE = math.cos(E), math.sin(E)
    r_p  = a * (1 - e*cosE)
    x_p  = a * (cosE - e)
    y_p  = a * math.sqrt(1 - e*e) * sinE
    rdot = (math.sqrt(mu*a) / r_p) * (-sinE)
    rfdot= (math.sqrt(mu*a) / r_p) * (math.sqrt(1 - e*e) * cosE)

    cO,sO = math.cos(Om), math.sin(Om)
    ci,si = math.cos(i),  math.sin(i)
    cw,sw = math.cos(w),  math.sin(w)

    # PQW -> IJK rotation
    R11 = cO*cw - sO*sw*ci; R12 = -cO*sw - sO*cw*ci; R13 = sO*si
    R21 = sO*cw + cO*sw*ci; R22 = -sO*sw + cO*cw*ci; R23 = -cO*si
    R31 = sw*si;             R32 = cw*si;             R33 = ci
    R   = np.array([[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]], dtype=float)

    r_pf = np.array([x_p, y_p, 0.0])
    v_pf = np.array([rdot, rfdot, 0.0])
    r = R @ r_pf
    v = R @ v_pf
    return r, v


# ----------------- Earth elements (fallback) -----------------
EARTH_ELEMENTS_J2000 = {
    "a_AU": 1.00000011, "e": 0.01671022, "i_deg": 0.00005,
    "raan_deg": -11.26064, "argp_deg": 102.94719,
    "ma_deg": 100.46435, "epoch_jd": 2451545.0
}


# ----------------- Polyline helpers -----------------
def _ensure_nd(poly_like: Any) -> Optional[np.ndarray]:
    """Coerce to np.ndarray (N,3) float64 or return None."""
    if poly_like is None:
        return None
    arr = np.asarray(poly_like, dtype=float)
    if arr.ndim != 2 or arr.shape[1] != 3:
        return None
    return arr

def try_lambert_polyline(r1_m: np.ndarray, r2_m: np.ndarray, tof_s: float, n: int) -> Optional[np.ndarray]:
    """
    Attempt to compute a Lambert arc and sample it into an (N,3) ndarray in meters.
    Always returns np.ndarray or None (never raw lists).
    """
    if _lambert is None:
        return None
    try:
        v1, v2 = _lambert(r1_m, r2_m, tof_s, mu=MU_SUN)  # type: ignore
        if _sample_poly is not None:
            poly_m = _sample_poly(r1_m, v1, mu=MU_SUN, tof_s=tof_s, n=int(max(2, n)))  # type: ignore
            return _ensure_nd(poly_m)
        # Fallback sampler: straight line in R^3 (typed & visual, not physical)
        t = np.linspace(0.0, 1.0, int(max(2, n)))[:, None]
        poly_m = r1_m[None, :] * (1.0 - t) + r2_m[None, :] * t
        return poly_m
    except Exception:
        return None

def straight_polyline(r1_m: np.ndarray, r2_m: np.ndarray, n: int) -> np.ndarray:
    t = np.linspace(0.0, 1.0, int(max(2, n)))[:, None]
    return r1_m[None, :] * (1.0 - t) + r2_m[None, :] * t


# ----------------- Core planning -----------------
def build_plan_for_one(row: Dict[str, Any], depart_jd: float, tof_days: float, poly_n: int, *, debug: bool=False) -> Optional[Dict[str, Any]]:
    sb = _sbdb_fetch_best(row, debug=debug)
    if sb is None:
        if debug: print("[planner] SBDB fetch failed")
        return None
    el_tgt = _parse_elements(sb, debug=debug)
    if el_tgt is None:
        if debug: print("[planner] element parsing failed")
        return None

    # Earth at departure (propagate from J2000)
    el_earth = dict(EARTH_ELEMENTS_J2000)
    r1_m, _ = rv_from_elements_at_time(
        el_earth["a_AU"], el_earth["e"], el_earth["i_deg"],
        el_earth["raan_deg"], el_earth["argp_deg"], el_earth["ma_deg"],
        el_earth["epoch_jd"], depart_jd, mu=MU_SUN
    )

    # Target at arrival
    arrive_jd = depart_jd + tof_days
    r2_m, _ = rv_from_elements_at_time(
        el_tgt["a_AU"], el_tgt["e"], el_tgt["i_deg"],
        el_tgt["raan_deg"], el_tgt["argp_deg"], el_tgt["ma_deg"],
        el_tgt["epoch_jd"], arrive_jd, mu=MU_SUN
    )

    tof_s = float(tof_days) * 86400.0
    poly_m = try_lambert_polyline(r1_m, r2_m, tof_s, n=poly_n)
    if poly_m is None:
        poly_m = straight_polyline(r1_m, r2_m, n=poly_n)

    plan: Dict[str, Any] = {
        "elements_earth": {
            "a_AU": el_earth["a_AU"], "e": el_earth["e"],
            "i_deg": el_earth["i_deg"], "raan_deg": el_earth["raan_deg"],
            "argp_deg": el_earth["argp_deg"], "ma_deg": el_earth["ma_deg"],
            "epoch_jd": el_earth["epoch_jd"]
        },
        "elements_target": {
            "a_AU": el_tgt["a_AU"], "e": el_tgt["e"],
            "i_deg": el_tgt["i_deg"], "raan_deg": el_tgt["raan_deg"],
            "argp_deg": el_tgt["argp_deg"], "ma_deg": el_tgt["ma_deg"],
            "epoch_jd": el_tgt["epoch_jd"]
        },
        "r1_au": (r1_m / AU_M).tolist(),
        "r2_au": (r2_m / AU_M).tolist(),
        "lambert_poly_xyz_au": (poly_m / AU_M).tolist(),
        "departure_jd": depart_jd,
        "arrival_jd": arrive_jd,
        "tof_days": float(tof_days),
    }
    return plan


# ----------------- Main -----------------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp",  default="data/hazardous_neos/latest.json",
                    help="Input JSON with potentially_hazardous_neos or neos list.")
    ap.add_argument("--out",  default="data/hazardous_neos/latest_intercept.json",
                    help="Output path (can overwrite --inp).")
    ap.add_argument("--with-sbdb", action="store_true",
                    help="Fetch elements from JPL SBDB (recommended).")
    ap.add_argument("--polyline-n", type=int, default=600,
                    help="Samples along transfer polyline.")
    ap.add_argument("--depart-days", type=float, default=90.0,
                    help="Days from now for departure.")
    ap.add_argument("--tof-days", type=float, default=180.0,
                    help="Time of flight (days).")
    ap.add_argument("--debug", action="store_true",
                    help="Verbose SBDB identifier + parsing diagnostics.")
    args = ap.parse_args()

    inp = Path(args.inp)
    if not inp.exists():
        print(f"[planner] Input not found: {inp}")
        return 2

    try:
        j = json.loads(inp.read_text(encoding="utf-8"))
    except Exception as e:
        print("[planner] Failed to read JSON:", e)
        return 2

    rows: List[Dict[str, Any]] = (j.get("potentially_hazardous_neos")
                                  or j.get("neos") or [])
    if not rows:
        print("[planner] No NEOs found in input.")
        return 0

    now_jd    = jd_from_unix(time.time())
    depart_jd = now_jd + float(args.depart_days)
    updated   = 0

    for row in rows:
        plan: Dict[str, Any] = row.setdefault("intercept_plan", {})

        if not args.with_sbdb:
            # Minimal placeholder so the viewer can still render something
            plan.setdefault("elements_earth", dict(EARTH_ELEMENTS_J2000))
            plan.setdefault("elements_target", {
                "a_AU": 1.2, "e": 0.1, "i_deg": 5.0,
                "raan_deg": 50.0, "argp_deg": 10.0,
                "ma_deg": 0.0, "epoch_jd": 2451545.0
            })
            plan.setdefault("departure_jd", depart_jd)
            plan.setdefault("tof_days", float(args.tof_days))
            plan.setdefault("arrival_jd", depart_jd + float(args.tof_days))
            continue

        try:
            new_plan = build_plan_for_one(
                row, depart_jd, float(args.tof_days), int(args.polyline_n), debug=args.debug
            )
            if new_plan:
                plan.update(new_plan)
                updated += 1
            else:
                name = row.get("name") or row.get("designation") or row.get("id") or "UNKNOWN"
                print(f"[planner] SBDB/Lambert failed for {name}; skipping.")
        except Exception as e:
            name = row.get("name") or row.get("designation") or row.get("id") or "UNKNOWN"
            print(f"[planner] {name}: error {e}")

    j["count"] = len(rows)
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(j, ensure_ascii=False, indent=2))
    print(f"[planner] Updated {updated}/{len(rows)} NEO(s) -> {outp}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
