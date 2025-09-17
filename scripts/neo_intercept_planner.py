#!/usr/bin/env python3
from __future__ import annotations
"""
NEO Intercept Planner (uses local scripts/orbital.py)

- Robust SBDB fetch (retries + 7-day cache) with numeric-ID first
- ALWAYS outputs a curved Lambert arc:
    * try universal-variable Lambert (prograde, retrograde)
    * if it fails, fall back to a shooting solver with your universal propagator
- Exports 3-D and 2-D polylines (AU), and departure/arrival markers
"""
import argparse, json, math, sys, time, re, hashlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# Use your local orbital helpers
ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from scripts.orbital import (  # type: ignore
    lambert_universal,
    kepler_universal_propagate,
    MU_SUN,
)

# ---------- constants ----------
AU_M = 149_597_870_700.0
DEG  = math.pi / 180.0

# ---------- time helpers ----------
def jd_from_unix(t: float) -> float:
    return t / 86400.0 + 2440587.5

def unix_from_jd(jd: float) -> float:
    return (jd - 2440587.5) * 86400.0

def utc_from_jd(jd: float) -> str:
    import datetime as _dt
    return _dt.datetime.fromtimestamp(unix_from_jd(jd), tz=_dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

# ---------- SBDB fetch (retries + cache) ----------
def _sbdb_fetch(des_or_name: str, *, debug: bool=False) -> Optional[Dict[str, Any]]:
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    cache_dir = Path("data/sbdb_cache"); cache_dir.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha1(des_or_name.encode("utf-8")).hexdigest()[:16]
    cpath = cache_dir / f"{key}.json"

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

# ---------- parsing helpers (bullet-proof) ----------
_num_re = re.compile(r"^[\s]*([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")

def _num(x: Any) -> float:
    """Robust float parser: None/'None'/'' → nan; parse leading numeric token."""
    if x is None:
        return float('nan')
    if isinstance(x, (int, float)):
        return float(x)
    s = str(x).strip()
    if s == "" or s.lower() in {"none", "nan", "null", "-"}:
        return float('nan')
    m = _num_re.match(s)
    try:
        return float(m.group(1)) if m else float(s)
    except Exception:
        return float('nan')

def _as_elem_map(el: Any) -> Dict[str, Any]:
    if isinstance(el, dict): return el
    if isinstance(el, list):
        out = {}
        for d in el:
            try: out[str(d.get("name"))] = d.get("value")
            except Exception: pass
        return out
    return {}

def _parse_elements(sbdb: Dict[str, Any], *, debug: bool=False) -> Optional[Dict[str, float]]:
    try:
        orbit = sbdb.get("orbit") or {}
        el = _as_elem_map(orbit.get("elements") or orbit)

        def g(*keys: str) -> float:
            for k in keys:
                if k in el:
                    return _num(el.get(k))
            return float('nan')

        # Epoch
        epoch_jd = _num(orbit.get("epoch_jd"))
        if math.isnan(epoch_jd):
            epoch_jd = g("epoch", "epoch_jd")
        if math.isnan(epoch_jd):
            epoch_jd = 2451545.0

        # Elements (AU / deg)
        a_AU   = g("a")
        e      = g("e")
        i_deg  = g("i")
        raan   = g("om", "Omega", "raan")
        argp   = g("w", "argp", "peri")
        ma_deg = g("ma", "M", "mean_anom")

        # Derive if missing
        q = g("q")
        Q = g("ad", "Q")
        if math.isnan(a_AU):
            if not math.isnan(q) and not math.isnan(e):
                a_AU = q / max(1e-12, (1.0 - e))
            elif not math.isnan(q) and not math.isnan(Q):
                a_AU = 0.5 * (q + Q)
        if math.isnan(e) and not math.isnan(q) and not math.isnan(Q) and (q + Q) > 0:
            e = (Q - q) / (Q + q)
        if math.isnan(ma_deg):
            n_deg_day = g("n")
            tp_jd = g("tp", "tp_jd", "t_p_jd", "tperi", "T", "T_jd")
            if not math.isnan(n_deg_day) and not math.isnan(tp_jd):
                ma_deg = (n_deg_day * (epoch_jd - tp_jd)) % 360.0

        # Fill angles if still NaN
        if math.isnan(i_deg):  i_deg  = 0.0
        if math.isnan(raan):   raan   = 0.0
        if math.isnan(argp):   argp   = 0.0
        if math.isnan(ma_deg): ma_deg = 0.0

        def wrap360(x: float) -> float:
            x %= 360.0
            return x + 360.0 if x < 0 else x

        return {
            "a_AU":     float(a_AU),
            "e":        float(e),
            "i_deg":    wrap360(float(i_deg)),
            "raan_deg": wrap360(float(raan)),
            "argp_deg": wrap360(float(argp)),
            "ma_deg":   wrap360(float(ma_deg)),
            "epoch_jd": float(epoch_jd),
        }
    except Exception as e:
        if debug: print("[parse] failed:", e)
        return None

# ---------- elements -> state ----------
def kepler_E_from_M(M: float, e: float, tol=1e-12, maxit=50) -> float:
    M = (M + math.pi) % (2*math.pi) - math.pi
    E = M if e < 0.8 else math.pi
    for _ in range(maxit):
        f  = E - e*math.sin(E) - M
        fp = 1 - e*math.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            break
    return E

def rv_from_elements_at_time(a_AU, e, i_deg, raan_deg, argp_deg, ma0_deg, epoch0_jd, t_jd, mu=MU_SUN) -> Tuple[np.ndarray, np.ndarray]:
    a = a_AU * AU_M
    i = i_deg * DEG; Om = raan_deg * DEG; w = argp_deg * DEG; M0 = ma0_deg * DEG
    n = math.sqrt(mu / (a**3))
    dt = (t_jd - epoch0_jd) * 86400.0
    M = M0 + n * dt
    E = kepler_E_from_M(M, e)

    cosE = math.cos(E); sinE = math.sin(E)
    fac  = math.sqrt(1 - e*e)
    r_p  = a * (1 - e*cosE)
    x_p  = a * (cosE - e)
    y_p  = a * fac * sinE
    rdot   = (math.sqrt(mu * a) / r_p) * e * sinE
    rfdot  = (math.sqrt(mu * a) / r_p) * fac * cosE

    cw, sw = math.cos(w), math.sin(w)
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(i),  math.sin(i)

    R11 = cO*cw - sO*sw*ci; R12 = -cO*sw - sO*cw*ci; R13 = sO*si
    R21 = sO*cw + cO*sw*ci; R22 = -sO*sw + cO*cw*ci; R23 = -cO*si
    R31 = sw*si;             R32 = cw*si;             R33 = ci
    R   = np.array([[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]], dtype=float)

    r_pf = np.array([x_p, y_p, 0.0])
    v_pf = np.array([rdot, rfdot, 0.0])
    r = R @ r_pf
    v = R @ v_pf
    return r, v

EARTH_ELEMENTS_J2000 = {
    "a_AU": 1.00000011, "e": 0.01671022, "i_deg": 0.00005,
    "raan_deg": -11.26064, "argp_deg": 102.94719,
    "ma_deg": 100.46435, "epoch_jd": 2451545.0
}

# ---------- SBDB candidates (numeric-first) ----------
def _make_candidates(row: Dict[str, Any]) -> List[str]:
    cands: List[str] = []
    for key in ("neo_reference_id", "id", "neo_ref_id"):
        v = (row.get(key) or "").strip()
        if v and v.isdigit():
            cands.append(v)
    for key in ("name", "designation"):
        val = (row.get(key) or "").strip()
        if not val: continue
        m = re.search(r"(\d{4,})", val)
        if m:
            num = m.group(1)
            if num not in cands:
                cands.insert(0, num)
        clean = val.replace("(", "").replace(")", "").strip()
        for s in (clean, clean.replace(" ", ""), val):
            if s and s not in cands:
                cands.append(s)
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
        else:
            if debug: print(f"[sbdb] failed: {s!r}")
    if debug: print("[sbdb] all candidates failed")
    return None

# ---------- Lambert polyline (UV → fallback shooting) ----------
def _force_endpoints(poly_m: np.ndarray, r1_m: np.ndarray, r2_m: np.ndarray) -> np.ndarray:
    poly_m = np.asarray(poly_m, dtype=float)
    poly_m[0, :]  = np.asarray(r1_m, dtype=float).reshape(3)
    poly_m[-1, :] = np.asarray(r2_m, dtype=float).reshape(3)
    return poly_m

def lambert_polyline_always(r1_m: np.ndarray, r2_m: np.ndarray, tof_s: float, n: int) -> np.ndarray:
    """
    1) try your lambert_universal (prograde, retrograde)
    2) on failure, use a shooting solver to find v1 s.t. r(tof; v1) = r2
    Returns (N,3) meters, endpoints forced to r1/r2.
    """
    r1_m = np.asarray(r1_m, dtype=float).reshape(3)
    r2_m = np.asarray(r2_m, dtype=float).reshape(3)
    tof_s = float(tof_s)
    n = int(max(2, n))

    # --- try universal-variable lambert (prograde / retrograde)
    last_err: Optional[Exception] = None
    for pro in (True, False):
        try:
            v1, _ = lambert_universal(r1_m, r2_m, tof_s, mu=MU_SUN, prograde=pro)
            ts = np.linspace(0.0, tof_s, n)
            pts = np.empty((n, 3), dtype=float)
            for i, t in enumerate(ts):
                r, _ = kepler_universal_propagate(r1_m, v1, float(t), MU_SUN)
                pts[i, :] = r
            return _force_endpoints(pts, r1_m, r2_m)
        except Exception as e:
            last_err = e

    # --- fallback: shooting solve in v1 space (Broyden-ish with finite differences)
    # Initial guess: chord velocity + a tangential component
    chord = (r2_m - r1_m) / max(1e-9, tof_s)
    # tangential unit at r1 (in-plane guess)
    zhat = np.cross(r1_m, r2_m)
    t_hat = np.cross(zhat, r1_m); t_n = np.linalg.norm(t_hat)
    if t_n > 0: t_hat /= t_n
    v_circ = math.sqrt(MU_SUN / max(1e-3, np.linalg.norm(r1_m))) * t_hat
    v = 0.6 * chord + 0.8 * v_circ

    def residual(v1: np.ndarray) -> np.ndarray:
        r, _ = kepler_universal_propagate(r1_m, v1, tof_s, MU_SUN)
        return (r - r2_m).astype(float)

    eps = 1e-3
    for _ in range(30):
        F = residual(v)
        if np.linalg.norm(F) < 1.0:  # 1 meter target
            break
        # finite-diff Jacobian
        J = np.zeros((3,3), dtype=float)
        for k in range(3):
            dv = np.zeros(3); dv[k] = eps
            J[:, k] = (residual(v + dv) - F) / eps
        # solve J * dv = -F
        try:
            dv, *_ = np.linalg.lstsq(J, -F, rcond=None)
        except Exception:
            dv = -F * 1e-6
        # damping
        v += 0.8 * dv

    # sample
    ts = np.linspace(0.0, tof_s, n)
    pts = np.empty((n, 3), dtype=float)
    for i, t in enumerate(ts):
        r, _ = kepler_universal_propagate(r1_m, v, float(t), MU_SUN)
        pts[i, :] = r
    return _force_endpoints(pts, r1_m, r2_m)

# ---------- plan build ----------
def build_plan_for_one(row: Dict[str, Any], depart_jd: float, tof_days: float, poly_n: int, *, debug: bool=False) -> Optional[Dict[str, Any]]:
    sb = _sbdb_fetch_best(row, debug=debug)
    if sb is None:
        if debug: print("[planner] SBDB fetch failed")
        return None
    el_tgt = _parse_elements(sb, debug=debug)
    if el_tgt is None:
        if debug: print("[planner] element parsing failed")
        return None

    # Earth at departure
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

    # Curved Lambert (with fallback shooting)
    tof_s = float(tof_days) * 86400.0
    poly_m = lambert_polyline_always(r1_m, r2_m, tof_s, n=int(poly_n))

    poly_au  = (poly_m / AU_M).astype(float)
    poly_xy  = poly_au[:, :2]

    plan: Dict[str, Any] = {
        "elements_earth": dict(el_earth),
        "elements_target": dict(el_tgt),
        "r1_au": (r1_m / AU_M).tolist(),
        "r2_au": (r2_m / AU_M).tolist(),
        "lambert_poly_xyz_au": poly_au.tolist(),
        "lambert_polyline_xyz_au": poly_au.tolist(),
        "lambert_polyline_xy_au": poly_xy.tolist(),
        "departure_jd": depart_jd,
        "arrival_jd": arrive_jd,
        "tof_days": float(tof_days),
        "departure_utc": utc_from_jd(depart_jd),
        "arrival_utc":   utc_from_jd(arrive_jd),
    }
    return plan

# ---------- main ----------
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp",  default="data/hazardous_neos/latest.json")
    ap.add_argument("--out",  default="data/hazardous_neos/latest_intercept.json")
    ap.add_argument("--with-sbdb", action="store_true", help="Fetch target elements from SBDB")
    ap.add_argument("--polyline-n", type=int, default=800)
    ap.add_argument("--depart-days", type=float, default=90.0)
    ap.add_argument("--tof-days", type=float, default=180.0)
    ap.add_argument("--debug", action="store_true")
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

    rows: List[Dict[str, Any]] = (j.get("potentially_hazardous_neos") or j.get("neos") or [])
    if not rows:
        print("[planner] No NEOs found in input.")
        return 0

    now_jd    = jd_from_unix(time.time())
    depart_jd = now_jd + float(args.depart_days)

    updated = 0
    for row in rows:
        plan: Dict[str, Any] = row.setdefault("intercept_plan", {})

        if not args.with_sbdb:
            # OFFLINE: synthesize target elements but still build a curved arc
            el_earth = dict(EARTH_ELEMENTS_J2000)
            el_tgt = plan.get("elements_target") or {
                "a_AU": 1.2, "e": 0.10, "i_deg": 5.0,
                "raan_deg": 50.0, "argp_deg": 10.0,
                "ma_deg": 0.0, "epoch_jd": 2451545.0
            }
            arrive_jd = depart_jd + float(args.tof_days)

            r1_m, _ = rv_from_elements_at_time(
                el_earth["a_AU"], el_earth["e"], el_earth["i_deg"],
                el_earth["raan_deg"], el_earth["argp_deg"], el_earth["ma_deg"],
                el_earth["epoch_jd"], depart_jd, mu=MU_SUN
            )
            r2_m, _ = rv_from_elements_at_time(
                el_tgt["a_AU"], el_tgt["e"], el_tgt["i_deg"],
                el_tgt["raan_deg"], el_tgt["argp_deg"], el_tgt["ma_deg"],
                el_tgt["epoch_jd"], arrive_jd, mu=MU_SUN
            )

            poly_m = lambert_polyline_always(r1_m, r2_m, float(args.tof_days)*86400.0, n=int(args.polyline_n))
            poly_au = (poly_m / AU_M).astype(float)

            plan.update({
                "elements_earth": dict(el_earth),
                "elements_target": dict(el_tgt),
                "r1_au": (r1_m / AU_M).tolist(),
                "r2_au": (r2_m / AU_M).tolist(),
                "lambert_poly_xyz_au": poly_au.tolist(),
                "lambert_polyline_xyz_au": poly_au.tolist(),
                "lambert_polyline_xy_au": poly_au[:, :2].tolist(),
                "departure_jd": depart_jd,
                "arrival_jd": arrive_jd,
                "tof_days": float(args.tof_days),
                "departure_utc": utc_from_jd(depart_jd),
                "arrival_utc":   utc_from_jd(arrive_jd),
            })
            updated += 1
            continue

        # SBDB path
        p = build_plan_for_one(row, depart_jd, float(args.tof_days), int(args.polyline_n), debug=args.debug)
        if p is not None:
            row["intercept_plan"] = p
            updated += 1

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(j, indent=2), encoding="utf-8")
    print(f"[planner] Updated {updated}/{len(rows)} NEO(s) -> {args.out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
