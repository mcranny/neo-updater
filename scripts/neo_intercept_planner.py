# scripts/neo_intercept_planner.py
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, cast
from scripts.test_lambert import lambert_universal
import numpy as np
from app.solar_system_viewer import (
    PLANETARY_DATA,
    planet_pos_au,
    datetime_to_jd,
    GAUSSIAN_GRAVITATIONAL_CONSTANT,
)

# --- Import the working Lambert from your test module.
# Make sure scripts/__init__.py exists so this is a proper package import.
try:
    from scripts.test_lambert import lambert_universal  # type: ignore[attr-defined]
except Exception as e:  # pragma: no cover
    print("[planner] FATAL: couldn't import lambert_universal from scripts.test_lambert")
    raise

# --- Pull math/utilities from the viewer without circular imports.
# Ensure app/__init__.py exists; this is what fixes Pylance "unresolved".
from app.solar_system_viewer import (  # type: ignore
    PLANETARY_DATA,
    planet_pos_au,             # heliocentric position from Kepler elements (AU)
    datetime_to_jd,            # UTC datetime -> JD
    GAUSSIAN_GRAVITATIONAL_CONSTANT,
)

MU_SUN = GAUSSIAN_GRAVITATIONAL_CONSTANT ** 2  # AU^3 / day^2

# =============================================================================
# Utilities (robust element normalization, 2-body propagation sampler)
# =============================================================================

def _to_f(x: Any) -> Optional[float]:
    if x is None:
        return None
    try:
        return float(x)
    except Exception:
        return None

def _require_float(x: Any, field: str) -> float:
    v = _to_f(x)
    if v is None:
        raise ValueError(f"Missing/non-numeric '{field}'")
    return v

def _walk(o: Any) -> Iterable[Dict[str, Any]]:
    if isinstance(o, dict):
        yield o
        for v in o.values():
            yield from _walk(v)
    elif isinstance(o, list):
        for it in o:
            yield from _walk(it)

def _normalize_elements(src: Mapping[str, Any]) -> Optional[Dict[str, float]]:
    """Map assorted SBDB-like keys into the baseline format the viewer expects."""
    def pick(*names: str) -> Optional[Any]:
        for n in names:
            v = src.get(n)
            if v is not None:
                return v
        return None

    a     = pick("a", "a_AU")
    e     = pick("e")
    i     = pick("i", "i_deg")
    raan  = pick("raan", "raan_deg", "Ω", "omega_node_deg")
    argp  = pick("argp", "argp_deg", "ω", "peri_arg_deg")

    if None in (a, e, i, raan, argp):
        return None

    out: Dict[str, float] = {
        "a":    _require_float(a,    "a"),
        "e":    _require_float(e,    "e"),
        "i":    _require_float(i,    "i"),
        "raan": _require_float(raan, "raan"),
        "argp": _require_float(argp, "argp"),
    }

    # Optional (for the moving NEO marker)
    ma    = pick("ma", "ma_deg", "M", "M_deg", "mean_anomaly", "mean_anomaly_deg")
    epoch = pick("epoch", "epoch_jd", "t_epoch", "epoch_tdb_jd", "epoch_tt_jd")
    if ma is not None:
        out["ma"] = _require_float(ma, "ma")
    if epoch is not None:
        out["epoch"] = _require_float(epoch, "epoch")

    # Perihelion time (lets us back out an epoch if only M is given)
    tp = pick("tp", "tp_jd", "t_peri", "t_peri_jd", "t_perihelion_jd")
    if tp is not None:
        out["_tp"] = _require_float(tp, "tp")

    if ("ma" in out) and ("epoch" not in out) and ("_tp" in out):
        a_au = out["a"]
        n = GAUSSIAN_GRAVITATIONAL_CONSTANT / (a_au ** 1.5)  # rad/day
        M = math.radians(out["ma"])
        out["epoch"] = out["_tp"] + (M / n)

    return out

# ---- Universal Kepler propagator (for sampling the transfer curve) ----

def _stumpff_C(z: float) -> float:
    if z > 0:
        s = math.sqrt(z)
        return (1.0 - math.cos(s)) / z
    if z < 0:
        s = math.sqrt(-z)
        return (1.0 - math.cosh(s)) / z
    return 0.5

def _stumpff_S(z: float) -> float:
    if z > 0:
        s = math.sqrt(z)
        return (s - math.sin(s)) / (s ** 3)
    if z < 0:
        s = math.sqrt(-z)
        return (math.sinh(s) - s) / (s ** 3)
    return 1.0 / 6.0

def propagate_universal(r0: np.ndarray, v0: np.ndarray, dt: float, mu: float = MU_SUN) -> Tuple[np.ndarray, np.ndarray]:
    """
    2-body propagation using universal variables.
    r0 [AU], v0 [AU/day], dt [day], mu [AU^3/day^2]
    """
    r0 = np.asarray(r0, float).reshape(3)
    v0 = np.asarray(v0, float).reshape(3)

    r0n = float(np.linalg.norm(r0))
    if r0n == 0.0:
        raise ValueError("propagate_universal: r0 is zero vector")
    v0n2 = float(np.dot(v0, v0))

    alpha = 2.0 / r0n - v0n2 / mu
    sqrt_mu = math.sqrt(mu)
    r0dotv0 = float(np.dot(r0, v0))

    # Initial guess for chi
    if abs(alpha) > 1e-12:
        chi = sqrt_mu * abs(alpha) * dt
    else:
        # near-parabolic
        h = np.cross(r0, v0)
        p = float(np.dot(h, h)) / mu
        s = 0.5 * (math.pi / 2.0 - math.atan(3.0 * math.sqrt(mu / (p ** 3)) * dt))
        w = math.atan(math.tan(s) ** (1.0 / 3.0))
        chi = math.sqrt(p) * (2.0 / math.tan(2.0 * w))

    # Newton iterations
    for _ in range(60):
        z = alpha * chi * chi
        C = _stumpff_C(z)
        S = _stumpff_S(z)
        r = chi * chi * C + r0dotv0 / sqrt_mu * chi * (1.0 - z * S) + r0n * (1.0 - z * C)
        F = r0n * r0dotv0 / sqrt_mu * chi * chi * C + (1.0 - alpha * r0n) * chi * chi * chi * S + r0n * chi - sqrt_mu * dt
        if abs(F) < 1e-12:
            break
        dF = r0n * r0dotv0 / sqrt_mu * chi * (1.0 - z * S) + (1.0 - alpha * r0n) * chi * chi * C + r0n
        chi -= F / dF

    # Lagrange f, g
    z = alpha * chi * chi
    C = _stumpff_C(z)
    S = _stumpff_S(z)
    f = 1.0 - (chi * chi / r0n) * C
    g = dt - (chi ** 3 / sqrt_mu) * S

    r_vec = f * r0 + g * v0
    rn = float(np.linalg.norm(r_vec))
    # fdot, gdot
    fdot = (sqrt_mu / (rn * r0n)) * (z * S - 1.0) * chi
    gdot = 1.0 - (chi * chi / rn) * C
    v_vec = fdot * r0 + gdot * v0
    return r_vec, v_vec

def sample_transfer_polyline(r0: np.ndarray, v0: np.ndarray, tof_days: float, n: int, mu: float = MU_SUN) -> np.ndarray:
    """Uniformly sample the conic from (r0,v0) over 'tof_days'. Returns (n,3) in AU."""
    n = max(2, int(n))
    ts = np.linspace(0.0, float(tof_days), n, dtype=float)
    pts = np.empty((n, 3), dtype=float)
    for i, t in enumerate(ts):
        r, _ = propagate_universal(r0, v0, t, mu)
        pts[i, :] = r
    return pts.astype(np.float32)

# =============================================================================
# Core planner
# =============================================================================

@dataclass
class Args:
    inp: str
    out: str
    depart_days: float
    tof_days: float
    polyline_n: int
    debug: bool

def _now_jd() -> float:
    import datetime as _dt
    return datetime_to_jd(_dt.datetime.utcnow())

def _load_latest(path: str) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    with p.open("r", encoding="utf-8") as f:
        return cast(Dict[str, Any], json.load(f))

def _save_out(path: str, payload: Dict[str, Any]) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(path + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)

def _get_targets(doc: Dict[str, Any]) -> List[Dict[str, Any]]:
    arr = doc.get("objects")
    if isinstance(arr, list) and arr:
        return [cast(Dict[str, Any], o) for o in arr if isinstance(o, dict)]
    cands: List[Dict[str, Any]] = []
    for d in _walk(doc):
        if "elements" in d and isinstance(d["elements"], dict):
            cands.append(d)
    return cands

def _earth_pos_au_at(jd: float) -> np.ndarray:
    return planet_pos_au(PLANETARY_DATA["Earth"]["elements"], jd)

def _target_pos_au_at(elems: Mapping[str, float], jd: float) -> np.ndarray:
    payload = {
        "a": elems["a"],
        "e": elems["e"],
        "i": elems["i"],
        "raan": elems["raan"],
        "argp": elems["argp"],
    }
    if "ma" in elems and "epoch" in elems:
        payload["ma"] = elems["ma"]
        payload["epoch"] = elems["epoch"]
    return planet_pos_au(payload, jd)

def plan_for_target(entry: Dict[str, Any], args: Args) -> Optional[Dict[str, Any]]:
    raw_elems = cast(Dict[str, Any], entry.get("elements", {}))
    elems = _normalize_elements(raw_elems)
    if not elems:
        if args.debug:
            print("[planner] skip: could not normalize target elements")
        return None

    depart_jd = _now_jd() + float(args.depart_days)
    tof_days = float(args.tof_days)
    if tof_days <= 0.0:
        if args.debug:
            print("[planner] Offline: TOF must be >0; coercing to 1 day")
        tof_days = 1.0
    arrive_jd = depart_jd + tof_days

    r1 = _earth_pos_au_at(depart_jd)
    r2 = _target_pos_au_at(elems, arrive_jd)

    if args.debug:
        print(f"[planner] r1 (Earth, AU) = {r1}")
        print(f"[planner] r2 (NEO,   AU) = {r2}")
        print(f"[planner] depart_jd={depart_jd:.9f}  arrive_jd={arrive_jd:.9f}  tof={tof_days:.3f} d")

    # Solve Lambert for AU/day velocities (short-way, prograde in test_lambert)
    try:
        v1_au_per_day, v2_au_per_day = lambert_universal(r1, r2, tof_days, MU_SUN)
    except Exception as e:
        if args.debug:
            print(f"[planner] Lambert failed: {e}")
        return {
            "target_elements": elems,
            "departure_jd": float(depart_jd),
            "arrival_jd": float(arrive_jd),
            "lambert_poly": None,
            "lambert_units": "au",
            "lambert_status": f"failed: {e}",
        }

    # Build a curved transfer polyline by propagating (r1, v1)
    n = max(50, int(args.polyline_n))
    poly_au = sample_transfer_polyline(r1, v1_au_per_day, tof_days, n, mu=MU_SUN)

    plan: Dict[str, Any] = {
        "target_elements": elems,               # normalized for the viewer
        "departure_jd": float(depart_jd),
        "arrival_jd": float(arrive_jd),
        "lambert_poly": poly_au.tolist(),      # AU
        "lambert_units": "au",
        "lambert_status": "ok",
        "v_depart_au_per_day": [float(x) for x in np.asarray(v1_au_per_day, float).reshape(3)],
        "v_arrive_au_per_day": [float(x) for x in np.asarray(v2_au_per_day, float).reshape(3)],
        "polyline_n": n,
    }
    return plan

# =============================================================================
# CLI
# =============================================================================

@dataclass
class ArgsNS:
    pass  # only for type clarity with argparse Namespace

def parse_cli(argv: List[str]) -> Args:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", required=True, help="Input JSON, e.g. data/hazardous_neos/latest.json")
    ap.add_argument("--out", required=True, help="Output JSON to write the planning result")
    ap.add_argument("--depart-days", type=float, default=90.0, dest="depart_days")
    ap.add_argument("--tof-days", type=float, default=180.0, dest="tof_days")
    ap.add_argument("--polyline-n", type=int, default=600, dest="polyline_n")
    ap.add_argument("--debug", action="store_true")
    ns = ap.parse_args(argv)
    return Args(
        inp=ns.inp,
        out=ns.out,
        depart_days=ns.depart_days,
        tof_days=ns.tof_days,
        polyline_n=ns.polyline_n,
        debug=bool(ns.debug),
    )

def main(argv: List[str]) -> int:
    args = parse_cli(argv)

    doc = _load_latest(args.inp)
    targets = _get_targets(doc)
    if not targets:
        print("[planner] No targets found in input")
        return 2

    out_payload: Dict[str, Any] = {"objects": []}
    updated = 0
    for t in targets:
        plan = plan_for_target(t, args)
        if plan is not None:
            out_payload["objects"].append(plan)
            updated += 1

    _save_out(args.out, out_payload)
    print(f"[planner] Updated {updated}/{len(targets)} NEO(s) -> {args.out}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
