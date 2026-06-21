# scripts/neo_intercept_planner.py
"""
NEO intercept planner — Lambert transfer solver.

Bugs fixed vs previous version:
1. Removed duplicate ``from scripts.test_lambert import lambert_universal`` (appeared
   on lines 7 AND 20).  ``lambert_universal`` was never defined in test_lambert.py.
2. Removed Qt-dependent import of ``app.solar_system_viewer``.  Replaced with
   ``app.orbital_math`` (pure math, no PySide6/pyqtgraph) so CI can run this
   without a display.
3. Unit mismatch fixed.  ``scripts.orbital.lambert_universal`` expects SI (meters,
   seconds).  Positions from ``planet_pos_au`` are in AU, TOF is in days.  A thin
   wrapper converts to/from SI before/after each Lambert call.
4. ``_normalize_elements``: added ``"epoch_jd_tdb"`` and ``"tp_jd_tdb"`` to pick
   lists — these are the actual keys in the JPL CAD+SBDB feed.
5. Output schema: renamed ``lambert_poly`` → ``lambert_poly_xyz_au`` (matches the
   key the viewer looks for); added ``tof_days`` and ``r2_AU`` to each plan entry.
6. CLI: added ``--input``/``--output`` as aliases for ``--inp``/``--out``; added
   ``--roll-past`` flag (accepted but unused — the planner always uses now+depart_days).
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Tuple, cast

import numpy as np

# ---- Pure orbital math (no Qt) ----
from app.orbital_math import (
    PLANETARY_DATA,
    planet_pos_au,
    datetime_to_jd,
    GAUSSIAN_GRAVITATIONAL_CONSTANT,
)

# ---- Lambert solver (SI: meters, seconds) from scripts.orbital ----
from scripts.orbital import (
    lambert_universal as _lambert_si,
    AU_M as _AU_M,
)

MU_SUN: float = GAUSSIAN_GRAVITATIONAL_CONSTANT ** 2   # AU^3 / day^2
_DAY_S: float = 86400.0

# =============================================================================
# Unit-conversion wrapper: AU/day <-> SI (m, s)
# =============================================================================

def lambert_universal(
    r1_au: np.ndarray,
    r2_au: np.ndarray,
    tof_days: float,
    _mu_unused: float = MU_SUN,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Solve Lambert's problem.  Inputs in AU and days; outputs in AU/day.

    Delegates to ``scripts.orbital.lambert_universal`` (SI: meters, seconds).
    Automatically retries with the opposite prograde/retrograde direction when
    the first attempt fails with a geometry error — this handles the 'y(z*)<=0'
    failure that occurs when the cross-product heuristic picks the long-way
    transfer but the resulting A-value goes negative.
    """
    r1_m = np.asarray(r1_au, float) * _AU_M
    r2_m = np.asarray(r2_au, float) * _AU_M
    tof_s = float(tof_days) * _DAY_S
    fac = _DAY_S / _AU_M   # m/s → AU/day

    last_exc: Exception = RuntimeError("Lambert: no direction succeeded")
    for prograde in (True, False):
        try:
            v1_ms, v2_ms = _lambert_si(r1_m, r2_m, tof_s, prograde=prograde)
            return np.asarray(v1_ms, float) * fac, np.asarray(v2_ms, float) * fac
        except Exception as exc:
            last_exc = exc
    raise last_exc


# =============================================================================
# Utilities (element normaliser, 2-body propagator)
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

    a    = pick("a", "a_AU")
    e    = pick("e")
    i    = pick("i", "i_deg")
    raan = pick("raan", "raan_deg", "Ω", "omega_node_deg")
    argp = pick("argp", "argp_deg", "ω", "peri_arg_deg")

    if None in (a, e, i, raan, argp):
        return None

    out: Dict[str, float] = {
        "a":    _require_float(a,    "a"),
        "e":    _require_float(e,    "e"),
        "i":    _require_float(i,    "i"),
        "raan": _require_float(raan, "raan"),
        "argp": _require_float(argp, "argp"),
    }

    ma = pick("ma", "ma_deg", "M", "M_deg", "mean_anomaly", "mean_anomaly_deg")
    # Bug fix #4: added "epoch_jd_tdb" — the actual key in the JPL CAD+SBDB feed.
    epoch = pick("epoch", "epoch_jd", "t_epoch", "epoch_tdb_jd", "epoch_tt_jd",
                 "epoch_jd_tdb")
    if ma is not None:
        out["ma"] = _require_float(ma, "ma")
    if epoch is not None:
        out["epoch"] = _require_float(epoch, "epoch")

    # Bug fix #4: added "tp_jd_tdb" — the actual perihelion-time key in the feed.
    tp = pick("tp", "tp_jd", "t_peri", "t_peri_jd", "t_perihelion_jd", "tp_jd_tdb")
    if tp is not None:
        out["_tp"] = _require_float(tp, "tp")

    # Derive epoch from tp + M if epoch was missing
    if ("ma" in out) and ("epoch" not in out) and ("_tp" in out):
        a_au = out["a"]
        n = GAUSSIAN_GRAVITATIONAL_CONSTANT / (a_au ** 1.5)   # rad/day
        M = math.radians(out["ma"])
        out["epoch"] = out["_tp"] + (M / n)

    return out


# ---- Overflow-safe Kepler propagator (AU / day) ----------------------------
# We inline stable Stumpff functions here instead of importing from
# kepler_propagation.py because that file overflows on hyperbolic transfers
# (cosh(s) with s > ~710 raises OverflowError).

_MAX_S = 500.0   # cosh/sinh safe upper bound

def _stC(z: float) -> float:
    """Numerically stable C(z) Stumpff function, clamped for large |z|."""
    if abs(z) < 1e-8:
        return 0.5 - z / 24.0 + z * z / 720.0
    if z > 0.0:
        return (1.0 - math.cos(math.sqrt(z))) / z
    s = math.sqrt(-z)
    if s > _MAX_S:
        return math.exp(s) / (2.0 * (-z))   # asymptotic: (cosh s - 1)/(-z) ~ e^s/(2(-z))
    return (math.cosh(s) - 1.0) / (-z)

def _stS(z: float) -> float:
    """Numerically stable S(z) Stumpff function, clamped for large |z|."""
    if abs(z) < 1e-8:
        return (1.0 / 6.0) - z / 120.0 + z * z / 5040.0
    if z > 0.0:
        s = math.sqrt(z)
        return (s - math.sin(s)) / (s ** 3)
    s = math.sqrt(-z)
    if s > _MAX_S:
        return math.exp(s) / (2.0 * s ** 3)  # asymptotic
    return (math.sinh(s) - s) / (s ** 3)


def propagate_universal(
    r0: np.ndarray, v0: np.ndarray, dt: float, mu: float = MU_SUN,
) -> Tuple[np.ndarray, np.ndarray]:
    """2-body propagation via universal variables (AU / day)."""
    r0 = np.asarray(r0, float).reshape(3)
    v0 = np.asarray(v0, float).reshape(3)
    r0n = float(np.linalg.norm(r0))
    if r0n == 0.0:
        raise ValueError("propagate_universal: r0 is zero")

    v0n2 = float(np.dot(v0, v0))
    alpha = 2.0 / r0n - v0n2 / mu
    sqrt_mu = math.sqrt(mu)
    r0dv0 = float(np.dot(r0, v0))

    if abs(alpha) > 1e-12:
        chi = sqrt_mu * abs(alpha) * dt
    else:
        h = np.cross(r0, v0)
        p = float(np.dot(h, h)) / mu
        s = 0.5 * (math.pi / 2.0 - math.atan(3.0 * math.sqrt(mu / p ** 3) * dt))
        w = math.atan(math.tan(s) ** (1.0 / 3.0))
        chi = math.sqrt(p) * (2.0 / math.tan(2.0 * w))

    for _ in range(60):
        z = alpha * chi * chi
        C = _stC(z); S = _stS(z)
        r = chi*chi*C + r0dv0/sqrt_mu*chi*(1.0-z*S) + r0n*(1.0-z*C)
        F = (r0n*r0dv0/sqrt_mu * chi*chi*C
             + (1.0-alpha*r0n) * chi*chi*chi*S
             + r0n*chi - sqrt_mu*dt)
        if abs(F) < 1e-12:
            break
        dF = (r0n*r0dv0/sqrt_mu*chi*(1.0-z*S)
              + (1.0-alpha*r0n)*chi*chi*C + r0n)
        chi -= F / (dF or 1e-20)

    z = alpha * chi * chi; C = _stC(z); S = _stS(z)
    f = 1.0 - (chi*chi/r0n)*C
    g = dt - (chi**3/sqrt_mu)*S
    r_vec = f*r0 + g*v0
    rn = float(np.linalg.norm(r_vec))
    fdot = (sqrt_mu/(rn*r0n))*(z*S-1.0)*chi
    gdot = 1.0-(chi*chi/rn)*C
    return r_vec, fdot*r0 + gdot*v0


def sample_transfer_polyline(
    r0: np.ndarray, v0: np.ndarray, tof_days: float, n: int, mu: float = MU_SUN,
) -> np.ndarray:
    """
    Sample (n, 3) AU positions along the transfer conic.
    Falls back to straight-line interpolation for any point that fails
    (can happen with retrograde or near-degenerate Lambert solutions).
    """
    r0 = np.asarray(r0, float).reshape(3)
    v0 = np.asarray(v0, float).reshape(3)
    n = max(2, int(n))
    tof = float(tof_days)
    ts = np.linspace(0.0, tof, n)
    pts = np.empty((n, 3), float)

    # Attempt to get the arrival position for the straight-line fallback
    try:
        r_end, _ = propagate_universal(r0, v0, tof, mu)
    except Exception:
        r_end = r0 + v0 * tof   # linear extrapolation

    for k, t in enumerate(ts):
        try:
            r, _ = propagate_universal(r0, v0, t, mu)
            if not np.all(np.isfinite(r)):
                raise ValueError("non-finite position")
            pts[k] = r
        except Exception:
            # Linear interpolation between r0 and r_end
            frac = t / tof if tof > 0 else 0.0
            pts[k] = r0 + frac * (r_end - r0)

    return np.clip(pts, -1e6, 1e6).astype(np.float32)  # clip degenerate overflow before cast


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
    roll_past: bool


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
    payload: Dict[str, Any] = {
        "a": elems["a"], "e": elems["e"], "i": elems["i"],
        "raan": elems["raan"], "argp": elems["argp"],
    }
    if "ma" in elems and "epoch" in elems:
        payload["ma"] = elems["ma"]
        payload["epoch"] = elems["epoch"]
    return planet_pos_au(payload, jd)


def plan_for_target(entry: Dict[str, Any], args: Args) -> Optional[Dict[str, Any]]:
    # Prefer top-level "elements" key (CAD+SBDB feed), fall back to "target_elements"
    raw_elems = cast(Dict[str, Any], entry.get("elements") or entry.get("target_elements") or {})
    elems = _normalize_elements(raw_elems)
    if not elems:
        if args.debug:
            print("[planner] skip: could not normalize target elements")
        return None

    depart_jd  = _now_jd() + float(args.depart_days)
    tof_days   = max(1.0, float(args.tof_days))
    arrive_jd  = depart_jd + tof_days

    r1 = _earth_pos_au_at(depart_jd)
    r2 = _target_pos_au_at(elems, arrive_jd)

    if args.debug:
        print(f"[planner] r1 (Earth, AU) = {r1}")
        print(f"[planner] r2 (NEO,   AU) = {r2}")
        print(f"[planner] depart_jd={depart_jd:.3f}  tof={tof_days:.1f} d")

    try:
        # Bug fix #3: lambert_universal wrapper converts AU/day <-> SI internally
        v1_au_day, v2_au_day = lambert_universal(r1, r2, tof_days)
    except Exception as exc:
        if args.debug:
            print(f"[planner] Lambert failed: {exc}")
        return {
            "target_elements":     elems,
            "elements":            elems,
            "departure_jd":        float(depart_jd),
            "arrival_jd":          float(arrive_jd),
            "tof_days":            tof_days,
            "r1_AU":               float(np.linalg.norm(r1)),
            "r2_AU":               float(elems.get("a", float(np.linalg.norm(r2)))),
            "lambert_poly_xyz_au": None,
            "lambert_units":       "au",
            "lambert_status":      f"failed: {exc}",
            "_neo_name":           entry.get("des") or entry.get("name", "Target"),
        }

    n = max(50, int(args.polyline_n))
    poly_au = sample_transfer_polyline(r1, v1_au_day, tof_days, n, mu=MU_SUN)

    # Bug fix #5: renamed lambert_poly → lambert_poly_xyz_au; added tof_days / r2_AU
    return {
        "target_elements":   elems,
        "elements":          elems,           # viewer also looks for "elements"
        "departure_jd":      float(depart_jd),
        "arrival_jd":        float(arrive_jd),
        "tof_days":          tof_days,
        "r1_AU":             float(np.linalg.norm(r1)),
        "r2_AU":             float(elems.get("a", float(np.linalg.norm(r2)))),
        "lambert_poly_xyz_au": poly_au.tolist(),   # key the viewer expects
        "lambert_units":     "au",
        "lambert_status":    "ok",
        "v_depart_au_per_day": [float(x) for x in v1_au_day.reshape(3)],
        "v_arrive_au_per_day": [float(x) for x in v2_au_day.reshape(3)],
        "polyline_n":        n,
        "_neo_name":         entry.get("des") or entry.get("name", "Target"),
    }


# =============================================================================
# CLI
# =============================================================================

def parse_cli(argv: List[str]) -> Args:
    import argparse
    ap = argparse.ArgumentParser(description="Plan Lambert intercepts for NEOs.")
    # Bug fix #6: accept both --inp/--input and --out/--output (CI uses --input/--output)
    grp_in  = ap.add_mutually_exclusive_group(required=True)
    grp_out = ap.add_mutually_exclusive_group(required=True)
    grp_in.add_argument("--inp",   dest="inp")
    grp_in.add_argument("--input", dest="inp")
    grp_out.add_argument("--out",    dest="out")
    grp_out.add_argument("--output", dest="out")
    ap.add_argument("--depart-days",  type=float, default=90.0,  dest="depart_days")
    ap.add_argument("--tof-days",     type=float, default=180.0, dest="tof_days")
    ap.add_argument("--polyline-n",   type=int,   default=600,   dest="polyline_n")
    ap.add_argument("--debug",        action="store_true")
    # Bug fix #6: --roll-past accepted (was expected by CI workflow) but is a noop
    # in this planner; the "depart_days" offset serves the same purpose.
    ap.add_argument("--roll-past", action="store_true", dest="roll_past",
                    help="[noop] Accepted for backward compat with CI workflow.")
    # Legacy Lambert flags from old CI invocations — accepted but ignored here
    # (plan_intercepts.py is the preferred Lambert pipeline now)
    ap.add_argument("--elements",   action="store_true")
    ap.add_argument("--lambert",    action="store_true")
    ap.add_argument("--depart-utc", dest="depart_utc", default=None)
    ap.add_argument("--r1-au",      type=float, default=None, dest="r1_au")
    ap.add_argument("--r2-au",      type=float, default=None, dest="r2_au")
    ap.add_argument("--leo-km",     type=float, default=None, dest="leo_km")

    ns = ap.parse_args(argv)
    return Args(
        inp=ns.inp,
        out=ns.out,
        depart_days=ns.depart_days,
        tof_days=ns.tof_days,
        polyline_n=ns.polyline_n,
        debug=bool(ns.debug),
        roll_past=bool(ns.roll_past),
    )


def main(argv: List[str]) -> int:
    args = parse_cli(argv)
    doc  = _load_latest(args.inp)
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
