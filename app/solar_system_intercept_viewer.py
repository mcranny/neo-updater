# app/solar_system_intercept_viewer.py
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, Iterable, cast, Literal

import numpy as np
from PySide6 import QtWidgets  # type: ignore
import pyqtgraph.opengl as gl  # type: ignore

# Baseline viewer & math (unchanged)
from solar_system_viewer import (
    SolarSystemViewer,
    elliptical_orbit_points,
    planet_pos_au,
    GAUSSIAN_GRAVITATIONAL_CONSTANT,  # rad/day
    datetime_to_jd,                    # ISO -> JD helper
)

# ============================ constants & helpers ============================

AU_KM = 1.495978707e8
AU_M  = 1.495978707e11
_MJD_OFFSET = 2400000.5

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
        raise ValueError(f"Missing or non-numeric '{field}' in JSON")
    return v

def _as_poly_xyz(arr_like: Any) -> Optional[np.ndarray]:
    """Convert a list-like to (N,3) float32; guard before np.asarray to avoid Any|None issues."""
    if arr_like is None or not isinstance(arr_like, (list, tuple, np.ndarray)):
        return None
    try:
        arr = np.asarray(arr_like, dtype=float)
    except Exception:
        return None
    if arr.ndim != 2 or arr.shape[1] not in (2, 3):
        return None
    if arr.shape[1] == 2:
        arr = np.c_[arr, np.zeros((arr.shape[0],), dtype=float)]
    if not np.all(np.isfinite(arr)):
        return None
    return arr.astype(np.float32)

def _walk(o: Any) -> Iterable[Dict[str, Any]]:
    if isinstance(o, dict):
        yield o
        for v in o.values():
            yield from _walk(v)
    elif isinstance(o, list):
        for it in o:
            yield from _walk(it)

def _load_json(path: str) -> Optional[Dict[str, Any]]:
    p = Path(path)
    if not p.exists():
        return None
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else None
    except Exception:
        return None

# ====================== intercept / latest readers ======================

def _find_poly_and_times(doc: Dict[str, Any]) -> Tuple[Optional[np.ndarray], Optional[float], Optional[float]]:
    # Prefer lambert poly keys
    for d in _walk(doc):
        for k, v in d.items():
            name = str(k).lower()
            if "lambert" in name and ("poly" in name or "path" in name):
                poly = _as_poly_xyz(v)
                dep = _to_f(d.get("departure_jd") or d.get("depart_jd"))
                arr = _to_f(d.get("arrival_jd") or d.get("arrive_jd"))
                return poly, dep, arr
    # Last resort: any Nx2/Nx3 numeric array
    for d in _walk(doc):
        for _, v in d.items():
            poly = _as_poly_xyz(v)
            if poly is not None:
                return poly, None, None
    return None, None, None

def _find_target_elements_raw(doc: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    for d in _walk(doc):
        for key in ("target_elements", "tgt_elements", "elements"):
            val = d.get(key)
            if isinstance(val, dict) and ("a" in val or "a_AU" in val):
                return val
    return None

# ---- epoch/ma resolution ----

def _parse_epoch_any(src: Dict[str, Any]) -> Optional[float]:
    # direct JD variants
    for k in ("epoch_jd", "epoch", "epoch_tdb_jd", "epoch_tt_jd", "epoch_jd_tdb", "epoch_jd_tt", "t_epoch", "t0_jd", "M_epoch_jd"):
        v = _to_f(src.get(k))
        if v is not None:
            return v
    # MJD variants
    for k in ("epoch_mjd", "epoch_mjd_tdb", "epoch_mjd_tt", "t0_mjd"):
        v = _to_f(src.get(k))
        if v is not None:
            return v + _MJD_OFFSET
    # ISO strings
    for k in ("epoch_iso", "epoch_str", "epoch_date"):
        s = src.get(k)
        if isinstance(s, str) and s:
            try:
                iso = s.rstrip("Z")
                from datetime import datetime, timezone
                dt = datetime.fromisoformat(iso)
                if dt.tzinfo is None:
                    dt = dt.replace(tzinfo=timezone.utc)
                return datetime_to_jd(dt)
            except Exception:
                pass
    return None

def _parse_ma_any(src: Dict[str, Any]) -> Optional[float]:
    for k in ("ma", "ma_deg", "M", "M_deg", "mean_anomaly", "mean_anomaly_deg"):
        v = _to_f(src.get(k))
        if v is not None:
            return v
    return None

def _parse_tp_any(src: Dict[str, Any]) -> Optional[float]:
    # Perihelion time
    for k in ("tp", "tp_jd", "t_peri", "t_peri_jd", "t_perihelion_jd"):
        v = _to_f(src.get(k))
        if v is not None:
            return v
    for k in ("tp_mjd", "t_peri_mjd"):
        v = _to_f(src.get(k))
        if v is not None:
            return v + _MJD_OFFSET
    return None

def _normalize_elements(src: Dict[str, Any]) -> Optional[Dict[str, float]]:
    """Normalize assorted SBDB-like keys into the baseline format."""
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

    ma = _parse_ma_any(src)
    if ma is not None:
        out["ma"] = ma
    ep = _parse_epoch_any(src)
    if ep is not None:
        out["epoch"] = ep
    tp = _parse_tp_any(src)
    if tp is not None:
        out["_tp"] = tp
    return out

def _enrich_ma_epoch(elems: Dict[str, float], latest_path: str) -> None:
    """
    Fill missing epoch/ma if possible:
      • if ma present, epoch missing, and tp present: epoch = tp + M/n
      • else merge both from latest.json if available
    """
    if ("ma" in elems) and ("epoch" not in elems) and ("_tp" in elems):
        a = elems["a"]
        n = GAUSSIAN_GRAVITATIONAL_CONSTANT / (a ** 1.5)  # rad/day
        M = math.radians(elems["ma"])
        elems["epoch"] = elems["_tp"] + (M / n)
        return

    if ("ma" in elems) and ("epoch" in elems):
        return

    doc = _load_json(latest_path)
    if not doc:
        return
    objs = doc.get("objects")
    if not isinstance(objs, list):
        return
    for entry in objs:
        if not isinstance(entry, dict):
            continue
        raw = entry.get("elements")
        if not isinstance(raw, dict):
            continue
        norm = _normalize_elements(raw)
        if not norm:
            continue
        if "ma" in norm and "epoch" in norm:
            elems.setdefault("ma", norm["ma"])
            elems.setdefault("epoch", norm["epoch"])
            return

def load_intercept_bits(path: str) -> Tuple[Optional[np.ndarray], Optional[float], Optional[float], Optional[Dict[str, float]], Optional[Dict[str, Any]]]:
    doc = _load_json(path)
    if not doc:
        return None, None, None, None, None
    poly, dep, arr = _find_poly_and_times(doc)
    raw = _find_target_elements_raw(doc)
    elems = _normalize_elements(raw) if raw else None
    return poly, dep, arr, elems, raw

# ============================ units detection ============================

LambertUnits = Literal["auto", "au", "km", "m"]

def _detect_units(poly: np.ndarray) -> Tuple[str, float, float]:
    """
    Guess units of the polyline by its median distance from origin.
    Returns (unit_str, scale_to_AU, median_radius_in_original_units).
    """
    p = np.asarray(poly, float)
    r = np.linalg.norm(p, axis=1)
    finite = r[np.isfinite(r)]
    med_r = float(np.median(finite)) if finite.size else 0.0
    # Heuristics: ~1 AU ~ O(1-40) in AU; ~1.5e8 in km; ~1.5e11 in m
    if med_r > 1e10:
        return "m", 1.0 / AU_M, med_r
    if med_r > 1e6:
        return "km", 1.0 / AU_KM, med_r
    return "au", 1.0, med_r

def _convert_poly_to_au(poly: np.ndarray, units: LambertUnits = "auto") -> Tuple[np.ndarray, str, float]:
    """
    Convert polyline to AU based on declared/detected units.
    Returns (poly_au, unit_used, med_radius_original_units).
    """
    p = np.asarray(poly, float)
    if units == "au":
        r = np.linalg.norm(p, axis=1)
        med = float(np.median(r[np.isfinite(r)])) if p.size else 0.0
        return p.astype(np.float32), "au", med
    if units == "km":
        r = np.linalg.norm(p, axis=1)
        med = float(np.median(r[np.isfinite(r)])) if p.size else 0.0
        return (p * (1.0 / AU_KM)).astype(np.float32), "km", med
    if units == "m":
        r = np.linalg.norm(p, axis=1)
        med = float(np.median(r[np.isfinite(r)])) if p.size else 0.0
        return (p * (1.0 / AU_M)).astype(np.float32), "m", med
    # auto
    unit_guess, scale, med = _detect_units(p)
    return (p * scale).astype(np.float32), unit_guess, med

# ============================ polyline sanitizer ============================

def sanitize_polyline(poly_au: np.ndarray, *, k: float = 8.0, max_step_au: Optional[float] = None) -> np.ndarray:
    """
    Adaptive cleaner for (N,3) polylines (in AU):
      - drop non-finite rows
      - threshold = max(max_step_au, k*median_step) if provided
      - cut at large hops and keep the longest contiguous segment
      - drop near-duplicate points
      - FAIL-SAFE: if result < 8 pts, return a 'soft' cleaned version (no cutting)
    """
    p = np.asarray(poly_au, float)
    if p.ndim != 2 or p.shape[1] != 3:
        return p

    # 1) drop non-finite
    mask = np.all(np.isfinite(p), axis=1)
    p = p[mask]
    if p.shape[0] < 2:
        return p

    # robust step stats in AU
    d = np.linalg.norm(np.diff(p, axis=0), axis=1)
    finite = d[np.isfinite(d)]
    med = float(np.median(finite[finite > 1e-12])) if finite.size else 0.0
    thr_auto = (k * med) if med > 0 else 1e-3  # ~150,000 km
    thr = max_step_au if (max_step_au is not None) else thr_auto

    # 2) cut at hops
    cuts = np.where(d > thr)[0]
    if cuts.size > 0:
        idx = np.r_[0, cuts + 1, p.shape[0]]
        spans = [(idx[i], idx[i + 1]) for i in range(len(idx) - 1)]
        a, b = max(spans, key=lambda ab: ab[1] - ab[0])
        p2 = p[a:b]
    else:
        p2 = p

    # 3) drop near-duplicates
    if p2.shape[0] >= 2:
        d2 = np.linalg.norm(np.diff(p2, axis=0), axis=1)
        keep = np.r_[True, d2 > 1e-9]
        p2 = p2[keep]

    # 4) fail-safe
    if p2.shape[0] < 8:
        d2 = np.linalg.norm(np.diff(p, axis=0), axis=1)
        keep = np.r_[True, d2 > 1e-9]
        p_soft = p[keep]
        print(f"[viewer] sanitize: fallback soft-clean (N {p.shape[0]} → {p_soft.shape[0]}), med={med:.3g} AU, thr={thr:.3g} AU")
        return p_soft

    print(f"[viewer] sanitize: N {p.shape[0]} → {p2.shape[0]} (med={med:.3g} AU, thr={thr:.3g} AU, cuts={cuts.size})")
    return p2

# ============================ viewer subclass ============================

class InterceptSolarSystemViewer(SolarSystemViewer):
    """
    Adds:
      • NEO orbit curve (heliocentric)
      • NEO moving marker (resolves epoch via multiple paths)
      • Optional Lambert arc + probe (units → AU, then sanitized)
      • Optionally hide outer planets (Jupiter..Neptune)
    """

    def __init__(
        self,
        *,
        intercept_path: Optional[str] = "data/hazardous_neos/latest_intercept.json",
        latest_feed_path: Optional[str] = "data/hazardous_neos/latest.json",
        start_at_departure: bool = False,
        show_lambert: bool = True,
        show_neo_orbit: bool = True,
        show_neo_marker: bool = True,
        show_outer: bool = True,
        sanitize: bool = True,
        sanitizer_k: float = 8.0,
        sanitizer_max_step: Optional[float] = None,
        lambert_units: LambertUnits = "auto",
        day_step: float = 1.0,
        refresh_ms: int = 25,
    ) -> None:
        super().__init__(day_step=day_step, refresh_ms=refresh_ms)

        # Optionally hide outer giants to reduce long chords at inner zoom
        if not show_outer:
            keep = {"Mercury", "Venus", "Earth", "Mars"}
            for name in list(self.orbit_items.keys()):
                if name not in keep:
                    item = self.orbit_items.pop(name, None)
                    if item:
                        self.view.removeItem(item)
            for name in list(self.planet_items.keys()):
                if name not in keep:
                    item = self.planet_items.pop(name, None)
                    if item:
                        self.view.removeItem(item)

        self._lambert: Optional[Dict[str, Any]] = None
        self._neo: Optional[Dict[str, Any]] = None
        self._sanitize = sanitize
        self._san_k = float(sanitizer_k)
        self._san_max = sanitizer_max_step  # Optional[float]
        self._units_pref = lambert_units

        # Load intercept bits and elements
        poly_raw, dep, arr, elems, raw = (None, None, None, None, None)
        if intercept_path:
            poly_raw, dep, arr, elems, raw = load_intercept_bits(intercept_path)

        # Enrich epoch/ma if needed
        if elems is not None and latest_feed_path:
            _enrich_ma_epoch(elems, latest_feed_path)

        # --- NEO orbit first ---
        if show_neo_orbit and (elems is not None):
            orbit = elliptical_orbit_points(
                {
                    "a":    _require_float(elems.get("a"),    "a"),
                    "e":    _require_float(elems.get("e"),    "e"),
                    "i":    _require_float(elems.get("i"),    "i"),
                    "raan": _require_float(elems.get("raan"), "raan"),
                    "argp": _require_float(elems.get("argp"), "argp"),
                },
                num_points=512,
            )
            neo_orbit_item = gl.GLLinePlotItem(
                pos=orbit, width=1.2, color=(1.0, 0.35, 0.35, 0.75),
                antialias=True, mode="line_strip"
            )
            self.view.addItem(neo_orbit_item)

            neo_marker = None
            if show_neo_marker and ("ma" in elems) and ("epoch" in elems):
                pos = planet_pos_au(
                    {
                        "a": elems["a"], "e": elems["e"], "i": elems["i"],
                        "raan": elems["raan"], "argp": elems["argp"],
                        "ma": elems["ma"], "epoch": elems["epoch"],
                    },
                    self.sim_jd,
                )
                neo_marker = gl.GLScatterPlotItem(
                    pos=pos.reshape(1, 3), size=9.0, color=(1.0, 0.55, 0.55, 1.0), pxMode=True
                )
                self.view.addItem(neo_marker)
            else:
                if show_neo_marker:
                    print("[viewer] NEO marker skipped: missing ma/epoch after enrichment")

            self._neo = {"elems": elems, "orbit_item": neo_orbit_item, "marker": neo_marker}

        # --- Lambert overlay (optional) ---
        if show_lambert and isinstance(poly_raw, np.ndarray):
            # Convert units → AU before sanitizing/drawing
            poly_au, unit_used, med_r_orig = _convert_poly_to_au(poly_raw, self._units_pref)
            print(f"[viewer] lambert units: {unit_used}→AU (median |r|={med_r_orig:.3g} {unit_used})")

            pdraw = poly_au
            if self._sanitize:
                pdraw = sanitize_polyline(poly_au, k=self._san_k, max_step_au=self._san_max)
            # If sanitizer pruned too much, draw raw so something is visible
            if pdraw.shape[0] < 2:
                print("[viewer] sanitize left <2 points; drawing raw polyline instead")
                pdraw = poly_au

            line = gl.GLLinePlotItem(
                pos=pdraw, width=1.0, color=(1.0, 0.55, 0.15, 1.0),
                antialias=False, mode="line_strip", glOptions="opaque",
            )
            self.view.addItem(line)

            probe = gl.GLScatterPlotItem(
                pos=pdraw[:1], size=10.0, color=(1.0, 0.9, 0.3, 1.0), pxMode=True
            )
            self.view.addItem(probe)

            self._lambert = {
                "poly": cast(np.ndarray, pdraw),
                "n": int(pdraw.shape[0]),
                "depart_jd": dep if isinstance(dep, float) else None,
                "arrive_jd": arr if isinstance(arr, float) else None,
                "line": line,
                "probe": probe,
                "_i": 0,
            }

            if start_at_departure and self._lambert["depart_jd"] is not None:
                self.sim_jd = cast(float, self._lambert["depart_jd"])
                print(f"[viewer] sim_jd set to departure JD {self.sim_jd}")

        print(
            f"[viewer] NEO elems: {bool(elems)} "
            f"| has ma/epoch: {('ma' in elems) and ('epoch' in elems) if elems else False} "
            f"| Lambert poly: {poly_raw.shape if isinstance(poly_raw, np.ndarray) else None}"
        )

    def _tick(self) -> None:
        super()._tick()

        # Move NEO marker (if present)
        if self._neo is not None and self._neo.get("marker") is not None:
            elems = cast(Dict[str, float], self._neo["elems"])
            pos = planet_pos_au(
                {
                    "a": elems['a'], "e": elems['e'], "i": elems['i'],
                    "raan": elems['raan'], "argp": elems['argp'],
                    "ma": elems['ma'], "epoch": elems['epoch'],
                },
                self.sim_jd,
            )
            cast(gl.GLScatterPlotItem, self._neo["marker"]).setData(pos=pos.reshape(1, 3))

        # Move Lambert probe (if shown)
        L = self._lambert
        if L:
            poly: np.ndarray = cast(np.ndarray, L["poly"])
            n: int = cast(int, L["n"])
            dep: Optional[float] = cast(Optional[float], L["depart_jd"])
            arr: Optional[float] = cast(Optional[float], L["arrive_jd"])
            if dep is not None and arr is not None and arr > dep:
                alpha = (self.sim_jd - dep) / (arr - dep)
                alpha = 0.0 if alpha < 0 else (1.0 if alpha > 1.0 else alpha)
                idx = int(round(alpha * (n - 1)))
            else:
                idx = (cast(int, L["_i"]) + 1) % n
                L["_i"] = idx
            cast(gl.GLScatterPlotItem, L["probe"]).setData(pos=poly[idx:idx + 1])


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--intercept", default="data/hazardous_neos/latest_intercept.json",
                    help="Intercept JSON (may include lambert + elements)")
    ap.add_argument("--latest", default="data/hazardous_neos/latest.json",
                    help="Fallback feed to fill ma/epoch if missing")
    ap.add_argument("--start-at-departure", action="store_true")
    ap.add_argument("--only-inner", action="store_true", help="Hide Jupiter..Neptune")
    ap.add_argument("--no-lambert", action="store_true")
    ap.add_argument("--no-neo-orbit", action="store_true")
    ap.add_argument("--no-neo-marker", action="store_true")
    ap.add_argument("--no-sanitize", action="store_true")
    ap.add_argument("--sanitizer-k", type=float, default=8.0,
                    help="Multiplier on median step (AU) for hop threshold")
    ap.add_argument("--sanitizer-max-step-au", type=float, default=None,
                    help="Absolute hop threshold in AU (overrides k*median if larger)")
    ap.add_argument("--lambert-units", choices=["auto", "au", "km", "m"], default="auto",
                    help="Units of lambert polyline in JSON")
    ap.add_argument("--day-step", type=float, default=1.0)
    ap.add_argument("--refresh-ms", type=int, default=25)
    args = ap.parse_args()

    app = QtWidgets.QApplication([])
    viewer = InterceptSolarSystemViewer(
        intercept_path=args.intercept,
        latest_feed_path=args.latest,
        start_at_departure=args.start_at_departure,
        show_lambert=not args.no_lambert,
        show_neo_orbit=not args.no_neo_orbit,
        show_neo_marker=not args.no_neo_marker,
        show_outer=not args.only_inner,
        sanitize=not args.no_sanitize,
        sanitizer_k=args.sanitizer_k,
        sanitizer_max_step=args.sanitizer_max_step_au,
        lambert_units=cast(LambertUnits, args.lambert_units),
        day_step=args.day_step,
        refresh_ms=args.refresh_ms,
    )
    viewer.show()
    return app.exec()

if __name__ == "__main__":
    raise SystemExit(main())
