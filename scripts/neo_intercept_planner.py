#!/usr/bin/env python3
"""
NEO Intercept Planner (fixed version)
=====================================

This script is adapted from the original neo-updater `neo_intercept_planner.py`.
It includes additional safeguards to ensure that Lambert transfers are
calculated robustly for each near-Earth object (NEO) in the input.

Key modifications:
- Bounds the time of flight (TOF) for each target to 50 % of its orbital period
  computed from its semi-major axis.  This avoids multi-revolution transfers.
- Retries the Lambert solver with progressively shorter TOFs if the first attempt
  fails; after three retries the target is skipped.
- Wraps the Lambert solver in a `safe_lambert_polyline` function to catch and
  suppress exceptions, preventing the entire script from crashing.
- Uses overflow-safe Stumpff functions in both the Lambert solver and the
  universal Kepler propagator.

Run this script as you would the original planner.  Use `--debug` to see
diagnostic messages about TOF bounds and retries.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
import re
import hashlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# -----------------------------------------------------------------------------
# Constants
AU_M = 149_597_870_700.0
DEG = math.pi / 180.0
MU_SUN = 1.32712440018e20  # m^3/s^2

# -----------------------------------------------------------------------------
# Time helpers
def jd_from_unix(t: float) -> float:
    return t / 86400.0 + 2440587.5

def unix_from_jd(jd: float) -> float:
    return (jd - 2440587.5) * 86400.0

def utc_from_jd(jd: float) -> str:
    import datetime as _dt
    return _dt.datetime.fromtimestamp(unix_from_jd(jd), tz=_dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S")

# -----------------------------------------------------------------------------
# SBDB fetch (retries + cache)
def _sbdb_fetch(des_or_name: str, *, debug: bool = False) -> Optional[Dict[str, Any]]:
    """Fetch NEO orbital data from JPL's SBDB (with caching and retries)."""
    import requests
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    cache_dir = Path("data/sbdb_cache")
    cache_dir.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha1(des_or_name.encode("utf-8")).hexdigest()[:16]
    cpath = cache_dir / f"{key}.json"
    try:
        # Return cached copy if fresh (<7 days old)
        if cpath.exists() and (time.time() - cpath.stat().st_mtime) < 7 * 86400:
            if debug:
                print(f"[sbdb] cache hit for {des_or_name!r}")
            return json.loads(cpath.read_text(encoding="utf-8"))
    except Exception:
        pass
    # Query SBDB API
    url = "https://ssd-api.jpl.nasa.gov/sbdb.api"
    params = {"sstr": des_or_name, "full-prec": "true", "soln-epoch": "true"}
    sess = requests.Session()
    retry = Retry(total=4, backoff_factor=0.6, status_forcelist=[429, 500, 502, 503, 504],
                  allowed_methods=["GET"], raise_on_status=False)
    sess.mount("https://", HTTPAdapter(max_retries=retry))
    try:
        r = sess.get(url, params=params, timeout=20)
    except Exception as e:
        if debug:
            print(f"[sbdb] request error for {des_or_name!r}: {e}")
        return None
    if r.status_code != 200:
        if debug:
            print(f"[sbdb] HTTP {r.status_code} for {des_or_name!r}")
        return None
    data = r.json()
    if not isinstance(data, dict) or "orbit" not in data:
        if debug:
            print(f"[sbdb] no 'orbit' in payload for {des_or_name!r}")
        return None
    try:
        cpath.write_text(json.dumps(data), encoding="utf-8")
    except Exception:
        pass
    return data

# -----------------------------------------------------------------------------
# Parsing helpers
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
    """Flatten element entries to a dictionary."""
    if isinstance(el, dict):
        return el
    if isinstance(el, list):
        out = {}
        for d in el:
            try:
                out[str(d.get("name"))] = d.get("value")
            except Exception:
                pass
        return out
    return {}

def _parse_elements(sbdb: Dict[str, Any], *, debug: bool = False) -> Optional[Dict[str, float]]:
    """Extract orbital elements from SBDB payload and return a dict."""
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
        a_AU = g("a")
        e = g("e")
        i_deg = g("i")
        raan = g("om", "Omega", "raan")
        argp = g("w", "argp", "peri")
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
                ma_deg = ((orbit.get("epoch_jd") or epoch_jd) - tp_jd) * n_deg_day
        return {
            "a_AU": a_AU,
            "e": e,
            "i_deg": i_deg,
            "raan_deg": raan,
            "argp_deg": argp,
            "ma_deg": ma_deg,
            "epoch_jd": epoch_jd,
        }
    except Exception:
        if debug:
            print("[planner] element parsing failed")
        return None

# -----------------------------------------------------------------------------
# Orbital propagation helpers
def kepler_E_from_M(M: float, e: float, tol: float = 1e-12, itmax: int = 60) -> float:
    """Solve Kepler’s equation for the eccentric anomaly."""
    M = float(M) % (2.0 * math.pi)
    if e < 1e-12:
        return M
    E = M if e < 0.8 else math.pi
    for _ in range(itmax):
        f = E - e * math.sin(E) - M
        fp = 1.0 - e * math.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            break
    return E

def rv_from_elements_at_time(a_AU: float, e: float, i_deg: float, raan_deg: float, argp_deg: float,
                             ma0_deg: float, epoch0_jd: float, t_jd: float, mu: float = MU_SUN
                             ) -> Tuple[np.ndarray, np.ndarray]:
    """Compute position and velocity at time t_jd given orbital elements."""
    a = a_AU * AU_M
    i = i_deg * DEG
    Om = raan_deg * DEG
    w = argp_deg * DEG
    M0 = ma0_deg * DEG
    n = math.sqrt(mu / (a ** 3))
    dt = (t_jd - epoch0_jd) * 86400.0
    M = M0 + n * dt
    E = kepler_E_from_M(M, e)
    cosE = math.cos(E)
    sinE = math.sin(E)
    fac = math.sqrt(1.0 - e * e)
    r_p = a * (1.0 - e * cosE)
    x_p = a * (cosE - e)
    y_p = a * fac * sinE
    rdot = (math.sqrt(mu * a) / r_p) * e * sinE
    rfdot = (math.sqrt(mu * a) / r_p) * fac * cosE
    cw = math.cos(w); sw = math.sin(w)
    cO = math.cos(Om); sO = math.sin(Om)
    ci = math.cos(i); si = math.sin(i)
    # Rotation matrix from perifocal to inertial
    R11 = cO * cw - sO * sw * ci
    R12 = -cO * sw - sO * cw * ci
    R13 = sO * si
    R21 = sO * cw + cO * sw * ci
    R22 = -sO * sw + cO * cw * ci
    R23 = -cO * si
    R31 = sw * si
    R32 = cw * si
    R33 = ci
    R = np.array([[R11, R12, R13], [R21, R22, R23], [R31, R32, R33]], dtype=float)
    r_pf = np.array([x_p, y_p, 0.0])
    v_pf = np.array([rdot, rfdot, 0.0])
    r = R @ r_pf
    v = R @ v_pf
    return r, v

# -----------------------------------------------------------------------------
# Lambert solver with overflow-safe Stumpff functions
def lambert_universal(r1, r2, tof_s: float, mu: float = MU_SUN, prograde: bool = True,
                      max_iter: int = 80, tol: float = 1e-8) -> Tuple[np.ndarray, np.ndarray]:
    """Lambert via universal variables (single revolution).  Returns (v1, v2)."""
    r1 = np.asarray(r1, dtype=float).reshape(3)
    r2 = np.asarray(r2, dtype=float).reshape(3)
    r1n = float(np.linalg.norm(r1))
    r2n = float(np.linalg.norm(r2))
    if r1n == 0.0 or r2n == 0.0:
        raise ValueError("Lambert: |r1| or |r2| is zero")
    sqrt_mu = math.sqrt(mu)
    # geometry
    cosd = float(np.dot(r1, r2) / (r1n * r2n + 1e-16))
    cosd = max(-1.0, min(1.0, cosd))
    dtheta = math.acos(cosd)
    # short- or long-way selection using out-of-plane direction
    cz = float(np.cross(r1, r2)[2])
    if prograde and cz < 0.0:
        dtheta = 2.0 * math.pi - dtheta
    if not prograde and cz > 0.0:
        dtheta = 2.0 * math.pi - dtheta
    denom = 1.0 - cosd
    if abs(denom) < 1e-16:
        raise ValueError("Lambert geometry singular (Δθ ~ 0)")
    A = math.sin(dtheta) * math.sqrt(max(0.0, r1n * r2n / denom))
    if abs(A) < 1e-14:
        raise ValueError("Lambert geometry singular: A≈0")
    # define y(z) and tof(z) with overflow-safe stumpff functions
    def stumpff_C(z: float) -> float:
        z = float(z); az = abs(z)
        if az < 1e-8:
            return 0.5 - z / 24.0 + (z * z) / 720.0
        if z > 0.0:
            s = math.sqrt(z)
            return (1.0 - math.cos(s)) / z
        s = math.sqrt(-z)
        if s > 50:
            # cosh(s) ~ 0.5 e^s for large s
            return (0.5 * math.exp(s) - 1.0) / (-z)
        return (math.cosh(s) - 1.0) / (-z)
    def stumpff_S(z: float) -> float:
        z = float(z); az = abs(z)
        if az < 1e-8:
            return (1.0 / 6.0) - z / 120.0 + (z * z) / 5040.0
        if z > 0.0:
            s = math.sqrt(z)
            return (s - math.sin(s)) / (s ** 3)
        s = math.sqrt(-z)
        if s > 50:
            # sinh(s) ~ 0.5 e^s for large s
            return (0.5 * math.exp(s) - s) / (s ** 3)
        return (math.sinh(s) - s) / (s ** 3)
    def y(z: float) -> float:
        C = stumpff_C(z)
        S = stumpff_S(z)
        Ceff = C if abs(C) > 1e-16 else 1e-16
        return r1n + r2n + A * ((z * S - 1.0) / Ceff)
    def tof(z: float) -> float:
        C = stumpff_C(z)
        S = stumpff_S(z)
        yy = y(z)
        if yy <= 0.0:
            return float('inf')
        Ceff = C if abs(C) > 1e-16 else 1e-16
        x = math.sqrt(max(0.0, yy / Ceff))
        return (x ** 3 * S + A * math.sqrt(yy)) / sqrt_mu
    # bracket and bisection on z
    z_lo, z_hi = -4.0 * math.pi ** 2, 4.0 * math.pi ** 2
    while y(z_lo) <= 0.0:
        z_lo += 0.5
        if z_lo > 0.0:
            break
    t_hi = tof(z_hi)
    tries = 0
    while (not math.isfinite(t_hi)) or (t_hi < tof_s):
        z_hi *= 2.0
        tries += 1
        t_hi = tof(z_hi)
        if tries > 40:
            break
    z = 0.0
    for _ in range(max_iter):
        t_z = tof(z)
        if abs(t_z - tof_s) < tol:
            break
        if t_z < tof_s:
            z_lo = z
        else:
            z_hi = z
        z = 0.5 * (z_lo + z_hi)
    yy = y(z)
    if yy <= 0.0:
        raise RuntimeError("Lambert: y(z*) <= 0 at solution")
    f = 1.0 - yy / r1n
    g = A * math.sqrt(yy / mu)
    gdot = 1.0 - yy / r2n
    if abs(g) < 1e-14:
        g = math.copysign(1e-14, g)
    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g
    return v1, v2

def kepler_universal_propagate(r0, v0, dt: float, mu: float = MU_SUN, itmax: int = 60, tol: float = 1e-9
                               ) -> Tuple[np.ndarray, np.ndarray]:
    """Propagate (r0, v0) by dt seconds via universal variables with overflow-safe functions."""
    r0 = np.asarray(r0, dtype=float)
    v0 = np.asarray(v0, dtype=float)
    r0n = float(np.linalg.norm(r0))
    if r0n == 0.0:
        raise ValueError("kepler_universal_propagate: |r0| is zero")
    dt = float(dt)
    vr0 = float(np.dot(r0, v0) / r0n)
    alpha = 2.0 / r0n - float(np.dot(v0, v0)) / mu
    sqrt_mu = math.sqrt(mu)
    # initial guess for chi
    if abs(alpha) > 1e-12:
        chi = math.copysign(1.0, dt) * math.sqrt(mu) * abs(alpha) * dt
    else:
        h = float(np.linalg.norm(np.cross(r0, v0)))
        p = h * h / mu
        s = 0.5 * (math.pi / 2.0 - math.atan(3.0 * math.sqrt(mu / (p ** 3)) * dt))
        w = math.atan(math.tan(s) ** (1.0 / 3.0))
        chi = math.sqrt(p) * 2.0 * math.tan(w)
    for _ in range(itmax):
        z = alpha * chi * chi
        # overflow-safe stumpff functions
        if abs(z) < 1e-8:
            Cz = 0.5 - z / 24.0 + (z * z) / 720.0
            Sz = (1.0 / 6.0) - z / 120.0 + (z * z) / 5040.0
        elif z > 0.0:
            s = math.sqrt(z)
            Cz = (1.0 - math.cos(s)) / z
            Sz = (s - math.sin(s)) / (s ** 3)
        else:
            s = math.sqrt(-z)
            if s > 50:
                Cz = (0.5 * math.exp(s) - 1.0) / (-z)
                Sz = (0.5 * math.exp(s) - s) / (s ** 3)
            else:
                Cz = (math.cosh(s) - 1.0) / (-z)
                Sz = (math.sinh(s) - s) / (s ** 3)
        r = chi * chi * Cz + vr0 * chi * (1.0 - z * Sz) + r0n * (1.0 - z * Cz)
        F = r - sqrt_mu * dt
        if abs(F) < tol:
            break
        dF = chi * (1.0 - z * Sz)
        chi -= F / (dF + 1e-16)
    z = alpha * chi * chi
    if abs(z) < 1e-8:
        Cz = 0.5 - z / 24.0 + (z * z) / 720.0
        Sz = (1.0 / 6.0) - z / 120.0 + (z * z) / 5040.0
    elif z > 0.0:
        s = math.sqrt(z)
        Cz = (1.0 - math.cos(s)) / z
        Sz = (s - math.sin(s)) / (s ** 3)
    else:
        s = math.sqrt(-z)
        if s > 50:
            Cz = (0.5 * math.exp(s) - 1.0) / (-z)
            Sz = (0.5 * math.exp(s) - s) / (s ** 3)
        else:
            Cz = (math.cosh(s) - 1.0) / (-z)
            Sz = (math.sinh(s) - s) / (s ** 3)
    f = 1.0 - (chi * chi * Cz) / r0n
    g = dt - (chi ** 3) * Sz / sqrt_mu
    r = f * r0 + g * v0
    rn = float(np.linalg.norm(r))
    fdot = (sqrt_mu / (rn * r0n)) * (z * Sz - 1.0) * chi
    gdot = 1.0 - (chi * chi * Cz) / rn
    v = fdot * r0 + gdot * v0
    return r, v

def _force_endpoints(poly_m: np.ndarray, r1_m: np.ndarray, r2_m: np.ndarray) -> np.ndarray:
    """Ensure the first and last points of the transfer polyline match the endpoints exactly."""
    poly_m = np.asarray(poly_m, dtype=float)
    poly_m[0, :] = np.asarray(r1_m, dtype=float).reshape(3)
    poly_m[-1, :] = np.asarray(r2_m, dtype=float).reshape(3)
    return poly_m

def lambert_polyline_always(r1_m: np.ndarray, r2_m: np.ndarray, tof_s: float, n: int) -> np.ndarray:
    """
    Compute a polyline for the transfer orbit between r1 and r2 over time tof_s.

    Tries the universal-variable Lambert solver (prograde then retrograde).  If
    both calls fail, falls back to a shooting solver.  Returns an (N, 3) array.
    """
    r1_m = np.asarray(r1_m, dtype=float).reshape(3)
    r2_m = np.asarray(r2_m, dtype=float).reshape(3)
    tof_s = float(tof_s)
    n = int(max(2, n))
    # Try universal-variable Lambert for both directions
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
    # Fallback: shooting solve in v1 space (Broyden-like with finite differences)
    chord = (r2_m - r1_m) / max(1e-9, tof_s)
    zhat = np.cross(r1_m, r2_m)
    t_hat = np.cross(zhat, r1_m)
    t_n = np.linalg.norm(t_hat)
    if t_n > 0:
        t_hat /= t_n
    v_circ = math.sqrt(MU_SUN / max(1e-3, np.linalg.norm(r1_m))) * t_hat
    v = 0.6 * chord + 0.8 * v_circ
    def residual(v1: np.ndarray) -> np.ndarray:
        r, _ = kepler_universal_propagate(r1_m, v1, tof_s, MU_SUN)
        return (r - r2_m).astype(float)
    eps = 1e-3
    for _ in range(30):
        F = residual(v)
        if np.linalg.norm(F) < 1.0:  # 1 m tolerance
            break
        J = np.zeros((3, 3), dtype=float)
        for k in range(3):
            dv = np.zeros(3)
            dv[k] = eps
            J[:, k] = (residual(v + dv) - F) / eps
        # Solve J * dv = -F
        try:
            dv, *_ = np.linalg.lstsq(J, -F, rcond=None)
        except Exception:
            dv = -F * 1e-6
        v += 0.8 * dv
    # Sample transfer
    ts = np.linspace(0.0, tof_s, n)
    pts = np.empty((n, 3), dtype=float)
    for i, t in enumerate(ts):
        r, _ = kepler_universal_propagate(r1_m, v, float(t), MU_SUN)
        pts[i, :] = r
    return _force_endpoints(pts, r1_m, r2_m)

def safe_lambert_polyline(r1_m: np.ndarray, r2_m: np.ndarray, tof_s: float, n: int, *, debug: bool = False
                          ) -> Optional[np.ndarray]:
    """Wrapper around lambert_polyline_always that catches exceptions and returns None."""
    try:
        return lambert_polyline_always(r1_m, r2_m, tof_s, n)
    except Exception as e:
        if debug:
            print(f"[planner] Lambert failed: {e}")
        return None

# -----------------------------------------------------------------------------
# Plan builder
# Reference Earth elements (J2000, approximate)
EARTH_ELEMENTS_J2000 = {
    "a_AU": 1.00000011,
    "e": 0.01671022,
    "i_deg": 0.00005,
    "raan_deg": -11.26064,
    "argp_deg": 102.94719,
    "ma_deg": 100.46435,
    "epoch_jd": 2451545.0,
}

def compute_orbital_period_days(a_AU: float) -> float:
    """Return approximate orbital period in days for semi-major axis a_AU (heliocentric)."""
    return 365.25 * (a_AU ** 1.5)

def build_plan_for_one(row: Dict[str, Any], depart_jd: float, base_tof_days: float, poly_n: int, *, debug: bool = False
                       ) -> Optional[Dict[str, Any]]:
    """Build an intercept plan for a single NEO row.  Returns None on failure."""
    # Fetch SBDB data and parse elements
    sb = _sbdb_fetch_best(row, debug=debug)
    if sb is None:
        if debug:
            print("[planner] SBDB fetch failed")
        return None
    el_tgt = _parse_elements(sb, debug=debug)
    if el_tgt is None:
        if debug:
            print("[planner] element parsing failed")
        return None
    # Earth at departure
    el_earth = dict(EARTH_ELEMENTS_J2000)
    r1_m, _ = rv_from_elements_at_time(el_earth["a_AU"], el_earth["e"], el_earth["i_deg"],
                                       el_earth["raan_deg"], el_earth["argp_deg"], el_earth["ma_deg"],
                                       el_earth["epoch_jd"], depart_jd, mu=MU_SUN)
    # Compute a bounded time-of-flight based on target orbital period
    tgt_period_days = compute_orbital_period_days(max(0.01, el_tgt["a_AU"]))
    max_tof_days = 0.5 * tgt_period_days
    tof_days = min(base_tof_days, max_tof_days)
    if debug:
        print(f"[planner] TOF bounded: requested={base_tof_days:.1f} d, max={max_tof_days:.1f} d, using={tof_days:.1f} d")
    arrive_jd = depart_jd + tof_days
    # Target at arrival
    r2_m, _ = rv_from_elements_at_time(el_tgt["a_AU"], el_tgt["e"], el_tgt["i_deg"],
                                       el_tgt["raan_deg"], el_tgt["argp_deg"], el_tgt["ma_deg"],
                                       el_tgt["epoch_jd"], arrive_jd, mu=MU_SUN)
    # Attempt Lambert with retries (halve TOF on failure)
    tof_s = float(tof_days) * 86400.0
    poly_m: Optional[np.ndarray] = None
    attempts = 0
    curr_tof_s = tof_s
    while poly_m is None and attempts < 3:
        poly_m = safe_lambert_polyline(r1_m, r2_m, curr_tof_s, n=int(poly_n), debug=debug)
        if poly_m is None:
            curr_tof_s *= 0.5
            attempts += 1
            if debug:
                print(f"[planner] Lambert retry {attempts}: reducing TOF to {curr_tof_s/86400.0:.1f} d")
    if poly_m is None:
        if debug:
            print("[planner] Failed to compute Lambert arc after retries")
        return None
    poly_au = (poly_m / AU_M).astype(float)
    poly_xy = poly_au[:, :2]
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
        "arrival_utc": utc_from_jd(arrive_jd),
    }
    return plan

def _make_candidates(row: Dict[str, Any]) -> List[str]:
    """Generate candidate identifiers/names to search SBDB."""
    cands: List[str] = []
    for key in ("neo_reference_id", "id", "neo_ref_id"):
        v = (row.get(key) or "").strip()
        if v and v.isdigit():
            cands.append(v)
    for key in ("name", "designation"):
        val = (row.get(key) or "").strip()
        if not val:
            continue
        m = re.search(r"(\d{4,})", val)
        if m:
            num = m.group(1)
            if num not in cands:
                cands.insert(0, num)
        clean = val.replace("(", "").replace(")", "").strip()
        for s in (clean, clean.replace(" ", ""), val):
            if s and s not in cands:
                cands.append(s)
    seen = set(); out = []
    for s in cands:
        if s and s not in seen:
            seen.add(s); out.append(s)
    return out

def _sbdb_fetch_best(row: Dict[str, Any], *, debug: bool = False) -> Optional[Dict[str, Any]]:
    """Try multiple candidate names/IDs until SBDB returns valid data."""
    cands = _make_candidates(row)
    if debug:
        print(f"[sbdb] candidates: {cands}")
    for s in cands:
        sb = _sbdb_fetch(s, debug=debug)
        if sb is not None:
            if debug:
                print(f"[sbdb] success: {s!r}")
            return sb
        else:
            if debug:
                print(f"[sbdb] failed: {s!r}")
    if debug:
        print("[sbdb] all candidates failed")
    return None

def main() -> int:
    ap = argparse.ArgumentParser(description="NEO intercept planner (Lambert, fixed version)")
    ap.add_argument("--inp", default="data/hazardous_neos/latest.json")
    ap.add_argument("--out", default="data/hazardous_neos/latest_intercept.json")
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
    now_jd = jd_from_unix(time.time())
    depart_jd = now_jd + float(args.depart_days)
    updated = 0
    for row in rows:
        plan: Dict[str, Any] = row.setdefault("intercept_plan", {})
        if not args.with_sbdb:
            # Offline: synthesise target elements but still build a curved arc
            el_earth = dict(EARTH_ELEMENTS_J2000)
            el_tgt = plan.get("elements_target") or {
                "a_AU": 1.2, "e": 0.10, "i_deg": 5.0,
                "raan_deg": 50.0, "argp_deg": 10.0,
                "ma_deg": 0.0, "epoch_jd": 2451545.0,
            }
            # Bound TOF using target period
            tgt_period = compute_orbital_period_days(max(0.01, el_tgt["a_AU"]))
            max_tof = 0.5 * tgt_period
            tof_days = min(float(args.tof_days), max_tof)
            if args.debug:
                print(f"[planner] Offline: TOF bounded to {tof_days:.1f} d")
            arrive_jd = depart_jd + tof_days
            r1_m, _ = rv_from_elements_at_time(el_earth["a_AU"], el_earth["e"], el_earth["i_deg"],
                                               el_earth["raan_deg"], el_earth["argp_deg"], el_earth["ma_deg"],
                                               el_earth["epoch_jd"], depart_jd, mu=MU_SUN)
            r2_m, _ = rv_from_elements_at_time(el_tgt["a_AU"], el_tgt["e"], el_tgt["i_deg"],
                                               el_tgt["raan_deg"], el_tgt["argp_deg"], el_tgt["ma_deg"],
                                               el_tgt["epoch_jd"], arrive_jd, mu=MU_SUN)
            tof_s = float(tof_days) * 86400.0
            poly_m = safe_lambert_polyline(r1_m, r2_m, tof_s, n=int(args.polyline_n), debug=args.debug)
            if poly_m is None:
                if args.debug:
                    print("[planner] Offline Lambert failed; skipping")
                continue
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
                "tof_days": float(tof_days),
                "departure_utc": utc_from_jd(depart_jd),
                "arrival_utc": utc_from_jd(arrive_jd),
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
