# scripts/hud_metrics.py
from __future__ import annotations

import math
from typing import Any

import numpy as np

from scripts.ephem_earth import earth_rv
from scripts.orbital import AU_M, MU_SUN, deg2rad, kepler_E_from_M, wrap_2pi

JD_J2000 = 2451545.0
DAY_S = 86400.0
AU_KM = AU_M / 1000.0
MOON_MEAN_DIST_KM = 384_400.0


def _asteroid_M_at_time(el: dict[str, Any], t_jd: float) -> float | None:
    """Compute mean anomaly at time t using tp if available, else ma@epoch."""
    a = float(el["a_AU"]) * AU_M
    n = math.sqrt(MU_SUN / (a**3))  # rad/s

    tp = el.get("tp_jd_tdb")
    if tp is not None:
        return wrap_2pi(n * ((t_jd - float(tp)) * DAY_S))

    # fallback: ma at epoch + n*dt
    ma_deg = el.get("ma_deg")
    epoch = el.get("epoch_jd_tdb")
    if ma_deg is None or epoch is None:
        return None
    M0 = deg2rad(float(ma_deg))
    dt = (t_jd - float(epoch)) * DAY_S
    return wrap_2pi(M0 + n * dt)


def asteroid_rv_from_elements(el: dict[str, Any], jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
    a = float(el["a_AU"]) * AU_M
    e = float(el["e"])
    inc = deg2rad(float(el["i_deg"]))
    raan = deg2rad(float(el["raan_deg"]))
    argp = deg2rad(float(el["argp_deg"]))

    M = _asteroid_M_at_time(el, jd_tdb)
    if M is None:
        raise ValueError("insufficient elements to compute mean anomaly (need tp or ma@epoch)")

    E = kepler_E_from_M(M, e)
    cE, sE = math.cos(E), math.sin(E)
    r_peri = np.array([a * (cE - e), a * (math.sqrt(1 - e * e) * sE), 0.0])
    r_mag = a * (1 - e * cE)
    fac = math.sqrt(MU_SUN * a) / r_mag
    v_peri = np.array([-fac * sE, fac * math.sqrt(1 - e * e) * cE, 0.0])

    cO, sO = math.cos(raan), math.sin(raan)
    ci, si = math.cos(inc), math.sin(inc)
    cw, sw = math.cos(argp), math.sin(argp)
    R = np.array(
        [
            [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
            [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
            [sw * si, cw * si, ci],
        ]
    )
    return R @ r_peri, R @ v_peri


def compute_hud(el: dict[str, Any], jd_tdb: float) -> dict[str, Any]:
    r_ast, _ = asteroid_rv_from_elements(el, jd_tdb)
    r_ear, _ = earth_rv(jd_tdb)

    earth_dist_au = float(np.linalg.norm(r_ast - r_ear) / AU_M)
    sun_dist_au = float(np.linalg.norm(r_ast) / AU_M)

    geo = (r_ast - r_ear) / 1000.0  # km
    geo_xy_km = [float(geo[0]), float(geo[1])]

    return {
        "earth_distance_au": earth_dist_au,
        "earth_distance_km": earth_dist_au * AU_KM,
        "sun_distance_au": sun_dist_au,
        "geo_xy_km": geo_xy_km,
        "moon_radius_km": MOON_MEAN_DIST_KM,
    }
