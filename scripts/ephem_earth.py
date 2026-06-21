"""Earth heliocentric state providers in the ecliptic J2000 frame."""

from __future__ import annotations

import math
import os
from pathlib import Path

import numpy as np

from scripts.orbital import AU_M, MU_SUN, deg2rad, kepler_E_from_M, oe_to_rv, wrap_2pi

JD_J2000 = 2451545.0
DAY_S = 86400.0
_OBLIQUITY_J2000_DEG = 23.439291111


def earth_rv_heliocentric_analytic(jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
    """Low-precision J2000 Kepler model, suitable for deterministic sizing."""
    semi_major_axis = 1.00000011 * AU_M
    eccentricity = 0.01671022
    inclination = deg2rad(0.00005)
    longitude_ascending_node = deg2rad(-11.26064)
    longitude_perihelion = deg2rad(102.94719)
    argument_perihelion = longitude_perihelion - longitude_ascending_node
    mean_longitude_j2000 = deg2rad(100.46435)
    mean_anomaly_j2000 = mean_longitude_j2000 - longitude_perihelion
    mean_motion = math.sqrt(MU_SUN / semi_major_axis**3)
    mean_anomaly = wrap_2pi(mean_anomaly_j2000 + mean_motion * (jd_tdb - JD_J2000) * DAY_S)
    eccentric_anomaly = kepler_E_from_M(mean_anomaly, eccentricity)
    true_anomaly = 2.0 * math.atan2(
        math.sqrt(1.0 + eccentricity) * math.sin(eccentric_anomaly / 2.0),
        math.sqrt(1.0 - eccentricity) * math.cos(eccentric_anomaly / 2.0),
    )
    return oe_to_rv(
        semi_major_axis,
        eccentricity,
        inclination,
        longitude_ascending_node,
        argument_perihelion,
        true_anomaly,
    )


def earth_rv_heliocentric_skyfield(jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
    """JPL DE state through Skyfield, rotated from ICRF to ecliptic J2000."""
    from skyfield.api import load

    timescale = load.timescale()
    time = timescale.tdb_jd(jd_tdb)
    configured_path = os.getenv("JPL_EPHEMERIS_PATH", "").strip()
    ephemeris = (
        load(str(Path(configured_path).expanduser())) if configured_path else load("de440s.bsp")
    )
    relative = ephemeris["earth"].at(time) - ephemeris["sun"].at(time)
    position_m = np.asarray(relative.position.km, dtype=float) * 1000.0
    velocity_mps = np.asarray(relative.velocity.km_per_s, dtype=float) * 1000.0
    epsilon = deg2rad(_OBLIQUITY_J2000_DEG)
    cosine, sine = math.cos(epsilon), math.sin(epsilon)
    equatorial_to_ecliptic = np.array([[1.0, 0.0, 0.0], [0.0, cosine, sine], [0.0, -sine, cosine]])
    return equatorial_to_ecliptic @ position_m, equatorial_to_ecliptic @ velocity_mps


def earth_rv(jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
    """Select the deterministic analytic model or optional JPL DE ephemeris."""
    mode = os.getenv("EPHEMERIS_MODE", "analytic").strip().lower()
    if mode == "skyfield":
        return earth_rv_heliocentric_skyfield(jd_tdb)
    if mode != "analytic":
        raise ValueError("EPHEMERIS_MODE must be 'analytic' or 'skyfield'")
    return earth_rv_heliocentric_analytic(jd_tdb)
