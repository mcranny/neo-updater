# scripts/ephem_earth.py
from __future__ import annotations
from typing import Tuple, Any
import numpy as np

# Obliquity of the ecliptic at J2000 (IAU 2006/2000A)
_OBLIQUITY_J2000_DEG = 23.439291111

def _eq_to_ecl_rot() -> np.ndarray:
    """Rotation matrix: ICRF/Equatorial -> Ecliptic J2000 (R_x(-ε))."""
    eps = np.deg2rad(_OBLIQUITY_J2000_DEG)
    c, s = float(np.cos(eps)), float(np.sin(eps))
    # R_x(-ε) = [[1,0,0],[0, c, s],[0, -s, c]]
    return np.array([[1.0, 0.0, 0.0],
                     [0.0,   c,   s],
                     [0.0,  -s,   c]], dtype=float)

_R_EQ2ECL = _eq_to_ecl_rot()

def earth_rv_heliocentric_skyfield(jd_tdb: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Return Earth's heliocentric state in **ecliptic J2000**, meters & m/s,
    using Skyfield + JPL DE ephemeris at TDB Julian Day `jd_tdb`.
    """
    from skyfield.api import load

    ts = load.timescale()
    try:
        t = ts.tdb_jd(jd_tdb)          # type: ignore[attr-defined]
    except Exception:
        t = ts.tdb(jd=jd_tdb)          # type: ignore[attr-defined]

    try:
        eph = load("de440s.bsp")
    except Exception:
        eph = load("de421.bsp")

    e: Any = eph["earth"]
    s: Any = eph["sun"]

    # Skyfield states are ICRF/Equatorial. Compute Earth wrt Sun, then rotate.
    rel = e.at(t) - s.at(t)            # type: ignore[attr-defined]

    r_km   = rel.position.km           # type: ignore[attr-defined]
    v_km_s = rel.velocity.km_per_s     # type: ignore[attr-defined]

    r_m   = np.array(r_km, dtype=float)   * 1000.0
    v_m_s = np.array(v_km_s, dtype=float) * 1000.0

    # Equatorial -> Ecliptic J2000
    r_ecl = _R_EQ2ECL @ r_m
    v_ecl = _R_EQ2ECL @ v_m_s
    return r_ecl, v_ecl
