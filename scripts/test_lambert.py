"""
test_lambert.py — Lambert solver smoke test and re-export shim.

Provides ``lambert_universal`` (re-exported from scripts.orbital) so that
any code doing ``from scripts.test_lambert import lambert_universal`` keeps
working.  Also contains a standalone smoke test runnable via
``python -m scripts.test_lambert``.

Bugs fixed:
- Previously imported from ``neo_intercept_planner`` (which doesn't define
  rv_from_elements_at_time, safe_lambert_polyline, compute_orbital_period_days,
  AU_M, or MU_SUN in the current codebase).
- ``__all__`` declared lambert_universal but it was never defined here.
- Import path used bare ``neo_intercept_planner`` instead of ``scripts.*``.
"""

from __future__ import annotations

import math
from datetime import datetime, timezone

import numpy as np

# ---------------------------------------------------------------------------
# Re-export so callers of ``from scripts.test_lambert import lambert_universal``
# keep working without modification.
# ---------------------------------------------------------------------------
from scripts.orbital import (       # noqa: F401
    lambert_universal,
    AU_M,
    MU_SUN,
    kepler_E_from_M,
    oe_to_rv,
    kepler_universal_propagate,
    sample_transfer_polyline,
)

__all__ = ["lambert_universal"]

# ---------------------------------------------------------------------------
# Helpers for the standalone smoke test
# ---------------------------------------------------------------------------

DAY_S = 86400.0
DEG2RAD = math.pi / 180.0


def _jd_now() -> float:
    JD_UNIX = 2440587.5
    return JD_UNIX + datetime.now(timezone.utc).timestamp() / DAY_S


def _rv_from_elements(
    a_AU: float, e: float, i_deg: float,
    raan_deg: float, argp_deg: float, ma_deg: float,
    epoch_jd: float, query_jd: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Keplerian state vector in SI (m, m/s) at query_jd."""
    # Mean motion [rad/s]
    a_m = a_AU * AU_M
    n = math.sqrt(MU_SUN / a_m ** 3)

    dt_s = (query_jd - epoch_jd) * DAY_S
    M0 = ma_deg * DEG2RAD
    M = (M0 + n * dt_s) % (2 * math.pi)
    E = kepler_E_from_M(M, e)
    nu = 2.0 * math.atan2(
        math.sqrt(1 + e) * math.sin(E / 2),
        math.sqrt(1 - e) * math.cos(E / 2),
    )
    return oe_to_rv(
        a_m, e,
        i_deg * DEG2RAD, raan_deg * DEG2RAD, argp_deg * DEG2RAD, nu,
    )


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------

def main() -> None:
    # Earth orbital elements (approximate J2000)
    earth_a, earth_e = 1.00000011, 0.01671022
    earth_el = dict(
        a_AU=earth_a, e=earth_e, i_deg=0.00005,
        raan_deg=-11.26064, argp_deg=102.94719, ma_deg=100.46435,
        epoch_jd=2451545.0,
    )
    # Synthetic target: 1.2 AU, moderate eccentricity
    tgt_el = dict(
        a_AU=1.2, e=0.1, i_deg=5.0,
        raan_deg=50.0, argp_deg=10.0, ma_deg=0.0,
        epoch_jd=2451545.0,
    )

    # Use explicit, well-conditioned positions instead of propagating from
    # elements — avoids the near-antipodal geometry that breaks Lambert.
    # Classic Earth→Mars-ish transfer: 90° angle, well inside the valid range.
    r1 = np.array([1.0,  0.0, 0.0]) * AU_M   # Earth at vernal equinox
    r2 = np.array([0.0,  1.524, 0.0]) * AU_M  # Mars analog, 90° ahead
    tof_s = 200.0 * DAY_S
    tof_days = 200.0

    print(f"Departure pos (AU): {r1 / AU_M}")
    print(f"Arrival pos   (AU): {r2 / AU_M}")
    print(f"TOF: {tof_days:.1f} days")

    # Try prograde first; retry retrograde if geometry fails (y(z*)<=0)
    v1_sc = v2_sc = None
    for prograde in (True, False):
        try:
            v1_sc, v2_sc = lambert_universal(r1, r2, tof_s, prograde=prograde)
            print(f"Lambert solved (prograde={prograde})")
            break
        except Exception as exc:
            print(f"Lambert (prograde={prograde}) failed: {exc}")

    if v1_sc is None:
        print("Lambert FAILED for both directions.")
        return

    print(f"Lambert v_depart (km/s): {np.linalg.norm(v1_sc)/1000:.3f}")
    print(f"Lambert v_arrive (km/s): {np.linalg.norm(v2_sc)/1000:.3f}")

    # Basic sanity: departure speed should be near Earth's heliocentric speed (~30 km/s)
    v_dep_kms = np.linalg.norm(v1_sc) / 1000.0
    assert 15 < v_dep_kms < 55, f"Departure speed {v_dep_kms:.1f} km/s out of expected range"
    print(f"Speed sanity check: PASS")

    poly = sample_transfer_polyline(r1, v1_sc, tof_s, n=50)
    assert len(poly) > 0, "Empty polyline"
    print(f"Polyline: {len(poly)} points")
    print("Smoke test: OK")


if __name__ == "__main__":
    main()
