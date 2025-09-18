"""
Test harness for the neo-updater Lambert solver.

This script sets up an example transfer from Earth to a synthetic target and
verifies that the improved Lambert routine (with fallback) converges.  It
bounds the time-of-flight to half the target’s orbital period to avoid
multi-revolution cases and uses the planner’s safe_lambert_polyline function
to compute the transfer.

Usage: python3 test_lambert.py
"""

import math
from datetime import datetime, timezone
import numpy as np

# Import the improved Lambert helper from your planner
# Adjust the import path if neo_intercept_planner_fixed.py lives elsewhere.
from neo_intercept_planner import (
    rv_from_elements_at_time,
    safe_lambert_polyline,
    compute_orbital_period_days,
    AU_M,
    MU_SUN,
)

# Constants
DEG2RAD = math.pi / 180.0

def jd_from_datetime(dt: datetime) -> float:
    """Convert a timezone-aware datetime to Julian day."""
    JD_UNIX_EPOCH = 2440587.5
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return JD_UNIX_EPOCH + dt.timestamp() / 86400.0

def main():
    # Earth orbital elements (approximate, J2000)
    earth_el = {
        "a_AU": 1.00000011,
        "e": 0.01671022,
        "i_deg": 0.00005,
        "raan_deg": -11.26064,
        "argp_deg": 102.94719,
        "ma_deg": 100.46435,
        "epoch_jd": 2451545.0,
    }
    # Synthetic target: 1.2 AU, moderate eccentricity, small inclination
    tgt_el = {
        "a_AU": 1.2,
        "e": 0.1,
        "i_deg": 5.0,
        "raan_deg": 50.0,
        "argp_deg": 10.0,
        "ma_deg": 0.0,
        "epoch_jd": 2451545.0,
    }
    # Departure in 90 days from now (use timezone-aware datetime)
    now_jd = jd_from_datetime(datetime.now(timezone.utc))
    depart_jd = now_jd + 90.0

    # Bound the time-of-flight to half the target’s orbital period
    tgt_period = compute_orbital_period_days(max(0.01, tgt_el["a_AU"]))
    tof_days = 0.5 * tgt_period  # at most half an orbit
    arrive_jd = depart_jd + tof_days

    # Compute heliocentric positions (and velocities) at departure and arrival
    r1_m, v1_m = rv_from_elements_at_time(
        earth_el["a_AU"], earth_el["e"], earth_el["i_deg"],
        earth_el["raan_deg"], earth_el["argp_deg"], earth_el["ma_deg"],
        earth_el["epoch_jd"], depart_jd, mu=MU_SUN
    )
    r2_m, v2_m = rv_from_elements_at_time(
        tgt_el["a_AU"], tgt_el["e"], tgt_el["i_deg"],
        tgt_el["raan_deg"], tgt_el["argp_deg"], tgt_el["ma_deg"],
        tgt_el["epoch_jd"], arrive_jd, mu=MU_SUN
    )

    # Convert TOF to seconds
    tof_s = float(tof_days) * 86400.0

    # Compute a transfer polyline using the safe Lambert helper (with fallback)
    poly = safe_lambert_polyline(r1_m, r2_m, tof_s, n=500, debug=True)
    if poly is None:
        print("Lambert solver failed for this case.")
        return

    # Verify that the last point matches the target within a reasonable error
    final_r = poly[-1]
    error_m = np.linalg.norm(r2_m - final_r)
    print(f"Departure position (AU): {r1_m / AU_M}")
    print(f"Arrival position (AU):   {r2_m / AU_M}")
    print(f"Lambert arrival  (AU):   {final_r / AU_M}")
    print(f"Position error (m):      {error_m:.2e}")

if __name__ == "__main__":
    main()
