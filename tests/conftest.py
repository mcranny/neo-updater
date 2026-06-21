from __future__ import annotations

from typing import Any

import pytest


@pytest.fixture
def sample_payload() -> dict[str, Any]:
    return {
        "schema": "asteroid_intercept_planner/2.0",
        "objects": [
            {
                "des": "TEST-1",
                "next_ca": {
                    "fullname": "Synthetic test asteroid",
                    "jd_tdb": 2461200.5,
                    "cd_tdb": "2026-Jul-09 00:00",
                    "dist_au": 0.0123,
                    "dist_min_au": 0.0122,
                    "dist_max_au": 0.0124,
                    "v_rel_kms": 12.5,
                    "v_inf_kms": 11.9,
                    "h": 22.1,
                    "diameter_km": 0.12,
                    "body": "Earth",
                },
                "elements": {
                    "epoch_jd_tdb": 2461000.5,
                    "a_AU": 1.2,
                    "e": 0.1,
                    "i_deg": 5.0,
                    "raan_deg": 50.0,
                    "argp_deg": 10.0,
                    "ma_deg": 20.0,
                    "tp_jd_tdb": 2460900.5,
                    "q_AU": 1.08,
                    "Q_AU": 1.32,
                    "period_d": 480.0,
                    "moid_au": 0.01,
                    "condition_code": "1",
                    "orbit_id": "42",
                    "source": "JPL SBDB",
                },
                "intercept": {
                    "depart_jd_tdb": 2461080.5,
                    "arrive_jd_tdb": 2461200.5,
                    "tof_days": 120,
                    "dv_depart_kms": 3.2,
                    "dv_arrive_kms": 2.1,
                    "dv_sum_mps": 5300.0,
                    "c3_km2_s2": 10.24,
                    "leo_dv_kms": 3.8,
                    "lambert_polyline_xy_au": [[1.0, 0.0], [0.7, 0.9]],
                },
            }
        ],
    }
