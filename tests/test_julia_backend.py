from __future__ import annotations

import json
import os
import shutil

from app.database import upsert_payload
from app.export import build_viewer_payload
from app.solar_system import load_missions
from scripts.plan_intercepts import LOCAL_JULIA, attach_intercepts


def test_julia_backend_payload_persists_and_loads(tmp_path):
    if not os.getenv("JULIA") and not LOCAL_JULIA.exists() and shutil.which("julia") is None:
        raise AssertionError("Julia backend is required for calculation tests")

    payload = {
        "schema": "asteroid_intercept_planner/2.0",
        "objects": [
            {
                "des": "TEST-JULIA",
                "next_ca": {
                    "fullname": "Synthetic Julia asteroid",
                    "jd_tdb": 2461200.5,
                    "cd_tdb": "2026-Jul-09 00:00",
                    "dist_au": 0.0123,
                    "diameter_km": 1.0,
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
                },
            }
        ],
    }

    result = attach_intercepts(
        payload,
        tof_days_grid=[120],
        arrive_offset_hours=[0],
        leo_alt_m=500e3,
        final_orbits=3.0,
    )
    intercept = result["objects"][0]["intercept"]
    assert intercept is not None
    for key in (
        "depart_jd_tdb",
        "arrive_jd_tdb",
        "tof_days",
        "dv_depart_kms",
        "dv_arrive_kms",
        "dv_sum_mps",
        "c3_km2_s2",
        "vinf_kms",
        "leo_dv_kms",
        "lambert_polyline_xyz_au",
    ):
        assert key in intercept
    assert intercept["capture"]["model"] == "density-sphere-v1"
    assert intercept["capture"]["capture_dv_kms"] >= 0.0
    assert intercept["capture"]["requested_final_orbits"] == 3.0
    assert intercept["capture"]["orbit_polyline_heliocentric_au"]

    database = tmp_path / "asteroids.db"
    assert upsert_payload(result, database)["plans"] == 1
    exported = build_viewer_payload(database)
    exported_intercept = exported["objects"][0]["intercept"]
    assert exported_intercept["capture"]["model"] == "density-sphere-v1"

    viewer_payload = tmp_path / "viewer.json"
    viewer_payload.write_text(json.dumps(exported), encoding="utf-8")
    missions = load_missions(viewer_payload)
    assert missions[0].capture["model"] == "density-sphere-v1"
    assert missions[0].final_orbit_polyline_au is not None
    assert missions[0].capture_duration_days > 0.0
    assert missions[0].total_duration_days > missions[0].tof_days
