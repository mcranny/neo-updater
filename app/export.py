"""Export database records into the viewer's portable JSON format."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from app.config import get_settings
from app.database import fetch_all, initialize_database


def build_viewer_payload(path: Path | None = None) -> dict[str, Any]:
    initialize_database(path)
    rows = fetch_all(
        """
        SELECT
            a.designation, a.fullname,
            ca.close_approach_jd_tdb, ca.close_approach_text, ca.distance_au,
            oe.epoch_jd_tdb, oe.semi_major_axis_au, oe.eccentricity,
            oe.inclination_deg, oe.raan_deg, oe.arg_periapsis_deg,
            oe.mean_anomaly_deg, oe.perihelion_time_jd_tdb,
            ip.departure_jd_tdb, ip.arrival_jd_tdb, ip.tof_days,
            ip.departure_dv_kms, ip.arrival_dv_kms, ip.total_dv_kms,
            ip.c3_km2_s2, ip.leo_departure_dv_kms, ip.polyline_json
        FROM asteroids a
        JOIN close_approaches ca ON ca.designation = a.designation
        JOIN orbital_elements oe ON oe.designation = a.designation
        JOIN intercept_plans ip ON ip.approach_id = ca.approach_id
        ORDER BY ca.close_approach_jd_tdb ASC, ip.total_dv_kms ASC
        """,
        path=path,
    )
    objects: list[dict[str, Any]] = []
    for row in rows:
        objects.append(
            {
                "des": row["designation"],
                "name": row["fullname"] or row["designation"],
                "next_ca": {
                    "jd_tdb": row["close_approach_jd_tdb"],
                    "cd_tdb": row["close_approach_text"],
                    "dist_au": row["distance_au"],
                },
                "elements": {
                    "epoch_jd_tdb": row["epoch_jd_tdb"],
                    "a_AU": row["semi_major_axis_au"],
                    "e": row["eccentricity"],
                    "i_deg": row["inclination_deg"],
                    "raan_deg": row["raan_deg"],
                    "argp_deg": row["arg_periapsis_deg"],
                    "ma_deg": row["mean_anomaly_deg"],
                    "tp_jd_tdb": row["perihelion_time_jd_tdb"],
                },
                "intercept": {
                    "depart_jd_tdb": row["departure_jd_tdb"],
                    "arrive_jd_tdb": row["arrival_jd_tdb"],
                    "tof_days": row["tof_days"],
                    "dv_depart_kms": row["departure_dv_kms"],
                    "dv_arrive_kms": row["arrival_dv_kms"],
                    "dv_total_kms": row["total_dv_kms"],
                    "c3_km2_s2": row["c3_km2_s2"],
                    "leo_dv_kms": row["leo_departure_dv_kms"],
                    "lambert_polyline_xyz_au": json.loads(row["polyline_json"] or "[]"),
                },
            }
        )
    return {
        "schema": "asteroid_intercept_planner/2.0",
        "source": "SQLite export",
        "count": len(objects),
        "objects": objects,
    }


def export_viewer_json(output: Path | None = None, database: Path | None = None) -> Path:
    output_path = Path(output or get_settings().export_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    payload = build_viewer_payload(database)
    output_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return output_path
