"""Fetch JPL close approaches, compute transfers, and update SQLite."""

from __future__ import annotations

import argparse
import json
import sys
from typing import Any

from app.config import Settings, get_settings
from app.database import begin_ingestion, finish_ingestion, upsert_payload
from app.export import export_viewer_json
from scripts.hud_metrics import compute_hud
from scripts.plan_intercepts import attach_intercepts
from scripts.sbdb_client import discover_earth_close_approaches


def collect_payload(settings: Settings) -> dict[str, Any]:
    records = discover_earth_close_approaches(
        date_min=settings.date_min,
        date_max=settings.date_max,
        dist_max=settings.distance_max_au,
        limit=settings.result_limit,
        attach_elements=True,
    )
    for obj in records:
        elements = obj.get("elements")
        approach = obj.get("next_ca")
        obj["hud"] = None
        if elements and approach and approach.get("jd_tdb") is not None:
            try:
                obj["hud"] = compute_hud(elements, float(approach["jd_tdb"]))
            except (KeyError, TypeError, ValueError, ArithmeticError):
                pass

    payload = {
        "schema": "asteroid_intercept_planner/2.0",
        "source": "JPL SSD CAD + SBDB",
        "query": {
            "date_min": settings.date_min,
            "date_max": settings.date_max,
            "dist_max": settings.distance_max_au,
            "limit": settings.result_limit,
        },
        "count": len(records),
        "objects": records,
    }
    tof_grid = list(
        range(
            settings.tof_min_days,
            settings.tof_max_days + 1,
            settings.tof_step_days,
        )
    )
    return attach_intercepts(
        payload,
        tof_days_grid=tof_grid,
        arrive_offset_hours=list(settings.arrival_offsets_hours),
        leo_alt_m=settings.leo_altitude_km * 1000.0,
    )


def run_ingestion(settings: Settings | None = None) -> dict[str, int]:
    settings = settings or get_settings()
    query = {
        "date_min": settings.date_min,
        "date_max": settings.date_max,
        "dist_max": settings.distance_max_au,
        "limit": settings.result_limit,
    }
    run_id = begin_ingestion(query, settings.database_path)
    try:
        payload = collect_payload(settings)
        counts = upsert_payload(payload, settings.database_path)
        finish_ingestion(
            run_id,
            status="success",
            object_count=counts["asteroids"],
            path=settings.database_path,
        )
        export_viewer_json(settings.export_path, settings.database_path)
        return counts
    except Exception as error:
        finish_ingestion(
            run_id,
            status="failed",
            error_message=str(error)[:2000],
            path=settings.database_path,
        )
        raise


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Update the asteroid SQLite database from JPL CAD and SBDB."
    )
    parser.add_argument("--print-json", action="store_true", help="Print import counts as JSON.")
    args = parser.parse_args()
    counts = run_ingestion()
    if args.print_json:
        print(json.dumps(counts, sort_keys=True))
    else:
        print("Updated database: " + ", ".join(f"{name}={count}" for name, count in counts.items()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
