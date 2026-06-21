from __future__ import annotations

import sqlite3

import pytest

from app.database import (
    dashboard_stats,
    initialize_database,
    list_asteroids,
    run_readonly_query,
    upsert_payload,
)
from app.export import build_viewer_payload


def test_schema_import_and_export_are_idempotent(tmp_path, sample_payload):
    database = tmp_path / "asteroids.db"
    initialize_database(database)
    first = upsert_payload(sample_payload, database)
    second = upsert_payload(sample_payload, database)

    assert first == {"asteroids": 1, "approaches": 1, "elements": 1, "plans": 1}
    assert second == first
    assert dashboard_stats(database)["asteroid_count"] == 1
    assert dashboard_stats(database)["approach_count"] == 1
    assert dashboard_stats(database)["plan_count"] == 1

    row = list_asteroids(path=database)[0]
    assert row["designation"] == "TEST-1"
    assert row["total_dv_kms"] == pytest.approx(5.3)

    exported = build_viewer_payload(database)
    assert exported["count"] == 1
    assert exported["objects"][0]["intercept"]["tof_days"] == 120
    assert exported["objects"][0]["intercept"]["lambert_polyline_xyz_au"] == [
        [1.0, 0.0],
        [0.7, 0.9],
    ]


def test_readonly_query_rejects_mutations(tmp_path, sample_payload):
    database = tmp_path / "asteroids.db"
    upsert_payload(sample_payload, database)
    columns, rows = run_readonly_query("SELECT designation FROM asteroids", path=database)
    assert columns == ["designation"]
    assert rows == [{"designation": "TEST-1"}]

    with pytest.raises(ValueError, match="Only SELECT"):
        run_readonly_query("DELETE FROM asteroids", path=database)
    with sqlite3.connect(database) as connection:
        assert connection.execute("SELECT COUNT(*) FROM asteroids").fetchone()[0] == 1
