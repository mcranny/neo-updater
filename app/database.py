"""SQLite persistence and read models for the asteroid application."""

from __future__ import annotations

import json
import sqlite3
from collections.abc import Iterator, Sequence
from contextlib import contextmanager
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

from app.config import get_settings


def utc_now() -> str:
    return datetime.now(UTC).isoformat(timespec="seconds")


def _dict_row(cursor: sqlite3.Cursor, row: tuple[Any, ...]) -> dict[str, Any]:
    return {column[0]: row[index] for index, column in enumerate(cursor.description)}


@contextmanager
def connect(path: Path | None = None) -> Iterator[sqlite3.Connection]:
    database_path = Path(path or get_settings().database_path)
    database_path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(database_path, timeout=30)
    connection.row_factory = _dict_row
    connection.execute("PRAGMA foreign_keys = ON")
    connection.execute("PRAGMA busy_timeout = 30000")
    try:
        yield connection
    finally:
        connection.close()


def initialize_database(path: Path | None = None) -> Path:
    database_path = Path(path or get_settings().database_path)
    schema = Path(__file__).with_name("schema.sql").read_text(encoding="utf-8")
    with connect(database_path) as connection:
        schema_version = int(connection.execute("PRAGMA user_version").fetchone()["user_version"])
        if schema_version < 2:
            connection.executescript(schema)
            connection.commit()
        elif schema_version < 3:
            existing = {
                row["name"]
                for row in connection.execute("PRAGMA table_info(intercept_plans)").fetchall()
            }
            if "capture_dv_kms" not in existing:
                connection.execute("ALTER TABLE intercept_plans ADD COLUMN capture_dv_kms REAL")
            if "stable_final_orbit" not in existing:
                connection.execute("ALTER TABLE intercept_plans ADD COLUMN stable_final_orbit INTEGER")
            if "capture_json" not in existing:
                connection.execute("ALTER TABLE intercept_plans ADD COLUMN capture_json TEXT")
            connection.executescript(schema)
            connection.commit()
    return database_path


def begin_ingestion(query: dict[str, Any], path: Path | None = None) -> int:
    initialize_database(path)
    with connect(path) as connection:
        cursor = connection.execute(
            """
            INSERT INTO ingestion_runs(started_at, status, source, query_json)
            VALUES (?, 'running', 'JPL SSD CAD + SBDB', ?)
            """,
            (utc_now(), json.dumps(query, sort_keys=True)),
        )
        connection.commit()
        return int(cursor.lastrowid)


def finish_ingestion(
    run_id: int,
    *,
    status: str,
    object_count: int = 0,
    error_message: str | None = None,
    path: Path | None = None,
) -> None:
    with connect(path) as connection:
        connection.execute(
            """
            UPDATE ingestion_runs
            SET finished_at = ?, status = ?, object_count = ?, error_message = ?
            WHERE run_id = ?
            """,
            (utc_now(), status, object_count, error_message, run_id),
        )
        connection.commit()


def upsert_payload(payload: dict[str, Any], path: Path | None = None) -> dict[str, int]:
    """Atomically persist an API payload and its computed transfer plans."""
    initialize_database(path)
    now = utc_now()
    counts = {"asteroids": 0, "approaches": 0, "elements": 0, "plans": 0}
    with connect(path) as connection:
        connection.execute("BEGIN IMMEDIATE")
        for obj in payload.get("objects", []):
            designation = str(obj["des"])
            ca = obj.get("next_ca") or {}
            elements = obj.get("elements") or {}
            connection.execute(
                """
                INSERT INTO asteroids(
                    designation, fullname, absolute_magnitude_h, diameter_km,
                    first_seen_at, updated_at
                ) VALUES (?, ?, ?, ?, ?, ?)
                ON CONFLICT(designation) DO UPDATE SET
                    fullname = COALESCE(excluded.fullname, asteroids.fullname),
                    absolute_magnitude_h = COALESCE(
                        excluded.absolute_magnitude_h, asteroids.absolute_magnitude_h
                    ),
                    diameter_km = COALESCE(excluded.diameter_km, asteroids.diameter_km),
                    updated_at = excluded.updated_at
                """,
                (
                    designation,
                    ca.get("fullname"),
                    ca.get("h") if ca.get("h") is not None else elements.get("H"),
                    ca.get("diameter_km"),
                    now,
                    now,
                ),
            )
            counts["asteroids"] += 1

            if elements and elements.get("a_AU") is not None:
                connection.execute(
                    """
                    INSERT INTO orbital_elements(
                        designation, epoch_jd_tdb, semi_major_axis_au, eccentricity,
                        inclination_deg, raan_deg, arg_periapsis_deg, mean_anomaly_deg,
                        perihelion_time_jd_tdb, perihelion_distance_au,
                        aphelion_distance_au, orbital_period_days, earth_moid_au,
                        condition_code, orbit_id, source, updated_at
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(designation) DO UPDATE SET
                        epoch_jd_tdb=excluded.epoch_jd_tdb,
                        semi_major_axis_au=excluded.semi_major_axis_au,
                        eccentricity=excluded.eccentricity,
                        inclination_deg=excluded.inclination_deg,
                        raan_deg=excluded.raan_deg,
                        arg_periapsis_deg=excluded.arg_periapsis_deg,
                        mean_anomaly_deg=excluded.mean_anomaly_deg,
                        perihelion_time_jd_tdb=excluded.perihelion_time_jd_tdb,
                        perihelion_distance_au=excluded.perihelion_distance_au,
                        aphelion_distance_au=excluded.aphelion_distance_au,
                        orbital_period_days=excluded.orbital_period_days,
                        earth_moid_au=excluded.earth_moid_au,
                        condition_code=excluded.condition_code,
                        orbit_id=excluded.orbit_id,
                        source=excluded.source,
                        updated_at=excluded.updated_at
                    """,
                    (
                        designation,
                        elements.get("epoch_jd_tdb"),
                        elements["a_AU"],
                        elements["e"],
                        elements["i_deg"],
                        elements["raan_deg"],
                        elements["argp_deg"],
                        elements.get("ma_deg"),
                        elements.get("tp_jd_tdb"),
                        elements.get("q_AU"),
                        elements.get("Q_AU"),
                        elements.get("period_d"),
                        elements.get("moid_au"),
                        elements.get("condition_code"),
                        elements.get("orbit_id"),
                        elements.get("source") or "JPL SBDB",
                        now,
                    ),
                )
                counts["elements"] += 1

            if ca.get("jd_tdb") is None or ca.get("dist_au") is None:
                continue
            connection.execute(
                """
                INSERT INTO close_approaches(
                    designation, body, close_approach_jd_tdb, close_approach_text,
                    distance_au, distance_min_au, distance_max_au,
                    relative_velocity_kms, hyperbolic_excess_velocity_kms,
                    time_uncertainty, ingested_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(designation, body, close_approach_jd_tdb) DO UPDATE SET
                    close_approach_text=excluded.close_approach_text,
                    distance_au=excluded.distance_au,
                    distance_min_au=excluded.distance_min_au,
                    distance_max_au=excluded.distance_max_au,
                    relative_velocity_kms=excluded.relative_velocity_kms,
                    hyperbolic_excess_velocity_kms=excluded.hyperbolic_excess_velocity_kms,
                    time_uncertainty=excluded.time_uncertainty,
                    ingested_at=excluded.ingested_at
                """,
                (
                    designation,
                    ca.get("body") or "Earth",
                    ca["jd_tdb"],
                    ca.get("cd_tdb"),
                    ca["dist_au"],
                    ca.get("dist_min_au"),
                    ca.get("dist_max_au"),
                    ca.get("v_rel_kms"),
                    ca.get("v_inf_kms"),
                    ca.get("t_sigma_fmt"),
                    now,
                ),
            )
            approach = connection.execute(
                """
                SELECT approach_id FROM close_approaches
                WHERE designation=? AND body=? AND close_approach_jd_tdb=?
                """,
                (designation, ca.get("body") or "Earth", ca["jd_tdb"]),
            ).fetchone()
            approach_id = int(approach["approach_id"])
            counts["approaches"] += 1

            plan = obj.get("intercept")
            if plan:
                departure = plan.get("depart_jd_tdb", plan.get("departure_jd_tdb"))
                arrival = plan.get("arrive_jd_tdb", plan.get("arrival_jd_tdb"))
                total_dv = float(plan["dv_sum_mps"]) / 1000.0
                capture = plan.get("capture") if isinstance(plan.get("capture"), dict) else None
                connection.execute(
                    """
                    INSERT INTO intercept_plans(
                        approach_id, method, departure_jd_tdb, arrival_jd_tdb,
                        tof_days, departure_dv_kms, arrival_dv_kms, total_dv_kms,
                        c3_km2_s2, leo_departure_dv_kms, capture_dv_kms,
                        stable_final_orbit, capture_json, polyline_json, computed_at
                    ) VALUES (?, 'lambert-universal-v1', ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(approach_id, method) DO UPDATE SET
                        departure_jd_tdb=excluded.departure_jd_tdb,
                        arrival_jd_tdb=excluded.arrival_jd_tdb,
                        tof_days=excluded.tof_days,
                        departure_dv_kms=excluded.departure_dv_kms,
                        arrival_dv_kms=excluded.arrival_dv_kms,
                        total_dv_kms=excluded.total_dv_kms,
                        c3_km2_s2=excluded.c3_km2_s2,
                        leo_departure_dv_kms=excluded.leo_departure_dv_kms,
                        capture_dv_kms=excluded.capture_dv_kms,
                        stable_final_orbit=excluded.stable_final_orbit,
                        capture_json=excluded.capture_json,
                        polyline_json=excluded.polyline_json,
                        computed_at=excluded.computed_at
                    """,
                    (
                        approach_id,
                        departure,
                        arrival,
                        plan["tof_days"],
                        plan["dv_depart_kms"],
                        plan["dv_arrive_kms"],
                        total_dv,
                        plan.get("c3_km2_s2"),
                        plan.get("leo_dv_kms"),
                        capture.get("capture_dv_kms") if capture else None,
                        int(bool(capture.get("stable_final_orbit"))) if capture else None,
                        json.dumps(capture) if capture else None,
                        json.dumps(
                            plan.get("lambert_polyline_xyz_au")
                            or plan.get("lambert_polyline_xy_au")
                            or []
                        ),
                        now,
                    ),
                )
                counts["plans"] += 1
        connection.commit()
    return counts


def fetch_all(
    query: str, params: Sequence[Any] = (), path: Path | None = None
) -> list[dict[str, Any]]:
    initialize_database(path)
    with connect(path) as connection:
        return connection.execute(query, params).fetchall()


def fetch_one(
    query: str, params: Sequence[Any] = (), path: Path | None = None
) -> dict[str, Any] | None:
    rows = fetch_all(query, params, path)
    return rows[0] if rows else None


def dashboard_stats(path: Path | None = None) -> dict[str, Any]:
    return (
        fetch_one(
            """
        SELECT
          (SELECT COUNT(*) FROM asteroids) AS asteroid_count,
          (SELECT COUNT(*) FROM close_approaches) AS approach_count,
          (SELECT COUNT(*) FROM intercept_plans) AS plan_count,
          (SELECT MIN(distance_au) FROM close_approaches) AS nearest_distance_au,
          (SELECT MIN(total_dv_kms) FROM intercept_plans) AS lowest_total_dv_kms
        """,
            path=path,
        )
        or {}
    )


def list_asteroids(
    search: str = "",
    sort: str = "close_approach_jd_tdb",
    order: str = "asc",
    limit: int = 500,
    path: Path | None = None,
) -> list[dict[str, Any]]:
    allowed_sort = {
        "designation": "designation",
        "date": "close_approach_jd_tdb",
        "distance": "distance_au",
        "velocity": "relative_velocity_kms",
        "dv": "total_dv_kms",
        "inclination": "inclination_deg",
    }
    sort_column = allowed_sort.get(sort, "close_approach_jd_tdb")
    sort_direction = "DESC" if order.lower() == "desc" else "ASC"
    like = f"%{search.strip()}%"
    return fetch_all(
        f"""
        SELECT * FROM latest_asteroid_summary
        WHERE (? = '' OR designation LIKE ? OR COALESCE(fullname, '') LIKE ?)
        ORDER BY {sort_column} {sort_direction}, designation ASC
        LIMIT ?
        """,
        (search.strip(), like, like, limit),
        path,
    )


def get_asteroid(designation: str, path: Path | None = None) -> dict[str, Any] | None:
    return fetch_one(
        """
        SELECT a.*, oe.*
        FROM asteroids a
        LEFT JOIN orbital_elements oe USING(designation)
        WHERE a.designation = ?
        """,
        (designation,),
        path,
    )


def asteroid_approaches(designation: str, path: Path | None = None) -> list[dict[str, Any]]:
    return fetch_all(
        """
        SELECT ca.*, ip.*
        FROM close_approaches ca
        LEFT JOIN intercept_plans ip USING(approach_id)
        WHERE ca.designation = ?
        ORDER BY ca.close_approach_jd_tdb DESC
        """,
        (designation,),
        path,
    )


def latest_ingestion(path: Path | None = None) -> dict[str, Any] | None:
    return fetch_one("SELECT * FROM ingestion_runs ORDER BY run_id DESC LIMIT 1", path=path)


def list_tables(path: Path | None = None) -> list[dict[str, Any]]:
    return fetch_all(
        """
        SELECT name, type FROM sqlite_master
        WHERE type IN ('table', 'view') AND name NOT LIKE 'sqlite_%'
        ORDER BY type, name
        """,
        path=path,
    )


def run_readonly_query(
    query: str, *, limit: int = 5000, path: Path | None = None
) -> tuple[list[str], list[dict[str, Any]]]:
    normalized = query.strip().rstrip(";")
    if not normalized.lower().startswith(("select ", "with ")):
        raise ValueError("Only SELECT or WITH queries are allowed.")
    if ";" in normalized:
        raise ValueError("Only one SQL statement may be executed at a time.")
    with connect(path) as connection:
        connection.execute("PRAGMA query_only = ON")
        cursor = connection.execute(f"SELECT * FROM ({normalized}) LIMIT ?", (limit,))
        rows = cursor.fetchall()
        columns = [column[0] for column in cursor.description]
    return columns, rows
