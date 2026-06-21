PRAGMA foreign_keys = ON;

CREATE TABLE IF NOT EXISTS ingestion_runs (
    run_id INTEGER PRIMARY KEY AUTOINCREMENT,
    started_at TEXT NOT NULL,
    finished_at TEXT,
    status TEXT NOT NULL CHECK (status IN ('running', 'success', 'failed')),
    source TEXT NOT NULL,
    query_json TEXT NOT NULL,
    object_count INTEGER NOT NULL DEFAULT 0,
    error_message TEXT
);

CREATE TABLE IF NOT EXISTS asteroids (
    designation TEXT PRIMARY KEY,
    fullname TEXT,
    absolute_magnitude_h REAL,
    diameter_km REAL,
    first_seen_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS orbital_elements (
    designation TEXT PRIMARY KEY REFERENCES asteroids(designation) ON DELETE CASCADE,
    epoch_jd_tdb REAL,
    semi_major_axis_au REAL NOT NULL,
    eccentricity REAL NOT NULL,
    inclination_deg REAL NOT NULL,
    raan_deg REAL NOT NULL,
    arg_periapsis_deg REAL NOT NULL,
    mean_anomaly_deg REAL,
    perihelion_time_jd_tdb REAL,
    perihelion_distance_au REAL,
    aphelion_distance_au REAL,
    orbital_period_days REAL,
    earth_moid_au REAL,
    condition_code TEXT,
    orbit_id TEXT,
    source TEXT NOT NULL,
    updated_at TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS close_approaches (
    approach_id INTEGER PRIMARY KEY AUTOINCREMENT,
    designation TEXT NOT NULL REFERENCES asteroids(designation) ON DELETE CASCADE,
    body TEXT NOT NULL DEFAULT 'Earth',
    close_approach_jd_tdb REAL NOT NULL,
    close_approach_text TEXT,
    distance_au REAL NOT NULL,
    distance_min_au REAL,
    distance_max_au REAL,
    relative_velocity_kms REAL,
    hyperbolic_excess_velocity_kms REAL,
    time_uncertainty TEXT,
    ingested_at TEXT NOT NULL,
    UNIQUE (designation, body, close_approach_jd_tdb)
);

CREATE TABLE IF NOT EXISTS intercept_plans (
    plan_id INTEGER PRIMARY KEY AUTOINCREMENT,
    approach_id INTEGER NOT NULL REFERENCES close_approaches(approach_id) ON DELETE CASCADE,
    method TEXT NOT NULL,
    departure_jd_tdb REAL NOT NULL,
    arrival_jd_tdb REAL NOT NULL,
    tof_days REAL NOT NULL,
    departure_dv_kms REAL NOT NULL,
    arrival_dv_kms REAL NOT NULL,
    total_dv_kms REAL NOT NULL,
    c3_km2_s2 REAL,
    leo_departure_dv_kms REAL,
    polyline_json TEXT,
    computed_at TEXT NOT NULL,
    UNIQUE (approach_id, method)
);

CREATE INDEX IF NOT EXISTS idx_approaches_date
    ON close_approaches(close_approach_jd_tdb);
CREATE INDEX IF NOT EXISTS idx_approaches_designation
    ON close_approaches(designation);
CREATE INDEX IF NOT EXISTS idx_plans_total_dv
    ON intercept_plans(total_dv_kms);

DROP VIEW IF EXISTS latest_asteroid_summary;
CREATE VIEW latest_asteroid_summary AS
SELECT
    a.designation,
    a.fullname,
    a.absolute_magnitude_h,
    a.diameter_km,
    ca.approach_id,
    ca.close_approach_jd_tdb,
    ca.close_approach_text,
    ca.distance_au,
    ca.relative_velocity_kms,
    oe.semi_major_axis_au,
    oe.eccentricity,
    oe.inclination_deg,
    oe.earth_moid_au,
    ip.plan_id,
    ip.tof_days,
    ip.departure_dv_kms,
    ip.arrival_dv_kms,
    ip.total_dv_kms,
    ip.c3_km2_s2,
    ip.leo_departure_dv_kms
FROM asteroids a
JOIN close_approaches ca ON ca.approach_id = (
    SELECT ca2.approach_id
    FROM close_approaches ca2
    WHERE ca2.designation = a.designation
    ORDER BY
        CASE WHEN ca2.close_approach_jd_tdb >= julianday('now') THEN 0 ELSE 1 END,
        CASE WHEN ca2.close_approach_jd_tdb >= julianday('now')
             THEN ca2.close_approach_jd_tdb END ASC,
        CASE WHEN ca2.close_approach_jd_tdb < julianday('now')
             THEN ca2.close_approach_jd_tdb END DESC
    LIMIT 1
)
LEFT JOIN orbital_elements oe ON oe.designation = a.designation
LEFT JOIN intercept_plans ip ON ip.approach_id = ca.approach_id;
