"""Environment-backed application configuration."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

try:
    from dotenv import load_dotenv
except ImportError:  # pragma: no cover - optional during minimal tooling
    load_dotenv = None


ROOT_DIR = Path(__file__).resolve().parents[1]
if load_dotenv is not None:
    load_dotenv(ROOT_DIR / ".env")


def _env_bool(name: str, default: bool = False) -> bool:
    value = os.getenv(name)
    if value is None:
        return default
    return value.strip().lower() in {"1", "true", "yes", "on"}


def _env_int(name: str, default: int) -> int:
    return int(os.getenv(name, str(default)))


def _env_float(name: str, default: float) -> float:
    return float(os.getenv(name, str(default)))


def _project_path(value: str | Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else ROOT_DIR / path


@dataclass(frozen=True)
class Settings:
    database_path: Path
    export_path: Path
    secret_key: str
    flask_host: str
    flask_port: int
    flask_debug: bool
    date_min: str
    date_max: str
    distance_max_au: str
    result_limit: int
    tof_min_days: int
    tof_max_days: int
    tof_step_days: int
    arrival_offsets_hours: tuple[int, ...]
    leo_altitude_km: float


def get_settings() -> Settings:
    offsets = tuple(
        int(part.strip())
        for part in os.getenv("INTERCEPT_ARRIVAL_OFFSETS_HOURS", "-12,-6,0,6,12").split(",")
        if part.strip()
    )
    return Settings(
        database_path=_project_path(os.getenv("ASTEROID_DATABASE_PATH", "data/asteroids.db")),
        export_path=_project_path(
            os.getenv(
                "ASTEROID_VIEWER_EXPORT_PATH",
                "data/latest_intercepts.json",
            )
        ),
        secret_key=os.getenv("SECRET_KEY", "local-development-only"),
        flask_host=os.getenv("FLASK_HOST", "127.0.0.1"),
        flask_port=_env_int("FLASK_PORT", 5000),
        flask_debug=_env_bool("FLASK_DEBUG", False),
        date_min=os.getenv("JPL_DATE_MIN", "now"),
        date_max=os.getenv("JPL_DATE_MAX", "+60"),
        distance_max_au=os.getenv("JPL_DISTANCE_MAX_AU", "0.05"),
        result_limit=_env_int("JPL_RESULT_LIMIT", 2000),
        tof_min_days=_env_int("INTERCEPT_TOF_MIN_DAYS", 30),
        tof_max_days=_env_int("INTERCEPT_TOF_MAX_DAYS", 180),
        tof_step_days=_env_int("INTERCEPT_TOF_STEP_DAYS", 10),
        arrival_offsets_hours=offsets,
        leo_altitude_km=_env_float("LEO_ALTITUDE_KM", 500.0),
    )
