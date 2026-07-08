"""Solar-system ephemeris helpers and viewer mission data models."""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from scripts.ephem_earth import earth_rv
from scripts.orbital import AU_M, deg2rad, kepler_E_from_M, oe_to_rv, wrap_2pi
from scripts.plan_intercepts import asteroid_rv_from_elements

JD_J2000 = 2451545.0
DAYS_PER_CENTURY = 36525.0


@dataclass(frozen=True)
class Planet:
    name: str
    color: tuple[float, float, float, float]
    marker_size: float
    elements: tuple[float, float, float, float, float, float]
    rates: tuple[float, float, float, float, float, float]


# JPL approximate heliocentric elements and rates per Julian century.
# Values are intended for 1800-2050 and use ecliptic J2000 angles in degrees.
PLANETS: tuple[Planet, ...] = (
    Planet(
        "Mercury",
        (0.72, 0.68, 0.62, 1.0),
        6.0,
        (0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593),
        (0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081),
    ),
    Planet(
        "Venus",
        (0.95, 0.72, 0.34, 1.0),
        7.0,
        (0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255),
        (0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418),
    ),
    Planet(
        "Earth",
        (0.26, 0.67, 1.0, 1.0),
        8.0,
        (1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0),
        (0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0),
    ),
    Planet(
        "Mars",
        (0.96, 0.35, 0.20, 1.0),
        7.0,
        (1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891),
        (0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343),
    ),
    Planet(
        "Jupiter",
        (0.88, 0.66, 0.43, 1.0),
        11.0,
        (5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909),
        (-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106),
    ),
    Planet(
        "Saturn",
        (0.91, 0.80, 0.52, 1.0),
        10.0,
        (9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448),
        (-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794),
    ),
    Planet(
        "Uranus",
        (0.46, 0.88, 0.91, 1.0),
        9.0,
        (19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503),
        (-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589),
    ),
    Planet(
        "Neptune",
        (0.28, 0.43, 0.96, 1.0),
        9.0,
        (30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574),
        (0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664),
    ),
)


@dataclass(frozen=True)
class Mission:
    designation: str
    name: str
    approach_text: str
    distance_au: float
    elements: dict[str, Any]
    departure_jd: float
    arrival_jd: float
    tof_days: float
    departure_dv_kms: float
    arrival_dv_kms: float
    total_dv_kms: float
    c3_km2_s2: float | None
    leo_dv_kms: float | None
    polyline_au: np.ndarray
    capture: dict[str, Any]
    final_orbit_polyline_au: np.ndarray | None

    @property
    def label(self) -> str:
        return f"{self.designation}  ·  {self.approach_text}  ·  Δv {self.total_dv_kms:.2f} km/s"

    @property
    def capture_duration_days(self) -> float:
        if self.final_orbit_polyline_au is None:
            return 0.0
        try:
            return max(0.0, float(self.capture.get("propagation_days") or 0.0))
        except (TypeError, ValueError):
            return 0.0

    @property
    def total_duration_days(self) -> float:
        return self.tof_days + self.capture_duration_days


def julian_centuries(jd_tdb: float) -> float:
    return (float(jd_tdb) - JD_J2000) / DAYS_PER_CENTURY


def planet_elements(planet: Planet, jd_tdb: float) -> dict[str, float]:
    century = julian_centuries(jd_tdb)
    values = tuple(
        base + rate * century for base, rate in zip(planet.elements, planet.rates, strict=True)
    )
    semi_major_axis, eccentricity, inclination, mean_longitude, longitude_peri, node = values
    return {
        "a_AU": semi_major_axis,
        "e": eccentricity,
        "i_deg": inclination,
        "raan_deg": node,
        "argp_deg": longitude_peri - node,
        "ma_deg": (mean_longitude - longitude_peri) % 360.0,
    }


def _state_from_elements(elements: dict[str, float]) -> tuple[np.ndarray, np.ndarray]:
    eccentric_anomaly = kepler_E_from_M(deg2rad(elements["ma_deg"]), elements["e"])
    true_anomaly = 2.0 * math.atan2(
        math.sqrt(1.0 + elements["e"]) * math.sin(eccentric_anomaly / 2.0),
        math.sqrt(1.0 - elements["e"]) * math.cos(eccentric_anomaly / 2.0),
    )
    return oe_to_rv(
        elements["a_AU"] * AU_M,
        elements["e"],
        deg2rad(elements["i_deg"]),
        deg2rad(elements["raan_deg"]),
        deg2rad(elements["argp_deg"]),
        wrap_2pi(true_anomaly),
    )


def planet_position_au(planet: Planet, jd_tdb: float) -> np.ndarray:
    if planet.name == "Earth":
        position, _ = earth_rv(jd_tdb)
    else:
        position, _ = _state_from_elements(planet_elements(planet, jd_tdb))
    return np.asarray(position, dtype=float) / AU_M


def orbit_curve_au(elements: dict[str, float], samples: int = 360) -> np.ndarray:
    anomalies = np.linspace(0.0, 2.0 * math.pi, max(64, samples), endpoint=True)
    points = [
        oe_to_rv(
            elements["a_AU"] * AU_M,
            elements["e"],
            deg2rad(elements["i_deg"]),
            deg2rad(elements["raan_deg"]),
            deg2rad(elements["argp_deg"]),
            anomaly,
        )[0]
        / AU_M
        for anomaly in anomalies
    ]
    return np.asarray(points, dtype=np.float32)


def asteroid_position_au(elements: dict[str, Any], jd_tdb: float) -> np.ndarray:
    position, _ = asteroid_rv_from_elements(elements, jd_tdb)
    return np.asarray(position, dtype=float) / AU_M


def asteroid_orbit_curve_au(elements: dict[str, Any], samples: int = 420) -> np.ndarray:
    required = ("a_AU", "e", "i_deg", "raan_deg", "argp_deg")
    normalized = {key: float(elements[key]) for key in required}
    return orbit_curve_au(normalized, samples)


def jd_to_datetime(jd_tdb: float) -> datetime:
    return datetime.fromtimestamp((float(jd_tdb) - 2440587.5) * 86400.0, UTC)


def format_jd(jd_tdb: float) -> str:
    return jd_to_datetime(jd_tdb).strftime("%Y-%m-%d  %H:%M UTC")


def _polyline_xyz(values: Any) -> np.ndarray:
    array = np.asarray(values or [], dtype=float)
    if array.ndim != 2 or array.shape[0] < 2 or array.shape[1] not in (2, 3):
        raise ValueError("mission has no valid transfer polyline")
    if array.shape[1] == 2:
        array = np.pad(array, ((0, 0), (0, 1)))
    if not np.all(np.isfinite(array)):
        raise ValueError("transfer polyline contains non-finite values")
    return array.astype(np.float32)


def load_missions(path: Path) -> list[Mission]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    missions: list[Mission] = []
    for obj in payload.get("objects") or []:
        intercept = obj.get("intercept") or {}
        elements = obj.get("elements") or {}
        departure_value = intercept.get("depart_jd_tdb")
        if departure_value is None:
            departure_value = intercept.get("departure_jd_tdb")
        arrival_value = intercept.get("arrive_jd_tdb")
        if arrival_value is None:
            arrival_value = intercept.get("arrival_jd_tdb")
        polyline = intercept.get("lambert_polyline_xyz_au")
        if not polyline:
            polyline = intercept.get("lambert_polyline_xy_au")
        capture = intercept.get("capture") if isinstance(intercept.get("capture"), dict) else {}
        final_orbit_polyline = capture.get("orbit_polyline_heliocentric_au") if capture else None
        try:
            departure = float(departure_value)
            arrival = float(arrival_value)
            mission = Mission(
                designation=str(obj.get("des") or obj.get("name") or "Unknown"),
                name=str(obj.get("name") or obj.get("des") or "Unknown"),
                approach_text=str((obj.get("next_ca") or {}).get("cd_tdb") or "Unknown approach"),
                distance_au=float((obj.get("next_ca") or {}).get("dist_au") or 0.0),
                elements=elements,
                departure_jd=departure,
                arrival_jd=arrival,
                tof_days=float(intercept["tof_days"]),
                departure_dv_kms=float(intercept["dv_depart_kms"]),
                arrival_dv_kms=float(intercept["dv_arrive_kms"]),
                total_dv_kms=float(
                    intercept.get("dv_total_kms") or intercept.get("dv_sum_mps", 0.0) / 1000.0
                ),
                c3_km2_s2=(
                    float(intercept["c3_km2_s2"])
                    if intercept.get("c3_km2_s2") is not None
                    else None
                ),
                leo_dv_kms=(
                    float(intercept["leo_dv_kms"])
                    if intercept.get("leo_dv_kms") is not None
                    else None
                ),
                polyline_au=_polyline_xyz(polyline),
                capture=capture,
                final_orbit_polyline_au=(
                    _polyline_xyz(final_orbit_polyline) if final_orbit_polyline else None
                ),
            )
        except (KeyError, TypeError, ValueError):
            continue
        missions.append(mission)
    if not missions:
        raise ValueError("No complete interception plans were found in the viewer export.")
    return missions
