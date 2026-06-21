from __future__ import annotations

import json

import numpy as np

from app.solar_system import (
    JD_J2000,
    PLANETS,
    asteroid_orbit_curve_au,
    format_jd,
    load_missions,
    orbit_curve_au,
    planet_elements,
    planet_position_au,
)


def test_planet_positions_and_orbits_are_finite_and_three_dimensional():
    expected_radius_ranges = {
        "Mercury": (0.25, 0.50),
        "Venus": (0.65, 0.80),
        "Earth": (0.95, 1.05),
        "Mars": (1.30, 1.70),
        "Jupiter": (4.8, 5.6),
        "Saturn": (8.8, 10.3),
        "Uranus": (18.0, 20.5),
        "Neptune": (29.0, 31.2),
    }
    for planet in PLANETS:
        position = planet_position_au(planet, 2461200.5)
        radius = float(np.linalg.norm(position))
        lower, upper = expected_radius_ranges[planet.name]
        assert position.shape == (3,)
        assert np.all(np.isfinite(position))
        assert lower < radius < upper
        curve = orbit_curve_au(planet_elements(planet, 2461200.5), samples=80)
        assert curve.shape == (80, 3)
        assert np.all(np.isfinite(curve))


def test_mission_loader_pads_legacy_2d_paths(tmp_path, sample_payload):
    path = tmp_path / "viewer.json"
    path.write_text(json.dumps(sample_payload), encoding="utf-8")
    missions = load_missions(path)
    assert len(missions) == 1
    assert missions[0].designation == "TEST-1"
    assert missions[0].polyline_au.shape == (2, 3)
    assert np.allclose(missions[0].polyline_au[:, 2], 0.0)
    curve = asteroid_orbit_curve_au(missions[0].elements, samples=80)
    assert curve.shape == (80, 3)


def test_julian_date_display():
    assert format_jd(JD_J2000) == "2000-01-01  12:00 UTC"
