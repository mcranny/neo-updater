from __future__ import annotations

import math

import numpy as np
import pytest

from scripts.orbital import (
    AU_M,
    MU_SUN,
    kepler_E_from_M,
    kepler_universal_propagate,
    lambert_universal,
    sample_transfer_polyline,
)


def test_kepler_equation_residual():
    mean_anomaly = 2.1
    eccentricity = 0.63
    eccentric_anomaly = kepler_E_from_M(mean_anomaly, eccentricity)
    assert eccentric_anomaly - eccentricity * math.sin(eccentric_anomaly) == pytest.approx(
        mean_anomaly, abs=1e-11
    )


def test_one_year_circular_orbit_returns_to_start():
    position = np.array([AU_M, 0.0, 0.0])
    velocity = np.array([0.0, math.sqrt(MU_SUN / AU_M), 0.0])
    period = 2.0 * math.pi * math.sqrt(AU_M**3 / MU_SUN)
    propagated, _ = kepler_universal_propagate(position, velocity, period)
    assert np.linalg.norm(propagated - position) < 1000.0


def test_lambert_solution_hits_requested_endpoint():
    departure = np.array([1.0, 0.0, 0.0]) * AU_M
    arrival = np.array([0.0, 1.524, 0.0]) * AU_M
    tof_seconds = 200.0 * 86400.0
    departure_velocity, _ = lambert_universal(departure, arrival, tof_seconds)
    propagated, _ = kepler_universal_propagate(departure, departure_velocity, tof_seconds)
    assert np.linalg.norm(propagated - arrival) < 25_000.0
    polyline = np.asarray(
        sample_transfer_polyline(departure, departure_velocity, tof_seconds, n=20)
    )
    assert polyline.shape == (21, 3)
    assert np.linalg.norm(polyline[-1] * AU_M - arrival) < 25_000.0


def test_lambert_rejects_nonpositive_time():
    with pytest.raises(ValueError, match="positive"):
        lambert_universal(np.ones(3), np.array([1.0, 2.0, 3.0]), 0.0)
