# scripts/kepler_propagation.py
from __future__ import annotations
import math
import numpy as np

# If you already import this constant elsewhere, you can import it instead.
GAUSSIAN_GRAVITATIONAL_CONSTANT = 0.01720209895  # rad/day
MU_SUN = GAUSSIAN_GRAVITATIONAL_CONSTANT ** 2    # AU^3 / day^2

# ---- Stumpff functions ----
def _stumpff_C(z: float) -> float:
    if z > 0:
        sz = math.sqrt(z)
        return (1.0 - math.cos(sz)) / z
    if z < 0:
        sz = math.sqrt(-z)
        return (1.0 - math.cosh(sz)) / z
    return 0.5

def _stumpff_S(z: float) -> float:
    if z > 0:
        sz = math.sqrt(z)
        return (sz - math.sin(sz)) / (sz ** 3)
    if z < 0:
        sz = math.sqrt(-z)
        return (math.sinh(sz) - sz) / (sz ** 3)
    return 1.0 / 6.0

# ---- Universal Kepler propagator (r0,v0 in AU, AU/day; dt in days) ----
def propagate_universal(r0: np.ndarray, v0: np.ndarray, dt: float, mu: float = MU_SUN) -> tuple[np.ndarray, np.ndarray]:
    """
    Returns (r, v) at t0+dt from heliocentric 2-body motion using universal variables.
    r0 [AU], v0 [AU/day], dt [day], mu [AU^3/day^2]
    """
    r0 = np.asarray(r0, dtype=float).reshape(3)
    v0 = np.asarray(v0, dtype=float).reshape(3)
    r0n = float(np.linalg.norm(r0))
    if r0n == 0.0:
        raise ValueError("r0 is zero vector")

    v0n2 = float(np.dot(v0, v0))
    alpha = 2.0 / r0n - v0n2 / mu  # 1/a

    # Initial guess for chi
    sqrt_mu = math.sqrt(mu)
    if abs(alpha) > 1e-12:  # elliptical / hyperbolic
        chi = sqrt_mu * abs(alpha) * dt
    else:  # near parabolic
        # Battin suggestion
        h = np.cross(r0, v0)
        p = (np.dot(h, h)) / mu
        s = 0.5 * (math.pi / 2.0 - math.atan(3.0 * math.sqrt(mu / (p ** 3)) * dt))
        w = math.atan(math.tan(s) ** (1.0 / 3.0))
        chi = math.sqrt(p) * (2.0 / math.tan(2.0 * w))

    # Newton iterations
    r0dotv0 = float(np.dot(r0, v0))
    for _ in range(60):
        z = alpha * chi * chi
        C = _stumpff_C(z)
        S = _stumpff_S(z)
        r = chi * chi * C + r0dotv0 / sqrt_mu * chi * (1.0 - z * S) + r0n * (1.0 - z * C)
        F = r0n * r0dotv0 / sqrt_mu * chi * chi * C + (1.0 - alpha * r0n) * chi * chi * chi * S + r0n * chi - sqrt_mu * dt
        if abs(F) < 1e-12:
            break
        dF = r0n * r0dotv0 / sqrt_mu * chi * (1.0 - z * S) + (1.0 - alpha * r0n) * chi * chi * C + r0n
        chi -= F / dF

    # f, g
    z = alpha * chi * chi
    C = _stumpff_C(z)
    S = _stumpff_S(z)
    f = 1.0 - (chi * chi / r0n) * C
    g = dt - (chi ** 3 / sqrt_mu) * S

    r_vec = f * r0 + g * v0
    rn = float(np.linalg.norm(r_vec))
    # fdot, gdot
    fdot = (sqrt_mu / (rn * r0n)) * (z * S - 1.0) * chi
    gdot = 1.0 - (chi * chi / rn) * C
    v_vec = fdot * r0 + gdot * v0
    return r_vec, v_vec

def sample_transfer_polyline(r0: np.ndarray, v0: np.ndarray, tof_days: float, n: int, mu: float = MU_SUN) -> np.ndarray:
    """
    Uniformly sample positions along the conic defined by (r0,v0) for 'tof_days'.
    Returns (n,3) in AU.
    """
    n = max(2, int(n))
    ts = np.linspace(0.0, float(tof_days), n, dtype=float)
    pts = np.empty((n, 3), dtype=float)
    for i, t in enumerate(ts):
        r, _ = propagate_universal(r0, v0, t, mu)
        pts[i, :] = r
    return pts.astype(np.float32)
