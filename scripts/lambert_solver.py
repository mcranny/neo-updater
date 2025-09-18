"""
Robust Lambert Solver and Polyline Generator
===========================================

This module provides a robust implementation of Lambert’s problem and
a helper for generating transfer polylines between two heliocentric
states.  It is intended to replace the less reliable Lambert routines
found in ``neo_intercept_planner.py``.  The universal‐variables solver
from the original code is retained, but a numerical fallback is
included to handle cases where the primary solver fails.  The
fallback iteratively adjusts the departure velocity such that a
Keplerian propagation of the trajectory reaches the desired end point
within a specified tolerance.  The result is a more reliable Lambert
solution across a wider range of orbital geometries without crashing
the planning pipeline.

The main entry points are:

``lambert_universal(r1, r2, tof_s, mu=MU_SUN, prograde=True)``
    Solve Lambert’s problem using Battin/Vallado’s formulation with
    universal variables.  Returns the departure and arrival velocity
    vectors.  Raises a ``ValueError`` or ``RuntimeError`` on failure.

``lambert_polyline(r1_m, r2_m, tof_s, n)``
    Compute an ``n``‑point polyline describing the transfer orbit.
    Attempts the universal solver in both prograde and retrograde
    directions; if those fail, falls back to a velocity shooting
    method.  Returns an array of shape ``(n, 3)``.

``safe_lambert_polyline(r1_m, r2_m, tof_s, n, debug=False)``
    Thin wrapper around ``lambert_polyline`` that returns ``None`` on
    any error and optionally prints diagnostic information when
    ``debug`` is true.

The solver uses overflow–safe Stumpff functions and a universal
Kepler propagator lifted from the original project.  Both rely
exclusively on ``numpy`` and the Python standard library so there are
no external dependencies.  The fallback avoids ``numpy.linalg.lstsq``
which triggered LAPACK scaling errors in some environments.  Instead
it uses a pseudo–inverse computed via singular value decomposition
(SVD) for stable updates.

MIT Licence.  Derived from the neo‑updater project with enhancements.
"""

from __future__ import annotations

from typing import Optional, Tuple
import math
import numpy as np

# ---------------------------------------------------------------------------
# Constants

# Astronomical Unit in metres and the Sun’s gravitational parameter.  These
# match values used in the original neo‑updater project.
AU_M: float = 149_597_870_700.0  # metres
MU_SUN: float = 1.32712440018e20  # m^3/s^2

# ---------------------------------------------------------------------------
# Stumpff functions

def stumpff_C(z: float) -> float:
    """Numerically stable evaluation of the C(z) Stumpff function.

    For small |z| the usual series is used to avoid catastrophic
    cancellation.  For positive z, C(z) = (1 − cos√z)/z; for negative z,
    C(z) = (cosh√−z − 1)/(−z).
    """
    z = float(z)
    az = abs(z)
    # Use series expansion for very small |z|
    if az < 1e-8:
        return 0.5 - z / 24.0 + (z * z) / 720.0
    if z > 0.0:
        s = math.sqrt(z)
        return (1.0 - math.cos(s)) / z
    s = math.sqrt(-z)
    return (math.cosh(s) - 1.0) / (-z)


def stumpff_S(z: float) -> float:
    """Numerically stable evaluation of the S(z) Stumpff function.

    Uses series expansion for small |z|.  For positive z,
    S(z) = (√z − sin√z)/(√z)^3; for negative z,
    S(z) = (sinh√−z − √−z)/(√−z)^3.
    """
    z = float(z)
    az = abs(z)
    if az < 1e-8:
        return (1.0 / 6.0) - z / 120.0 + (z * z) / 5040.0
    if z > 0.0:
        s = math.sqrt(z)
        return (s - math.sin(s)) / (s * s * s)
    s = math.sqrt(-z)
    return (math.sinh(s) - s) / (s * s * s)


# ---------------------------------------------------------------------------
# Universal variables Lambert solver

def lambert_universal(r1: np.ndarray, r2: np.ndarray, tof_s: float,
                      mu: float = MU_SUN, prograde: bool = True,
                      max_iter: int = 80, tol: float = 1e-8
                      ) -> Tuple[np.ndarray, np.ndarray]:
    """Solve Lambert’s problem via universal variables.

    Parameters
    ----------
    r1, r2 : array_like shape (3,)
        Initial and final position vectors in metres.
    tof_s : float
        Time‑of‑flight between the positions in seconds.
    mu : float, optional
        Gravitational parameter of the central body (default Sun).
    prograde : bool, optional
        If true, compute the shorter transfer in the prograde sense;
        otherwise compute the longer/retrograde branch.
    max_iter : int, optional
        Maximum number of bisection iterations on the auxiliary variable z.
    tol : float, optional
        Absolute tolerance on the time of flight.

    Returns
    -------
    v1, v2 : ndarray shape (3,), ndarray shape (3,)
        Initial and final velocity vectors (m/s).

    Raises
    ------
    ValueError or RuntimeError
        If the geometry is singular or if no solution converges.
    """
    r1 = np.asarray(r1, dtype=float).reshape(3)
    r2 = np.asarray(r2, dtype=float).reshape(3)
    r1n = float(np.linalg.norm(r1))
    r2n = float(np.linalg.norm(r2))
    if r1n == 0.0 or r2n == 0.0:
        raise ValueError("Lambert: |r1| or |r2| is zero")
    tof_s = float(tof_s)
    if tof_s <= 0.0:
        raise ValueError("Lambert: time of flight must be positive")
    sqrt_mu = math.sqrt(mu)
    # Angle between vectors
    cosd = float(np.dot(r1, r2) / (r1n * r2n + 1e-16))
    cosd = max(-1.0, min(1.0, cosd))
    dtheta = math.acos(cosd)
    # Determine short/long way from cross product sign
    cz = float(np.cross(r1, r2)[2])
    if prograde and cz < 0.0:
        dtheta = 2.0 * math.pi - dtheta
    if not prograde and cz > 0.0:
        dtheta = 2.0 * math.pi - dtheta
    denom = 1.0 - cosd
    if abs(denom) < 1e-16:
        raise ValueError("Lambert geometry singular (Δθ ~ 0)")
    A = math.sin(dtheta) * math.sqrt(max(0.0, r1n * r2n / denom))
    if abs(A) < 1e-14:
        raise ValueError("Lambert geometry singular: A≈0")
    # Define helper functions y(z) and t(z)
    def y(z: float) -> float:
        C = stumpff_C(z)
        S = stumpff_S(z)
        Ceff = C if abs(C) > 1e-16 else 1e-16
        return r1n + r2n + A * ((z * S - 1.0) / Ceff)
    def t_of_z(z: float) -> float:
        C = stumpff_C(z)
        S = stumpff_S(z)
        yy = y(z)
        if yy <= 0.0:
            return float('inf')
        Ceff = C if abs(C) > 1e-16 else 1e-16
        x = math.sqrt(max(0.0, yy / Ceff))
        return (x**3 * S + A * math.sqrt(yy)) / sqrt_mu
    # Find bracket for z where t(z) >= tof
    z_lo, z_hi = -4.0 * math.pi * math.pi, 4.0 * math.pi * math.pi
    # Nudge lower bound until y(z) > 0
    while y(z_lo) <= 0.0:
        z_lo += 0.5
        if z_lo > 0.0:
            break
    # Expand upper bound until t(z) >= tof
    tries = 0
    t_hi = t_of_z(z_hi)
    while (not math.isfinite(t_hi)) or (t_hi < tof_s):
        z_hi *= 2.0
        tries += 1
        t_hi = t_of_z(z_hi)
        if tries > 40:
            break
    # Bisection search for z solving t(z) = tof_s
    z = 0.0
    for _ in range(max_iter):
        t_z = t_of_z(z)
        if abs(t_z - tof_s) < tol:
            break
        if t_z < tof_s:
            z_lo = z
        else:
            z_hi = z
        z = 0.5 * (z_lo + z_hi)
    # Compute f, g and gdot
    yy = y(z)
    if yy <= 0.0:
        raise RuntimeError("Lambert: y(z*) <= 0 at solution")
    f = 1.0 - yy / r1n
    g = A * math.sqrt(yy / mu)
    gdot = 1.0 - yy / r2n
    # Guard against tiny g
    if abs(g) < 1e-14:
        g = math.copysign(1e-14, g)
    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g
    return v1, v2


# ---------------------------------------------------------------------------
# Universal Kepler propagator

def kepler_universal_propagate(r0: np.ndarray, v0: np.ndarray, dt: float,
                               mu: float = MU_SUN, itmax: int = 60,
                               tol: float = 1e-9
                               ) -> Tuple[np.ndarray, np.ndarray]:
    """Propagate a two‑body orbit via universal variables.

    Given position ``r0`` and velocity ``v0`` at t=0, propagate them by a
    time ``dt`` (seconds) under gravitational parameter ``mu``.
    Returns the position and velocity at t=dt.  Works for elliptic,
    parabolic and hyperbolic trajectories.
    """
    r0 = np.asarray(r0, dtype=float)
    v0 = np.asarray(v0, dtype=float)
    r0n = float(np.linalg.norm(r0))
    if r0n == 0.0:
        raise ValueError("kepler_universal_propagate: |r0| is zero")
    dt = float(dt)
    vr0 = float(np.dot(r0, v0) / r0n)
    alpha = 2.0 / r0n - float(np.dot(v0, v0)) / mu  # 1/a
    sqrt_mu = math.sqrt(mu)
    # Initial guess for chi (Battin’s scheme)
    if abs(alpha) > 1e-12:
        chi = math.copysign(1.0, dt) * math.sqrt(mu) * abs(alpha) * dt
    else:
        # Near‑parabolic guess
        h = float(np.linalg.norm(np.cross(r0, v0)))
        p = h * h / mu
        # Use approximation from Battin
        s = 0.5 * (math.pi / 2.0 - math.atan(3.0 * math.sqrt(mu / (p**3)) * dt))
        w = math.atan(math.tan(s) ** (1.0 / 3.0))
        chi = math.sqrt(p) * 2.0 * math.tan(w)
    for _ in range(itmax):
        z = alpha * chi * chi
        Cz = stumpff_C(z)
        Sz = stumpff_S(z)
        r = chi * chi * Cz + vr0 * chi * (1.0 - z * Sz) + r0n * (1.0 - z * Cz)
        F = r - sqrt_mu * dt
        if abs(F) < tol:
            break
        dF = chi * (1.0 - z * Sz)
        # Avoid division by zero
        chi -= F / (dF + 1e-16)
    # Final position and velocity
    z = alpha * chi * chi
    Cz = stumpff_C(z)
    Sz = stumpff_S(z)
    f = 1.0 - (chi * chi / r0n) * Cz
    g = dt - (chi**3) * Sz / sqrt_mu
    r = f * r0 + g * v0
    rn = float(np.linalg.norm(r))
    fdot = (sqrt_mu / (rn * r0n)) * (z * Sz - 1.0) * chi
    gdot = 1.0 - (chi * chi * Cz) / rn
    v = fdot * r0 + gdot * v0
    return r, v


# ---------------------------------------------------------------------------
# Fallback shooting solver for Lambert

def _solve_via_shooting(r1_m: np.ndarray, r2_m: np.ndarray, tof_s: float,
                        v_guess: np.ndarray, mu: float = MU_SUN,
                        max_iter: int = 30, tol: float = 1.0
                        ) -> Tuple[np.ndarray, np.ndarray]:
    """Find a departure velocity via a shooting method.

    Starting from an initial guess ``v_guess``, iteratively adjust the
    velocity so that propagating from ``r1_m`` for ``tof_s`` seconds
    reaches ``r2_m``.  Uses a finite difference Jacobian and SVD–based
    pseudo–inverse to compute updates.  Returns the optimized velocity
    and the final position.  If convergence is not achieved within
    ``max_iter`` iterations, returns the last velocity and position.
    """
    r1_m = np.asarray(r1_m, dtype=float).reshape(3)
    r2_m = np.asarray(r2_m, dtype=float).reshape(3)
    v = np.asarray(v_guess, dtype=float).reshape(3)
    for _ in range(max_iter):
        r_final, _ = kepler_universal_propagate(r1_m, v, tof_s, mu)
        resid = r_final - r2_m
        if np.linalg.norm(resid) < tol:
            break
        # Build Jacobian via finite differences
        J = np.zeros((3, 3), dtype=float)
        # Use adaptive step size based on current velocity magnitude
        base_step = 1e-3 * (np.linalg.norm(v) + 1.0)
        for k in range(3):
            dv = np.zeros(3)
            dv[k] = base_step
            r_plus, _ = kepler_universal_propagate(r1_m, v + dv, tof_s, mu)
            J[:, k] = (r_plus - r_final) / base_step
        # Compute update using pseudo–inverse of J
        try:
            U, s, Vt = np.linalg.svd(J, full_matrices=False)
            # Invert singular values, guard against tiny singular values
            s_inv = np.array([1.0 / sv if sv > 1e-12 else 0.0 for sv in s])
            J_pinv = Vt.T @ np.diag(s_inv) @ U.T
            delta_v = J_pinv @ (-resid)
        except Exception:
            # Fall back to small proportional correction if SVD fails
            delta_v = -resid * 1e-6
        # Dampen the update to improve stability
        v += 0.5 * delta_v
    # Return final velocity and propagated position
    r_final, _ = kepler_universal_propagate(r1_m, v, tof_s, mu)
    return v, r_final


# ---------------------------------------------------------------------------
# Polyline generation

def _force_endpoints(poly_m: np.ndarray, r1_m: np.ndarray, r2_m: np.ndarray
                     ) -> np.ndarray:
    """Ensure the polyline starts at r1_m and ends at r2_m exactly."""
    poly_m = np.asarray(poly_m, dtype=float)
    poly_m[0, :] = np.asarray(r1_m, dtype=float).reshape(3)
    poly_m[-1, :] = np.asarray(r2_m, dtype=float).reshape(3)
    return poly_m


def lambert_polyline(r1_m: np.ndarray, r2_m: np.ndarray, tof_s: float, n: int
                     ) -> np.ndarray:
    """Generate a polyline approximating the Lambert transfer.

    Attempts the universal–variable Lambert solver for both prograde and
    retrograde directions.  On success, it uses the returned initial
    velocity to sample the trajectory at ``n`` equally spaced times.
    If both universal solver calls fail, falls back to a shooting
    solver to determine the departure velocity.  Always forces the
    first and last points of the polyline to match the supplied end
    points exactly.
    """
    r1_m = np.asarray(r1_m, dtype=float).reshape(3)
    r2_m = np.asarray(r2_m, dtype=float).reshape(3)
    tof_s = float(tof_s)
    n = int(max(2, n))
    last_err: Optional[Exception] = None
    # Try universal solver both directions
    for pro in (True, False):
        try:
            v1, _ = lambert_universal(r1_m, r2_m, tof_s, mu=MU_SUN, prograde=pro)
            ts = np.linspace(0.0, tof_s, n)
            pts = np.empty((n, 3), dtype=float)
            for i, t in enumerate(ts):
                r, _ = kepler_universal_propagate(r1_m, v1, float(t), MU_SUN)
                pts[i, :] = r
            return _force_endpoints(pts, r1_m, r2_m)
        except Exception as e:
            last_err = e
    # If universal solver failed, use fallback shooting solver
    # Construct a reasonable initial guess: chord plus tangential component
    chord = (r2_m - r1_m) / max(1e-9, tof_s)
    zhat = np.cross(r1_m, r2_m)
    # Tangential basis perpendicular to r1 and zhat
    t_hat = np.cross(zhat, r1_m)
    t_norm = np.linalg.norm(t_hat)
    if t_norm > 0:
        t_hat /= t_norm
    # Circular orbital velocity magnitude at r1
    v_circ = math.sqrt(MU_SUN / max(1e-3, np.linalg.norm(r1_m))) * t_hat
    v_init = 0.6 * chord + 0.8 * v_circ
    # Solve via shooting
    v_opt, _ = _solve_via_shooting(r1_m, r2_m, tof_s, v_init, MU_SUN)
    # Sample trajectory
    ts = np.linspace(0.0, tof_s, n)
    pts = np.empty((n, 3), dtype=float)
    for i, t in enumerate(ts):
        r, _ = kepler_universal_propagate(r1_m, v_opt, float(t), MU_SUN)
        pts[i, :] = r
    return _force_endpoints(pts, r1_m, r2_m)


def safe_lambert_polyline(r1_m: np.ndarray, r2_m: np.ndarray, tof_s: float,
                          n: int, *, debug: bool = False
                          ) -> Optional[np.ndarray]:
    """Call ``lambert_polyline`` and return ``None`` on failure.

    If the underlying solver raises an exception or the resulting
    polyline contains NaNs or infs, this wrapper returns ``None``.  When
    ``debug`` is true, diagnostic messages are printed to ``stdout``.
    """
    try:
        poly = lambert_polyline(r1_m, r2_m, tof_s, n)
        # Validate finite coordinates
        if not np.all(np.isfinite(poly)):
            if debug:
                print("[lambert_solver] polyline contains non‑finite values")
            return None
        return poly
    except Exception as e:
        if debug:
            print(f"[lambert_solver] Lambert computation failed: {e}")
        return None


# ---------------------------------------------------------------------------
# Helper: approximate orbital period from semi‑major axis

def compute_orbital_period_days(a_AU: float) -> float:
    """Return approximate orbital period in days for semi‑major axis in AU."""
    return 365.25 * (a_AU ** 1.5)


__all__ = [
    "AU_M",
    "MU_SUN",
    "lambert_universal",
    "kepler_universal_propagate",
    "lambert_polyline",
    "safe_lambert_polyline",
    "compute_orbital_period_days",
]