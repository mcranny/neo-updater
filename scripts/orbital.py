# scripts/orbital.py
from __future__ import annotations
import math
import numpy as np
from typing import Tuple, Union

Number = Union[float, np.floating]

AU_M  = 149_597_870_700.0
MU_SUN = 1.32712440018e20  # m^3/s^2


# ---------------- Basic helpers ----------------
def deg2rad(x: Number) -> float:
    return float(x) * math.pi / 180.0

def wrap_2pi(x: Number) -> float:
    y = math.fmod(float(x), 2.0 * math.pi)
    return y + 2.0 * math.pi if y < 0.0 else y


def kepler_E_from_M(M: Number, e: float, tol: float = 1e-12, itmax: int = 60) -> float:
    M = wrap_2pi(float(M))
    if e < 1e-12:
        return M
    E = M if e < 0.8 else math.pi
    for _ in range(itmax):
        f = E - e * math.sin(E) - M
        fp = 1.0 - e * math.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            break
    return float(E)


def oe_to_rv(a_m: Number, e: float, i: Number, raan: Number, argp: Number, nu: Number,
             mu: float = MU_SUN) -> Tuple[np.ndarray, np.ndarray]:
    """Orbital elements -> (r, v) in IJK; a in meters; angles in radians."""
    a_m  = float(a_m); i = float(i); raan = float(raan); argp = float(argp); nu = float(nu)
    p = a_m * (1.0 - e * e)

    # Perifocal
    cnu, snu = math.cos(nu), math.sin(nu)
    r_pf = np.array([p * cnu / (1.0 + e * cnu),
                     p * snu / (1.0 + e * cnu),
                     0.0], dtype=float)
    v_pf = np.array([-math.sqrt(mu / p) * snu,
                      math.sqrt(mu / p) * (e + cnu),
                      0.0], dtype=float)

    ci, si = math.cos(i), math.sin(i)
    cO, sO = math.cos(raan), math.sin(raan)
    co, so = math.cos(argp), math.sin(argp)

    # PQW -> IJK rotation
    R = np.array([
        [ cO * co - sO * so * ci, -cO * so - sO * co * ci,  sO * si],
        [ sO * co + cO * so * ci, -sO * so + cO * co * ci, -cO * si],
        [             so * si,               co * si,            ci]
    ], dtype=float)

    r_ijk = R @ r_pf
    v_ijk = R @ v_pf
    return r_ijk, v_ijk


# ---------------- Robust Stumpff functions ----------------
def stumpff_C(z: Number) -> float:
    """Numerically stable C(z)."""
    z = float(z)
    az = abs(z)
    if az < 1e-8:
        # C(z) ≈ 1/2 - z/24 + z^2/720
        return float(0.5 - z/24.0 + (z*z)/720.0)
    if z > 0.0:
        s = math.sqrt(z)
        return float((1.0 - math.cos(s)) / z)
    s = math.sqrt(-z)
    return float((math.cosh(s) - 1.0) / (-z))


def stumpff_S(z: Number) -> float:
    """Numerically stable S(z)."""
    z = float(z)
    az = abs(z)
    if az < 1e-8:
        # S(z) ≈ 1/6 - z/120 + z^2/5040
        return float((1.0/6.0) - z/120.0 + (z*z)/5040.0)
    if z > 0.0:
        s = math.sqrt(z)
        return float((s - math.sin(s)) / (s**3))
    s = math.sqrt(-z)
    return float((math.sinh(s) - s) / (s**3))


# ---------------- Lambert (universal variables, single-rev) ----------------
def lambert_universal(r1, r2, tof_s: Number, mu: float = MU_SUN,
                      prograde: bool = True, max_iter: int = 80, tol: float = 1e-8):
    """
    Lambert via universal variables (short-way/long-way via 'prograde').
    r1, r2 in meters; tof_s in seconds. Returns (v1, v2) in m/s.

    Correct formulation (Battin/Vallado):
      y(z) = r1 + r2 + A * ((z*S(z) - 1) / C(z))
      x    = sqrt(y/C)
      t(z) = (x^3*S + A*sqrt(y)) / sqrt(mu)
    """
    r1 = np.asarray(r1, dtype=float).reshape(3)
    r2 = np.asarray(r2, dtype=float).reshape(3)
    r1n = float(np.linalg.norm(r1))
    r2n = float(np.linalg.norm(r2))
    if r1n == 0.0 or r2n == 0.0:
        raise ValueError("Lambert: |r1| or |r2| is zero")
    tof_s = float(tof_s)
    sqrt_mu = math.sqrt(mu)

    # geometry
    cosd = float(np.dot(r1, r2) / (r1n * r2n + 1e-16))
    cosd = max(-1.0, min(1.0, cosd))
    dtheta = math.acos(cosd)

    # short-/long-way selection using out-of-plane direction
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

    def y(z: float) -> float:
        C = stumpff_C(z)
        S = stumpff_S(z)
        Ceff = C if abs(C) > 1e-16 else 1e-16
        return r1n + r2n + A * ((z * S - 1.0) / Ceff)

    def tof(z: float) -> float:
        C = stumpff_C(z)
        S = stumpff_S(z)
        yy = y(z)
        if yy <= 0.0:
            return float('inf')
        Ceff = C if abs(C) > 1e-16 else 1e-16
        x = math.sqrt(max(0.0, yy / Ceff))
        return (x**3 * S + A * math.sqrt(yy)) / sqrt_mu

    # bracket and bisection on z
    z_lo, z_hi = -4.0 * math.pi**2, 4.0 * math.pi**2
    while y(z_lo) <= 0.0:  # nudge until valid
        z_lo += 0.5
        if z_lo > 0.0:
            break
    # expand upper bound until t >= tof
    t_hi = tof(z_hi); tries = 0
    while (not math.isfinite(t_hi)) or (t_hi < tof_s):
        z_hi *= 2.0; tries += 1
        t_hi = tof(z_hi)
        if tries > 40: break

    z = 0.0
    for _ in range(max_iter):
        t_z = tof(z)
        if abs(t_z - tof_s) < tol:
            break
        if t_z < tof_s:
            z_lo = z
        else:
            z_hi = z
        z = 0.5 * (z_lo + z_hi)

    # f, g, gdot and velocities
    yy = y(z)
    if yy <= 0.0:
        raise RuntimeError("Lambert: y(z*) <= 0 at solution")
    f    = 1.0 - yy / r1n
    g    = A * math.sqrt(yy / mu)
    gdot = 1.0 - yy / r2n
    if abs(g) < 1e-14:
        g = math.copysign(1e-14, g)

    v1 = (r2 - f * r1) / g
    v2 = (gdot * r2 - r1) / g
    return v1, v2

# ---------------- Universal Kepler propagator ----------------
def kepler_universal_propagate(r0, v0, dt: Number, mu: float = MU_SUN,
                               itmax: int = 60, tol: float = 1e-9):
    """
    Propagate (r0, v0) by dt seconds via universal variables.
    Returns (r, v). Works for elliptic/hyperbolic/parabolic.
    """
    r0 = np.asarray(r0, dtype=float)
    v0 = np.asarray(v0, dtype=float)

    r0n = float(np.linalg.norm(r0))
    if r0n == 0.0:
        raise ValueError("kepler_universal_propagate: |r0| is zero")

    dt = float(dt)
    vr0   = float(np.dot(r0, v0) / r0n)
    alpha = 2.0 / r0n - float(np.dot(v0, v0)) / mu  # 1/a
    sqrt_mu = math.sqrt(mu)

    # initial guess for chi
    if abs(alpha) > 1e-12:
        chi = math.copysign(1.0, dt) * math.sqrt(mu) * abs(alpha) * dt
    else:
        # near-parabolic guess (Battin-ish)
        h = float(np.linalg.norm(np.cross(r0, v0)))
        p = h * h / mu
        s = 0.5 * (math.pi/2.0 - math.atan(3.0 * math.sqrt(mu / (p**3)) * dt))
        w = math.atan(math.tan(s) ** (1.0/3.0))
        chi = math.sqrt(p) * 2.0 * math.tan(w)

    for _ in range(itmax):
        z  = alpha * chi * chi
        Cz = stumpff_C(z); Sz = stumpff_S(z)
        r  = chi*chi*Cz + vr0*chi*(1.0 - z*Sz) + r0n*(1.0 - z*Cz)
        F  = r - sqrt_mu * dt
        if abs(F) < tol:
            break
        dF = chi * (1.0 - z * Sz)
        chi -= F / (dF + 1e-16)

    z  = alpha * chi * chi
    Cz = stumpff_C(z); Sz = stumpff_S(z)
    f  = 1.0 - (chi*chi*Cz) / r0n
    g  = dt - (chi**3) * Sz / sqrt_mu
    r  = f * r0 + g * v0

    rn = float(np.linalg.norm(r))
    fdot = (sqrt_mu / (rn * r0n)) * (z * Sz - 1.0) * chi
    gdot = 1.0 - (chi*chi*Cz) / rn
    v  = fdot * r0 + gdot * v0
    return r, v


# ---------------- Sampling helper for viewer ----------------
def sample_transfer_polyline(r1, v1, tof_s: Number, n: int = 220, mu: float = MU_SUN):
    """
    Sample points along the transfer from (r1, v1) over tof_s seconds.
    Returns a list of (x_AU, y_AU) for plotting (ecliptic XY projection).
    """
    r1 = np.asarray(r1, dtype=float)
    v1 = np.asarray(v1, dtype=float)
    tof_s = float(tof_s)

    pts = []
    for k in range(n + 1):
        t = (k / n) * tof_s
        r, _ = kepler_universal_propagate(r1, v1, t, mu)
        pts.append((r[0] / AU_M, r[1] / AU_M))
    return pts
