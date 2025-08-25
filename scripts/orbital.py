# scripts/orbital.py
from __future__ import annotations
import math
import numpy as np
from typing import Tuple

AU_M = 149_597_870_700.0
MU_SUN = 1.32712440018e20  # m^3/s^2

def deg2rad(x): return x * math.pi/180.0
def wrap_2pi(x): 
    y = math.fmod(x, 2*math.pi)
    return y + 2*math.pi if y < 0 else y

def kepler_E_from_M(M, e, tol=1e-12, itmax=60):
    M = wrap_2pi(M)
    if e < 1e-12: return M
    E = M if e < 0.8 else math.pi
    for _ in range(itmax):
        f = E - e*math.sin(E) - M
        fp = 1 - e*math.cos(E)
        dE = -f/fp
        E += dE
        if abs(dE) < tol: break
    return E

def oe_to_rv(a_m, e, i, raan, argp, nu, mu=MU_SUN) -> Tuple[np.ndarray,np.ndarray]:
    # Perifocal
    p = a_m*(1 - e*e)
    r_pf = np.array([p*math.cos(nu)/(1+e*math.cos(nu)),
                     p*math.sin(nu)/(1+e*math.cos(nu)),
                     0.0])
    v_pf = np.array([-math.sqrt(mu/p)*math.sin(nu),
                      math.sqrt(mu/p)*(e+math.cos(nu)),
                      0.0])
    ci, si = math.cos(i), math.sin(i)
    cO, sO = math.cos(raan), math.sin(raan)
    co, so = math.cos(argp), math.sin(argp)
    # Rotation matrix PQW->IJK
    R = np.array([
        [ cO*co - sO*so*ci, -cO*so - sO*co*ci,  sO*si],
        [ sO*co + cO*so*ci, -sO*so + cO*co*ci, -cO*si],
        [            so*si,              co*si,    ci]
    ])
    return R @ r_pf, R @ v_pf

def stumpff_C(z):
    if z > 0: 
        s = math.sqrt(z); return (1 - math.cos(s))/z
    if z < 0:
        s = math.sqrt(-z); return (math.cosh(s)-1)/(-z)
    return 0.5

def stumpff_S(z):
    if z > 0:
        s = math.sqrt(z); return (s - math.sin(s))/(s**3)
    if z < 0:
        s = math.sqrt(-z); return (math.sinh(s)-s)/(s**3)
    return 1/6

def lambert_universal(r1, r2, tof_s, mu=MU_SUN, prograde=True, max_iter=80, tol=1e-8):
    """Single-rev universal-variables Lambert (Vallado-like). Returns v1,v2."""
    r1 = np.array(r1, dtype=float); r2 = np.array(r2, dtype=float)
    r1n, r2n = np.linalg.norm(r1), np.linalg.norm(r2)
    cosd = float(np.dot(r1, r2)/(r1n*r2n))
    cross = np.cross(r1, r2)
    dtheta = math.acos(max(-1.0, min(1.0, cosd)))
    if prograde and cross[2] < 0:   # pick the short-way, prograde
        dtheta = 2*math.pi - dtheta
    A = math.sin(dtheta) * math.sqrt(r1n*r2n/(1 - cosd + 1e-15))
    if abs(A) < 1e-12:
        raise ValueError("Geometry singular: A≈0")

    def y(z):
        Cz = stumpff_C(z)
        Sz = stumpff_S(z)
        return r1n + r2n + A*((z*Sz - 1)/math.sqrt(Cz))

    def F(z):
        yz = y(z)
        if yz < 0:  # not physical
            return 1e9
        Cz = stumpff_C(z); Sz = stumpff_S(z)
        term1 = ((yz/Cz)**1.5) * Sz
        term2 = A * math.sqrt(yz)
        return term1 + term2 - math.sqrt(mu)*tof_s

    # Simple bracket + secant
    z0, z1 = -4*math.pi**2, 4*math.pi**2
    f0, f1 = F(z0), F(z1)
    if abs(f0) < tol: z_star = z0
    elif abs(f1) < tol: z_star = z1
    else:
        z_star = 0.0
        for _ in range(max_iter):
            fz = F(z_star)
            if abs(fz) < tol: break
            # choose secant step against closer endpoint
            if abs(f0 - fz) > abs(f1 - fz):
                z_star = z_star - fz*(z1 - z_star)/(f1 - fz + 1e-16)
            else:
                z_star = z_star - fz*(z0 - z_star)/(f0 - fz + 1e-16)
            z_star = float(np.clip(z_star, z0, z1))

    yz = y(z_star); Cz = stumpff_C(z_star); Sz = stumpff_S(z_star)
    f = 1 - yz/r1n
    g = A * math.sqrt(yz/mu)
    gdot = 1 - yz/r2n
    v1 = (r2 - f*r1)/g
    v2 = (gdot*r2 - r1)/g
    return v1, v2

def kepler_universal_propagate(r0, v0, dt, mu=MU_SUN, itmax=60, tol=1e-9):
    """
    Propagate (r0,v0) by dt seconds using universal variables.
    Returns (r, v). Works for elliptic/hyperbolic/parabolic.
    """
    r0 = np.array(r0, dtype=float); v0 = np.array(v0, dtype=float)
    r0n = float(np.linalg.norm(r0))
    if r0n == 0: raise ValueError("r0 norm is zero")
    vr0 = float(np.dot(r0, v0) / r0n)
    alpha = 2.0 / r0n - float(np.dot(v0, v0)) / mu  # 1/a

    # initial guess for chi
    if abs(alpha) > 1e-12:
        chi = np.sign(dt) * np.sqrt(mu) * abs(alpha) * dt
    else:
        # near-parabolic
        h = np.linalg.norm(np.cross(r0, v0))
        p = h*h / mu
        s = 0.5 * (math.pi/2 - math.atan(3*np.sqrt(mu/(p**3))*dt))
        w = math.atan(math.tan(s)**(1/3))
        chi = np.sqrt(p) * 2 * math.tan(w)

    for _ in range(itmax):
        z  = alpha * chi * chi
        Cz = stumpff_C(z); Sz = stumpff_S(z)
        r  = chi*chi*Cz + vr0*chi*(1 - z*Sz) + r0n*(1 - z*Cz)
        F  = r - np.sqrt(mu) * dt
        if abs(F) < tol: break
        dF = chi*(1 - z*Sz)
        chi -= F / (dF + 1e-16)

    z  = alpha * chi * chi
    Cz = stumpff_C(z); Sz = stumpff_S(z)
    f  = 1 - (chi*chi*Cz) / r0n
    g  = dt - (chi**3)*Sz / np.sqrt(mu)
    r  = f*r0 + g*v0
    rn = float(np.linalg.norm(r))
    fdot = (np.sqrt(mu)/(rn*r0n)) * (z*Sz - 1) * chi
    gdot = 1 - (chi*chi*Cz) / rn
    v  = fdot*r0 + gdot*v0
    return r, v

def sample_transfer_polyline(r1, v1, tof_s, n=220, mu=MU_SUN):
    """
    Sample points along the transfer from (r1,v1) over tof_s seconds.
    Returns a list of (x_AU, y_AU) for plotting (ecliptic XY projection).
    """
    pts = []
    for k in range(n+1):
        t = (k/n) * float(tof_s)
        r, _ = kepler_universal_propagate(r1, v1, t, mu)
        pts.append((r[0]/AU_M, r[1]/AU_M))
    return pts