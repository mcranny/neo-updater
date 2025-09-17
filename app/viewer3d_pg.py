#!/usr/bin/env python3
# app/viewer3d_pg.py
from __future__ import annotations
import json, math, sys, time
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg
import pyqtgraph.opengl as gl

# ---------- constants ----------
AU = 1.0
AU_M = 149_597_870_700.0
MU_SUN = 1.32712440018e20  # m^3/s^2
DEG = math.pi / 180.0
DEFAULT_JSON = Path("data/hazardous_neos/latest_intercept.json")

# ---------- small kepler helpers (viewer-local; no extra imports) ----------
def _kepler_E_from_M(M: float, e: float, tol=1e-12, maxit=60) -> float:
    M = (M + math.pi) % (2*math.pi) - math.pi
    E = M if e < 0.8 else math.pi
    for _ in range(maxit):
        f = E - e*math.sin(E) - M
        fp = 1.0 - e*math.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            break
    return E

def _rv_from_elements_at_time(a_AU, e, i_deg, raan_deg, argp_deg, ma0_deg, epoch0_jd, t_jd) -> Tuple[np.ndarray, np.ndarray]:
    a = float(a_AU) * AU_M
    e = float(e)
    i = float(i_deg) * DEG; Om = float(raan_deg) * DEG; w = float(argp_deg) * DEG; M0 = float(ma0_deg) * DEG
    n = math.sqrt(MU_SUN / (a**3))
    dt = (float(t_jd) - float(epoch0_jd)) * 86400.0
    M = M0 + n * dt
    E = _kepler_E_from_M(M, e)

    cosE = math.cos(E); sinE = math.sin(E)
    fac  = math.sqrt(max(1e-16, 1 - e*e))
    r_p  = a * (1 - e*cosE)
    x_p  = a * (cosE - e)
    y_p  = a * fac * sinE
    rdot   = (math.sqrt(MU_SUN * a) / max(1e-12, r_p)) * e * sinE
    rfdot  = (math.sqrt(MU_SUN * a) / max(1e-12, r_p)) * fac * cosE

    cw, sw = math.cos(w), math.sin(w)
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(i),  math.sin(i)

    R11 = cO*cw - sO*sw*ci; R12 = -cO*sw - sO*cw*ci; R13 = sO*si
    R21 = sO*cw + cO*sw*ci; R22 = -sO*sw + cO*cw*ci; R23 = -cO*si
    R31 = sw*si;             R32 = cw*si;             R33 = ci
    R   = np.array([[R11,R12,R13],[R21,R22,R23],[R31,R32,R33]], dtype=float)

    r_pf = np.array([x_p, y_p, 0.0])
    v_pf = np.array([rdot, rfdot, 0.0])
    r = R @ r_pf
    v = R @ v_pf
    return r, v

def _ellipse_poly_from_elements(elements: Dict[str, Any], samples=512) -> np.ndarray:
    """Generate a static orbit curve in AU using osculating elements (orientation is correct)."""
    a = float(elements.get("a_AU", 1.0)) * AU_M
    e = float(elements.get("e", 0.0))
    i = float(elements.get("i_deg", 0.0)) * DEG
    Om = float(elements.get("raan_deg", 0.0)) * DEG
    w  = float(elements.get("argp_deg", 0.0)) * DEG

    # sweep true anomaly ν
    nu = np.linspace(0.0, 2.0*math.pi, int(max(64, samples)), endpoint=True)
    p = a * (1.0 - e*e)
    cnu, snu = np.cos(nu), np.sin(nu)
    r_pf = np.vstack([p * cnu / (1.0 + e * cnu),
                      p * snu / (1.0 + e * cnu),
                      np.zeros_like(nu)])
    # PQW->IJK
    ci, si = math.cos(i), math.sin(i)
    cO, sO = math.cos(Om), math.sin(Om)
    co, so = math.cos(w),  math.sin(w)
    R = np.array([
        [ cO*co - sO*so*ci, -cO*so - sO*co*ci,  sO*si],
        [ sO*co + cO*so*ci, -sO*so + cO*co*ci, -cO*si],
        [            so*si,             co*si,      ci]
    ], dtype=float)
    r_ijk = (R @ r_pf).T / AU_M
    return r_ijk.astype(np.float32)

# ---------- time helpers ----------
def jd_from_unix(t: float) -> float:
    return t/86400.0 + 2440587.5

# ---------- loader ----------
def load_plan(path: Path) -> Dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    items = data.get("potentially_hazardous_neos") or data.get("neos") or []
    if not items:
        raise RuntimeError("No NEOs found in file.")
    neo = items[0]
    plan = dict(neo.get("intercept_plan") or {})
    plan["_neo_name"] = neo.get("name", "Target")
    return plan

# ---------- main widget ----------
class Viewer(QtWidgets.QWidget):
    def __init__(self, plan: Dict[str, Any]):
        super().__init__()
        self.plan = plan
        self.setWindowTitle("NEO Updater — 3D viewer")

        lay = QtWidgets.QVBoxLayout(self)
        self.view = gl.GLViewWidget()
        self.view.setBackgroundColor((14,16,20))
        self.view.opts["distance"] = 6.0
        lay.addWidget(self.view, 1)

        # HUD
        self.hud = QtWidgets.QLabel("")
        self.hud.setStyleSheet("color:#ddd; font: 12px 'Inter', 'Helvetica'; padding:4px")
        lay.addWidget(self.hud, 0)

        # axes/grid
        gxy = gl.GLGridItem()
        gxz = gl.GLGridItem(); gxz.rotate(90, 1, 0, 0)
        gyz = gl.GLGridItem(); gyz.rotate(90, 0, 1, 0)
        for g in (gxy, gxz, gyz):
            g.setSize(8, 8); g.setSpacing(1, 1); g.setDepthValue(1e6)
            g.setColor((0.5,0.5,0.6,0.25))
            self.view.addItem(g)

        # Sun
        self.sun = gl.GLScatterPlotItem(pos=np.array([[0,0,0]], dtype=float),
                                        size=12, color=(1.0,0.85,0.4,1.0), pxMode=True)
        self.view.addItem(self.sun)

        # Orbits (static curves for context)
        earth_el = plan.get("elements_earth") or {}
        targ_el  = plan.get("elements_target") or plan.get("elements") or {}
        earth_curve = _ellipse_poly_from_elements(earth_el, samples=1024)
        targ_curve  = _ellipse_poly_from_elements(targ_el,  samples=1024)

        self.earth_ring = gl.GLLinePlotItem(pos=earth_curve, mode='line_strip', antialias=True, width=2,
                                            color=(0.25,0.65,1.0,1.0))
        self.targ_ring  = gl.GLLinePlotItem(pos=targ_curve,  mode='line_strip', antialias=True, width=2,
                                            color=(0.65,1.0,0.55,1.0))
        self.view.addItem(self.earth_ring)
        self.view.addItem(self.targ_ring)

        # Lambert polyline (CURVED) — segment to avoid huge jumps
        lam3 = plan.get("lambert_poly_xyz_au") or plan.get("lambert_polyline_xyz_au")
        lam2 = plan.get("lambert_polyline_xy_au")
        lam = None
        if lam3:
            lam = np.asarray(lam3, dtype=float)
            if lam.ndim == 2 and lam.shape[1] == 2:
                lam = np.pad(lam, ((0,0),(0,1)), mode='constant')
        elif lam2:
            lam = np.pad(np.asarray(lam2, dtype=float), ((0,0),(0,1)), mode='constant')

        self.transfer_segs = []
        self.lambert = None
        self.num_lam = 0
        if lam is not None and lam.shape[0] >= 2:
            # split when successive points jump more than 2 AU
            jumps = np.linalg.norm(np.diff(lam, axis=0), axis=1)
            split_idx = np.where(jumps > 2.0)[0]
            start = 0
            for idx in list(split_idx) + [lam.shape[0]-1]:
                seg = lam[start:idx+1]
                if seg.shape[0] >= 2:
                    item = gl.GLLinePlotItem(pos=seg.astype(np.float32),
                                             mode='line_strip', antialias=True, width=3,
                                             color=(1.0,0.82,0.45,1.0))
                    self.view.addItem(item)
                    self.transfer_segs.append(item)
                start = idx+1
            self.lambert = lam.astype(np.float32)
            self.num_lam = int(self.lambert.shape[0])

        # Moving markers
        self.earth_dot = gl.GLScatterPlotItem(size=8, color=(0.25,0.65,1.0,1.0), pxMode=True)
        self.targ_dot  = gl.GLScatterPlotItem(size=8, color=(0.65,1.0,0.55,1.0), pxMode=True)
        self.sc_dot    = gl.GLScatterPlotItem(size=8, color=(1.0,0.65,0.35,1.0), pxMode=True)
        self.view.addItem(self.earth_dot)
        self.view.addItem(self.targ_dot)
        self.view.addItem(self.sc_dot)

        # Dep/Arr markers (first/last Lambert points)
        r1 = np.array(plan.get("r1_au") or [0,0,0], dtype=float)
        r2 = np.array(plan.get("r2_au") or [0,0,0], dtype=float)
        self.dep = gl.GLScatterPlotItem(pos=r1[None,:], size=10, color=(0.2,0.8,1.0,1.0), pxMode=True)
        self.arr = gl.GLScatterPlotItem(pos=r2[None,:], size=10, color=(0.6,1.0,0.6,1.0), pxMode=True)
        self.view.addItem(self.dep); self.view.addItem(self.arr)

        # Time bookkeeping
        self.depart_jd = float(plan.get("departure_jd", jd_from_unix(time.time())))
        self.arrive_jd = float(plan.get("arrival_jd", self.depart_jd + float(plan.get("tof_days", 180.0))))
        self.tof_days  = float(plan.get("tof_days", max(1.0, self.arrive_jd - self.depart_jd)))
        self.t_days    = 0.0
        self.playing   = True
        self.speed     = 1.0
        self.loop      = True  # loop after TOF so Earth/NEO keep moving

        # Cache elements (used in propagation)
        self.el_earth = earth_el
        self.el_targ  = targ_el

        # timer
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._tick)
        self.timer.start(16)

    # ---- animation step ----
    def _tick(self):
        if not self.playing:
            return
        # advance time
        step = max(1e-4, self.tof_days/1200.0) * self.speed
        self.t_days += step
        if self.t_days >= self.tof_days:
            if self.loop:
                self.t_days = 0.0
            else:
                self.t_days = self.tof_days
                self.playing = False
        frac = self.t_days / max(1e-9, self.tof_days)
        jd = self.depart_jd + self.t_days

        # Earth / Target positions (AU)
        rE, _ = _rv_from_elements_at_time(self.el_earth["a_AU"], self.el_earth["e"], self.el_earth["i_deg"],
                                          self.el_earth["raan_deg"], self.el_earth["argp_deg"], self.el_earth["ma_deg"],
                                          self.el_earth["epoch_jd"], jd)
        rT, _ = _rv_from_elements_at_time(self.el_targ["a_AU"], self.el_targ["e"], self.el_targ["i_deg"],
                                          self.el_targ["raan_deg"], self.el_targ["argp_deg"], self.el_targ["ma_deg"],
                                          self.el_targ["epoch_jd"], jd)
        rE /= AU_M; rT /= AU_M
        self.earth_dot.setData(pos=rE.reshape(1,3))
        self.targ_dot.setData(pos=rT.reshape(1,3))

        # Spacecraft along Lambert
        if self.lambert is not None and self.num_lam >= 2:
            idx = int(frac * (self.num_lam - 1))
            idx = min(self.num_lam - 1, max(0, idx))
            self.sc_dot.setData(pos=self.lambert[idx:idx+1])
            self.hud.setText(f"{self.plan.get('_neo_name','Target')}  —  t = {self.t_days:6.1f} / {self.tof_days:.1f} days   (speed x{self.speed:.1f}, loop={'on' if self.loop else 'off'})")
        else:
            self.sc_dot.setData(pos=np.zeros((1,3), dtype=float))
            self.hud.setText("No Lambert polyline found in plan.")

    # ---- keys ----
    def keyPressEvent(self, ev: QtGui.QKeyEvent) -> None:
        k = ev.key()
        if k == QtCore.Qt.Key_Space:
            self.playing = not self.playing
        elif k in (QtCore.Qt.Key_Plus, QtCore.Qt.Key_Equal):
            self.speed = min(20.0, self.speed + 0.5)
        elif k in (QtCore.Qt.Key_Minus, QtCore.Qt.Key_Underscore):
            self.speed = max(0.5, self.speed - 0.5)
        elif k == QtCore.Qt.Key_0:
            self.t_days = 0.0; self.playing = True
        elif k == QtCore.Qt.Key_L:
            self.loop = not self.loop
        ev.accept()

# ---------- entry ----------
def main():
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_JSON
    plan = load_plan(path)
    app = QtWidgets.QApplication(sys.argv)
    w = Viewer(plan)
    w.resize(1100, 760)
    w.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
