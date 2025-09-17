#!/usr/bin/env python3
# pip install PySide6 pyqtgraph numpy PyOpenGL
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np
from PySide6 import QtCore, QtWidgets
from pyqtgraph.opengl import GLViewWidget, GLLinePlotItem, GLScatterPlotItem, GLMeshItem, MeshData

AU_M  = 149_597_870_700.0
MU_SUN = 1.32712440018e20
DEG = math.pi / 180.0
DEFAULT_JSON = Path("data/hazardous_neos/latest_intercept.json")

def load_first_plan(path: Path) -> dict:
    blob = json.loads(path.read_text(encoding="utf-8"))
    items = blob.get("potentially_hazardous_neos") or blob.get("neos") or []
    if not items: raise RuntimeError("No NEOs found in JSON.")
    neo = items[0]
    plan = neo.get("intercept_plan") or {}
    plan["_neo_name"] = neo.get("name", "Target")
    return plan

def kepler_E_from_M(M: float, e: float, tol=1e-12, maxit=50) -> float:
    M = (M + math.pi) % (2*math.pi) - math.pi
    E = M if e < 0.8 else math.pi
    for _ in range(maxit):
        f  = E - e*math.sin(E) - M
        fp = 1 - e*math.cos(E)
        dE = -f/fp
        E += dE
        if abs(dE) < tol: break
    return E

def rv_from_elements_at_time(a_au, e, i_deg, raan_deg, argp_deg, ma0_deg, epoch0_jd, t_jd, mu=MU_SUN):
    a = a_au * AU_M
    i = i_deg * DEG; Om = raan_deg * DEG; w = argp_deg * DEG; M0 = ma0_deg * DEG
    n = math.sqrt(mu / (a**3))
    dt = (t_jd - epoch0_jd) * 86400.0
    M = M0 + n * dt
    E = kepler_E_from_M(M, e)
    cosE, sinE = math.cos(E), math.sin(E)
    r_p = a * (1 - e*cosE)
    x_p = a * (cosE - e)
    y_p = a * math.sqrt(1 - e*e) * sinE
    cw, sw = math.cos(w), math.sin(w)
    cO, sO = math.cos(Om), math.sin(Om)
    ci, si = math.cos(i), math.sin(i)
    x1 = cw*x_p - sw*y_p
    y1 = sw*x_p + cw*y_p
    x2 = x1
    y2 =  ci*y1
    z2 =  si*y1
    x =  cO*x2 - sO*y2
    y =  sO*x2 + cO*y2
    z =  z2
    return np.array([x, y, z]), None

def sample_orbit_xyz(elements: dict, jd0: float, days: float, n: int) -> np.ndarray:
    t = np.linspace(jd0, jd0 + days, n)
    pts = []
    for tj in t:
        r,_ = rv_from_elements_at_time(
            elements["a_AU"], elements["e"], elements["i_deg"],
            elements["raan_deg"], elements["argp_deg"], elements["ma_deg"],
            elements["epoch_jd"], tj
        )
        pts.append((r / AU_M))
    return np.vstack(pts)

def line_item(points_au: np.ndarray, color=(1,1,1,1), width=2.0) -> GLLinePlotItem:
    return GLLinePlotItem(pos=points_au.astype(np.float32), color=color, width=width, antialias=True, mode='line_strip')

def sphere_item(radius_au: float, color=(1,1,0,1)) -> GLMeshItem:
    md = MeshData.sphere(rows=24, cols=24, radius=radius_au)
    return GLMeshItem(meshdata=md, smooth=True, color=color, shader='shaded', drawEdges=False)

class OrbitViewer3D(QtWidgets.QMainWindow):
    def __init__(self, json_path: Path = DEFAULT_JSON):
        super().__init__()
        self.setWindowTitle("NEO Updater — 3D Viewer (Pure Python)")
        self.resize(1200, 800)
        self.view = GLViewWidget()
        self.view.opts['distance'] = 5.0
        self.view.setBackgroundColor((14,15,19))
        cw = QtWidgets.QWidget(); lay = QtWidgets.QVBoxLayout(cw); lay.setContentsMargins(0,0,0,0)
        lay.addWidget(self.view); self.setCentralWidget(cw)
        self.hud = QtWidgets.QLabel(""); self.hud.setStyleSheet("color:#ddd; padding:6px;"); lay.addWidget(self.hud)

        self.plan = load_first_plan(json_path)

        self.earth = self.plan.get("elements_earth") or {
            "a_AU": 1.00000011, "e": 0.01671022, "i_deg": 0.00005,
            "raan_deg": -11.26064, "argp_deg": 102.94719,
            "ma_deg": 100.46435, "epoch_jd": 2451545.0
        }
        self.target = self.plan.get("elements_target") or self.plan.get("elements") or self.plan.get("sbdb_elements")

        dep_jd = float(self.plan.get("departure_jd") or self.plan.get("depart_epoch_jd") or 2451545.0)
        self.earth_curve  = sample_orbit_xyz(self.earth,  dep_jd - 365, 730, 800)

        self.target_curve = None
        if self.target:
            try:
                self.target_curve = sample_orbit_xyz(self.target, dep_jd - 365, 730, 800)
            except Exception:
                self.target_curve = None  # tolerate bad/partial elements

        lambert = (self.plan.get("lambert_poly_xyz_au")
                   or self.plan.get("lambert_polyline_xyz_au")
                   or self.plan.get("lambert_polyline_xy_au"))
        if lambert and len(lambert) and len(lambert[0]) == 2:
            lambert = [[x, y, 0.0] for x, y in lambert]
        self.transfer = np.array(lambert, dtype=np.float32) if lambert else None

        self._add_scene()
        self.sc_idx = 0
        self.timer = QtCore.QTimer(self); self.timer.timeout.connect(self._tick); self.timer.start(16)

    def _add_scene(self):
        self.view.addItem(sphere_item(0.05, (1,0.84,0.4,1)))
        self.view.addItem(line_item(self.earth_curve,  (0.25,0.75,1,1), 2.2))
        if self.target_curve is not None:
            self.view.addItem(line_item(self.target_curve, (0.65,1,0.5,1), 2.2))
        if self.transfer is not None and len(self.transfer) >= 2:
            self.view.addItem(line_item(self.transfer, (0.96,0.76,0.47,1), 2.5))
        self.sc = GLScatterPlotItem(pos=np.array([[0,0,0]], dtype=np.float32), size=8, color=(0.95,0.55,0.35,1), pxMode=True)
        self.view.addItem(self.sc)
        for ax in [(2,0,0),(-2,0,0),(0,2,0),(0,-2,0),(0,0,2),(0,0,-2)]:
            self.view.addItem(GLLinePlotItem(pos=np.array([[0,0,0], ax], dtype=np.float32), color=(0.6,0.6,0.7,0.7), width=1.0, antialias=True, mode='line_strip'))

    def _tick(self):
        if self.transfer is None or len(self.transfer) < 2:
            self.hud.setText(f"{self.plan.get('_neo_name')} — no transfer polyline found.")
            return
        self.sc_idx = (self.sc_idx + 1) % len(self.transfer)
        self.sc.setData(pos=self.transfer[self.sc_idx:self.sc_idx+1])
        frac = self.sc_idx / max(1, len(self.transfer)-1)
        self.hud.setText(f"{self.plan.get('_neo_name')} — along transfer: {frac*100:.1f}%")

def main():
    app = QtWidgets.QApplication(sys.argv)
    w = OrbitViewer3D(DEFAULT_JSON)
    w.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
