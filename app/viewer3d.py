from __future__ import annotations
"""
NEO Intercept 3D Viewer (pyqtgraph + OpenGL)
- Demo-friendly and tolerant of field-name differences
- Converts meters→AU automatically
- Draws axes, 1 AU ring, Earth (blue) & target (green) orbits
- Draws transfer (orange) and animates a probe along it
- NEW: time-synced animation — Earth/target move with a simple ephemeris
  and the probe reaches r2 at the same arrival time

Usage:
    python app/viewer3d.py path/to/intercept.json
"""

import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple

import numpy as np

# Qt / pyqtgraph
from PySide6 import QtCore, QtGui, QtWidgets
import pyqtgraph as pg
import pyqtgraph.opengl as gl

# ---------- Constants ----------
AU_M = 149_597_870_700.0  # meters per AU
DEG2RAD = math.pi / 180.0
SIDEREAL_YEAR_DAYS = 365.256363004
MU_SUN = 1.32712440018e20  # m^3/s^2

# ---------- Unit & drawing helpers ----------
def _as_au(pts: Iterable[Iterable[float]] | None) -> np.ndarray:
    if pts is None:
        return np.zeros((0, 3), dtype=np.float32)
    try:
        a = np.asarray(pts, dtype=np.float64)
    except Exception:
        return np.zeros((0, 3), dtype=np.float32)
    if a.ndim == 1:
        a = a.reshape(1, -1)
    if a.size == 0:
        return np.zeros((0, 3), dtype=np.float32)
    if a.shape[1] == 2:
        a = np.c_[a, np.zeros(len(a))]
    with np.errstate(over="ignore", invalid="ignore"):
        mx = float(np.nanmax(np.abs(a)))
    if np.isfinite(mx) and mx > 10.0:  # likely meters
        a = a / AU_M
    return a.astype(np.float32)


def make_line(
    pts: Iterable[Iterable[float]],
    width: float = 2.0,
    color: Tuple[float, float, float, float] = (1, 1, 1, 1),
) -> gl.GLLinePlotItem:
    arr = _as_au(pts)
    return gl.GLLinePlotItem(pos=arr, width=width, color=color, antialias=True, mode="line_strip")


def make_scatter(
    pts: Iterable[Iterable[float]] | np.ndarray,
    size: float = 8.0,
    color: Tuple[float, float, float, float] = (1, 1, 1, 1),
) -> gl.GLScatterPlotItem:
    arr = _as_au(pts)
    if arr.ndim == 1:
        arr = arr.reshape(1, 3)
    return gl.GLScatterPlotItem(pos=arr, size=size, color=color, pxMode=True)


def circle_orbit_points(R_au: float, n: int = 512) -> np.ndarray:
    th = np.linspace(0.0, 2.0 * np.pi, n, endpoint=True, dtype=np.float64)
    x = R_au * np.cos(th)
    y = R_au * np.sin(th)
    z = np.zeros_like(x)
    return np.c_[x, y, z].astype(np.float32)


def rot3(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]], dtype=np.float64)


def rot1(theta: float) -> np.ndarray:
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]], dtype=np.float64)


def neo_orbit_points_from_elements(elems: Dict[str, Any], n: int = 512) -> np.ndarray:
    a_au = float(elems.get("a_AU") or elems.get("a_au") or elems.get("a") or 1.2)
    e = float(elems.get("e", 0.0))
    i = float(elems.get("i_deg") or elems.get("i") or 0.0) * DEG2RAD
    raan = float(elems.get("raan_deg") or elems.get("RAAN_deg") or elems.get("raan") or 0.0) * DEG2RAD
    argp = float(elems.get("argp_deg") or elems.get("omega_deg") or elems.get("argp") or 0.0) * DEG2RAD
    nu = np.linspace(0.0, 2.0 * np.pi, n, dtype=np.float64)
    p = a_au * (1.0 - e * e)
    r = p / (1.0 + e * np.cos(nu))
    xp = r * np.cos(nu)
    yp = r * np.sin(nu)
    zp = np.zeros_like(xp)
    rp = np.vstack((xp, yp, zp))
    R = rot3(raan) @ rot1(i) @ rot3(argp)
    r_eci = (R @ rp).T
    return r_eci.astype(np.float32)


def _get_ci(d: dict, *names):
    if not isinstance(d, dict):
        return None
    low = {k.lower(): k for k in d.keys()}
    for name in names:
        k = low.get(str(name).lower())
        if k is not None:
            return d[k]
    return None

# ---------- Simple ephemerides ----------
def earth_pos_au_from_jd(jd: float) -> np.ndarray:
    th = 2.0 * np.pi * ((jd - 2451545.0) / SIDEREAL_YEAR_DAYS)
    return np.array([np.cos(th), np.sin(th), 0.0], dtype=np.float32)


def _kepler_E_from_M(M: float, e: float, tol=1e-12, maxit=50) -> float:
    M = (M + math.pi) % (2 * math.pi) - math.pi
    E = M if e < 0.8 else math.pi
    for _ in range(maxit):
        f = E - e * math.sin(E) - M
        fp = 1 - e * math.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            break
    return E


def target_pos_au_from_elements_at_jd(elems: Dict[str, Any], t_jd: float) -> np.ndarray:
    a_au = float(elems.get("a_AU") or elems.get("a_au") or elems.get("a") or 1.2)
    e = float(elems.get("e", 0.0))
    i = float(elems.get("i_deg") or elems.get("i") or 0.0) * DEG2RAD
    raan = float(elems.get("raan_deg") or elems.get("RAAN_deg") or elems.get("raan") or 0.0) * DEG2RAD
    argp = float(elems.get("argp_deg") or elems.get("omega_deg") or elems.get("argp") or 0.0) * DEG2RAD
    ma0_deg = float(elems.get("ma_deg", 0.0))
    epoch0_jd = float(elems.get("epoch_jd", t_jd))

    a_m = a_au * AU_M
    n = math.sqrt(MU_SUN / (a_m ** 3))  # rad/s
    dt = (t_jd - epoch0_jd) * 86400.0
    M = (ma0_deg * DEG2RAD) + n * dt
    E = _kepler_E_from_M(M, e)
    cosE, sinE = math.cos(E), math.sin(E)
    r_p = a_m * (1 - e * cosE)
    x_p = a_m * (cosE - e)
    y_p = a_m * (math.sqrt(1 - e * e) * sinE)

    cO, sO = math.cos(raan), math.sin(raan)
    ci, si = math.cos(i), math.sin(i)
    cw, sw = math.cos(argp), math.sin(argp)
    R = np.array(
        [
            [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
            [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
            [sw * si, cw * si, ci],
        ],
        dtype=float,
    )
    r_pf = np.array([x_p, y_p, 0.0])
    r = (R @ r_pf) / AU_M
    return r.astype(np.float32)


# ---------- Data loading ----------
def load_intercept(path: Path) -> Dict[str, Any]:
    with path.open("r") as f:
        j = json.load(f)
    arr = j.get("potentially_hazardous_neos") or j.get("neos") or []
    if not arr:
        return {}
    row = arr[0]
    return row.get("intercept_plan", {})


def extract_polyline(plan: Dict[str, Any]) -> Optional[np.ndarray]:
    pts = _get_ci(plan, "lambert_polyline_xyz_au", "lambert_polyline_xy_au", "lambert_polyline", "lambert_polyline_m")
    if pts is None:
        r1 = _get_ci(plan, "r1_au", "r1_AU", "r1", "r_depart", "r_depart_au", "r_depart_m")
        r2 = _get_ci(plan, "r2_au", "r2_AU", "r2", "r_arrive", "r_arrive_au", "r_arrive_m")
        if r1 is not None and r2 is not None:
            return _as_au([r1, r2])
        return None
    return _as_au(pts)


def extract_r1_r2(plan: Dict[str, Any]) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    r1 = _get_ci(plan, "r1_au", "r1_AU", "r1", "r_depart", "r_depart_au", "r_depart_m")
    r2 = _get_ci(plan, "r2_au", "r2_AU", "r2", "r_arrive", "r_arrive_au", "r_arrive_m")
    return (_as_au(r1) if r1 is not None else None, _as_au(r2) if r2 is not None else None)


# ---------- Main window ----------
class Viewer3D(QtWidgets.QMainWindow):
    def __init__(self, json_path: Path):
        super().__init__()
        self.setWindowTitle("NEO Intercept 3D Viewer")

        # GL view
        self.view = gl.GLViewWidget()
        self.view.setMinimumSize(900, 700)
        try:
            self.view.setBackgroundColor(QtGui.QColor("#0e0f13"))
        except Exception:
            pass
        self.setCentralWidget(self.view)

        # Axes + 1 AU reference ring
        ax = gl.GLAxisItem(); ax.setSize(1, 1, 1)
        self.view.addItem(ax)
        ring = circle_orbit_points(1.0, 256)
        self.view.addItem(gl.GLLinePlotItem(pos=ring, width=1.2, color=(0.2, 0.6, 1.0, 0.35), antialias=True, mode="line_strip"))
        self.view.setCameraPosition(distance=3.2, azimuth=45, elevation=22)

        # State
        self.items: Dict[str, Any] = {}
        self.polyline: Optional[np.ndarray] = None
        self.anim_idx: int = 0

        # Simulation time (days)
        plan = load_intercept(json_path)
        self.depart_jd = float(plan.get("depart_epoch_jd", 2460000.5))
        self.tof_days = float(plan.get("tof_days", 180.0))
        self.arrive_jd = float(plan.get("arrive_epoch_jd", self.depart_jd + self.tof_days))
        self.sim_jd = self.depart_jd
        self.sim_days_per_tick = float(plan.get("sim_days_per_tick", 2.0))  # speed

        # Build scene
        self._build_scene(plan)

        # Timer
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._tick)
        self.timer.setInterval(25)
        self.timer.start()

        # Toolbar
        tb = self.addToolBar("Playback")
        act_play = QtGui.QAction("Play", self)
        act_stop = QtGui.QAction("Pause", self)
        act_play.triggered.connect(lambda: self.timer.start())
        act_stop.triggered.connect(lambda: self.timer.stop())
        tb.addAction(act_play)
        tb.addAction(act_stop)

    def _build_scene(self, plan: Dict[str, Any]) -> None:
        # Earth orbit + moving dot driven by sim time
        self.earth_curve = circle_orbit_points(1.0, 720)
        self.view.addItem(make_line(self.earth_curve, width=1.5, color=(0.4, 0.7, 1.0, 0.9)))
        self.items["earth_dot"] = make_scatter([[1.0, 0.0, 0.0]], size=7, color=(0.35, 0.8, 1.0, 1.0))
        self.view.addItem(self.items["earth_dot"])

        # Target orbit curve for display
        self.elems = plan.get("elements") or {}
        if self.elems:
            tgt_pts = neo_orbit_points_from_elements(self.elems, n=720)
        else:
            r2 = _get_ci(plan, "r2_au", "r2_AU", "r2")
            R = float(np.linalg.norm(_as_au(r2)[0])) if r2 is not None else 1.2
            tgt_pts = circle_orbit_points(R, 720)
        self.view.addItem(make_line(tgt_pts, width=2.0, color=(0.6, 1.0, 0.3, 0.9)))
        # Target moving dot (position computed from elements or along curve by phase)
        self.items["target_dot"] = make_scatter(tgt_pts[0], size=7, color=(0.8, 1.0, 0.5, 1.0))
        self.view.addItem(self.items["target_dot"])

        # Transfer polyline
        self.polyline = extract_polyline(plan)
        if self.polyline is not None and len(self.polyline) >= 2:
            self.view.addItem(make_line(self.polyline, width=2.0, color=(1.0, 0.7, 0.2, 0.95)))
            # Probe marker
            self.items["probe"] = make_scatter(self.polyline[0], size=9, color=(1.0, 0.85, 0.3, 1.0))
            self.view.addItem(self.items["probe"])

        # r1 / r2 markers
        r1, r2 = extract_r1_r2(plan)
        if r1 is not None and len(r1) > 0:
            self.view.addItem(make_scatter(r1[0], size=9, color=(0.2, 0.8, 1.0, 1.0)))
        if r2 is not None and len(r2) > 0:
            self.view.addItem(make_scatter(r2[0], size=9, color=(0.9, 0.3, 0.3, 1.0)))

    # ---------- animation tick ----------
    def _tick(self) -> None:
        # advance sim time
        self.sim_jd += self.sim_days_per_tick
        if self.sim_jd > self.arrive_jd:
            # loop
            self.sim_jd = self.depart_jd

        # phase in [0,1]
        phase = (self.sim_jd - self.depart_jd) / max(self.tof_days, 1e-6)
        phase = float(np.clip(phase, 0.0, 1.0))

        # Earth dot (simple circular ephemeris)
        epos = earth_pos_au_from_jd(self.sim_jd)
        edot = self.items.get("earth_dot")
        if isinstance(edot, gl.GLScatterPlotItem):
            edot.setData(pos=epos.reshape(1, 3))

        # Target dot (physics if elements available; otherwise move by phase along curve)
        tdot = self.items.get("target_dot")
        if isinstance(tdot, gl.GLScatterPlotItem):
            if self.elems:
                tpos = target_pos_au_from_elements_at_jd(self.elems, self.sim_jd)
            else:
                # no elements; approximate along displayed curve by phase
                curve = getattr(self, "earth_curve", None)
                tpos = curve[int(phase * (len(curve) - 1))] if curve is not None else np.zeros(3)
            tdot.setData(pos=tpos.reshape(1, 3))

        # Probe along polyline synchronized to phase
        if self.polyline is not None and len(self.polyline) >= 2:
            idx = int(round(phase * (len(self.polyline) - 1)))
            idx = max(0, min(idx, len(self.polyline) - 1))
            ppos = self.polyline[idx]
            probe = self.items.get("probe")
            if isinstance(probe, gl.GLScatterPlotItem):
                probe.setData(pos=ppos.reshape(1, 3))


# ---------- Entrypoint ----------
def main(argv: List[str]) -> int:
    if len(argv) < 2:
        print("Usage: python app/viewer3d.py <path/to/intercept.json>")
        return 2
    json_path = Path(argv[1]).expanduser().resolve()
    if not json_path.exists():
        print(f"Error: file not found: {json_path}")
        return 2
    app = QtWidgets.QApplication(sys.argv)
    win = Viewer3D(json_path)
    win.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
