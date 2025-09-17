#!/usr/bin/env python3
"""
NEO Updater — Qt 2-D Viewer (Lambert-aware + orientation fix, Qt6 enum-safe)

- Accepts 3-D or 2-D Lambert polylines from the planner
- Correctly orients the target orbit by applying Δ(Ω+ω) between Earth and target
  when Lambert is NOT provided (so orbits no longer look “same-angled”)
- With Lambert, anchors Earth to departure heading and phases target to reach
  the arrival heading at TOF

Usage:
    python app/neo_viewer_qt.py  [optional_json_path]
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
from typing import Any

from PySide6.QtCore import Qt, QTimer, QRectF, QPointF
from PySide6.QtGui import (
    QPalette, QColor, QBrush, QPen, QFont, QPainterPath, QPainter
)
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QPushButton, QSlider, QGraphicsView, QGraphicsScene, QGraphicsEllipseItem,
    QGraphicsPathItem
)

AU = 1.0
DEFAULT_INTERCEPT = Path("data/hazardous_neos/latest_intercept.json")

def mean_motion(a_AU: float) -> float:
    """rad/year for circularized 2-D rings (relative timing only)."""
    return (max(0.05, a_AU) ** -1.5)

def load_first_plan(path: Path) -> dict:
    blob = json.loads(path.read_text(encoding="utf-8"))
    items = blob.get("potentially_hazardous_neos") or blob.get("neos") or []
    if not items: raise RuntimeError("No NEOs found in JSON.")
    neo = items[0]
    p = neo.get("intercept_plan") or {}
    p["_neo_name"] = neo.get("name", "Target")
    return p

def _mag_au(val: Any, default: float = 1.0) -> float:
    """
    Convert a JSON value (float | [x,y] | [x,y,z] | np.ndarray | str) to a scalar AU.
    Safe for Pylance: all casts are internal; callers get a float.
    """
    try:
        # numpy array path (optional)
        try:
            import numpy as _np  # noqa: WPS433 (local import OK)
            if isinstance(val, _np.ndarray):
                arr = _np.asarray(val, dtype=float).ravel()
                if arr.size == 0:
                    return float(default)
                if arr.size == 1:
                    return float(arr[0])
                x = float(arr[0]); y = float(arr[1] if arr.size > 1 else 0.0); z = float(arr[2] if arr.size > 2 else 0.0)
                return float(math.hypot(x, math.hypot(y, z)))
        except Exception:
            pass

        # list/tuple path
        if isinstance(val, (list, tuple)):
            if not val:
                return float(default)
            x = float(val[0])
            y = float(val[1]) if len(val) > 1 else 0.0
            z = float(val[2]) if len(val) > 2 else 0.0
            return float(math.hypot(x, math.hypot(y, z)))

        # scalar path
        return float(val)
    except Exception:
        return float(default)

class Canvas(QGraphicsView):
    def __init__(self, plan: dict):
        super().__init__()
        # Use Qt6-safe RenderHint enum (lints clean)
        self.setRenderHint(QPainter.RenderHint.Antialiasing, True)

        pal = self.palette()
        pal.setColor(QPalette.ColorRole.Base, QColor(20,22,26))
        pal.setColor(QPalette.ColorRole.Text, QColor(220,220,220))
        self.setPalette(pal)

        self._scene = QGraphicsScene(self)
        self.setScene(self._scene)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)

        # static items
        self.earthOrbit = QGraphicsEllipseItem()
        self.neoOrbit   = QGraphicsEllipseItem()
        self.transfer   = QGraphicsEllipseItem()  # legacy Hohmann ellipse (hidden in Lambert mode)
        self.xferFull   = QGraphicsPathItem()     # Lambert polyline (full path)

        for it,color,w in [
            (self.earthOrbit, QColor(70,180,255), 2.2),
            (self.neoOrbit,   QColor(160,255,140), 2.2),
            (self.transfer,   QColor(245,200,120), 2.5),
        ]:
            it.setPen(QPen(color, w))
            it.setBrush(QBrush(Qt.BrushStyle.NoBrush))
            self._scene.addItem(it)

        self.xferFull.setPen(QPen(QColor(245,200,120), 2.5))
        self._scene.addItem(self.xferFull)

        # moving dots & labels
        self.earthDot = QGraphicsEllipseItem()
        self.neoDot   = QGraphicsEllipseItem()
        self.scDot    = QGraphicsEllipseItem()
        for it,color in [
            (self.earthDot, QColor(70,180,255)),
            (self.neoDot,   QColor(160,255,140)),
            (self.scDot,    QColor(245,140,90)),
        ]:
            it.setBrush(QBrush(color))
            it.setPen(QPen(Qt.PenStyle.NoPen))
            it.setRect(QRectF(-5,-5,10,10))
            self._scene.addItem(it)

        self.lblEarth = self._mk_label("Earth")
        self.lblNeo   = self._mk_label("Target")
        self.lblSc    = self._mk_label("SC", bold=True)
        self.note     = self._mk_label("", italic=True)

        # trail for SC along Lambert arc
        self.scTrack = QGraphicsPathItem()
        self.scTrack.setPen(QPen(QColor(245,160,100), 2.0, Qt.PenStyle.DashLine))
        self._scene.addItem(self.scTrack)

        self.setPlan(plan)

    def _mk_label(self, text, bold=False, italic=False):
        lbl = self._scene.addSimpleText(text)
        f = QFont("Inter", 10)
        f.setBold(bold); f.setItalic(italic)
        lbl.setFont(f)
        lbl.setBrush(QBrush(QColor(220,220,220)))
        return lbl

    # --------- plan & geometry ----------
    def setPlan(self, plan: dict):
        self.plan = plan

        # radii (this 2-D view uses circularized radii; if vectors are provided, use their magnitude)
        r1v = plan.get("r1_AU") or plan.get("r1_au") or 1.0
        r2v = plan.get("r2_AU") or plan.get("r2_au") or 1.0
        self.r1_phys = _mag_au(r1v, default=1.0)
        self.r2_phys = _mag_au(r2v, default=1.0)

        self.tof_days = float(plan.get("tof_days", 180.0))

        # mean motions
        self.n1 = mean_motion(self.r1_phys)
        self.n2 = mean_motion(self.r2_phys)

        # default legacy phase for target (π if outward, 0 if inward)
        tof_yrs = self.tof_days / 365.25
        rendez  = math.pi if self.r2_phys >= self.r1_phys else 0.0
        self.thetaE0 = 0.0
        self.thetaN0 = (rendez - self.n2 * tof_yrs) % (2*math.pi)

        # --- Orientation fix when NO Lambert: apply Δ(Ω+ω) ---
        no_lambert = not (plan.get("lambert_polyline_xy_au") or plan.get("lambert_polyline_xyz_au") or plan.get("lambert_poly_xyz_au"))
        if no_lambert:
            e = plan.get("elements_earth") or {}
            n = plan.get("elements_target") or plan.get("elements") or {}
            try:
                L_earth  = float(e.get("raan_deg",0.0)) + float(e.get("argp_deg",0.0))
                L_target = float(n.get("raan_deg",0.0)) + float(n.get("argp_deg",0.0))
                dL = (L_target - L_earth) * math.pi/180.0
                self.thetaN0 = (self.thetaN0 + dL) % (2*math.pi)
            except Exception:
                pass

        # geometry (drawn circles / ellipse)
        self.r1_draw = self.r1_phys
        self.r2_draw = self.r2_phys
        if abs(self.r2_draw - self.r1_draw) < 0.03:
            self.r2_draw += 0.12  # visual separation only
            self.note.setText("Note: target ring visually offset for clarity.")

        a = 0.5 * (self.r1_draw + self.r2_draw)
        b = math.sqrt(max(1e-6, self.r1_draw * self.r2_draw))
        e = math.sqrt(max(0.0, 1.0 - (b*b)/(a*a)))
        self.a_draw, self.b_draw, self.e_draw = a, b, e

        # Lambert support
        lam3 = plan.get("lambert_poly_xyz_au") or plan.get("lambert_polyline_xyz_au")
        lam2 = plan.get("lambert_polyline_xy_au")
        if lam3:
            pts = [(float(x), float(y)) for (x,y,_) in lam3]
        elif lam2:
            pts = [(float(x), float(y)) for (x,y) in lam2]
        else:
            pts = []
        self.lambert_pts = pts
        self.use_lambert = len(pts) >= 2

        if self.use_lambert:
            th_dep = math.atan2(pts[0][1], pts[0][0])
            th_arr = math.atan2(pts[-1][1], pts[-1][0])
            self.thetaE0 = th_dep
            self.thetaN0 = (th_arr - self.n2 * (self.tof_days / 365.25)) % (2*math.pi)

        # draw static geometry
        self._update_geometry()

        # time state
        self.t_days = 0.0
        self.playing = False
        self.speed = 1.0

        # precompute arc-length fractions for trail
        if self.use_lambert:
            self._lam_s = [0.0]
            total = 0.0
            for i in range(1, len(pts)):
                dx = pts[i][0] - pts[i-1][0]; dy = pts[i][1] - pts[i-1][1]
                total += math.hypot(dx, dy)
                self._lam_s.append(total)
            if total > 0:
                self._lam_s = [s/total for s in self._lam_s]
        else:
            self._lam_s = None

        self._tick()

    def _scale(self) -> float:
        vw=max(self.viewport().width(),800); vh=max(self.viewport().height(),500); margin=80.0
        R=max(self.r1_draw,self.r2_draw); pix=0.5*min(vw,vh)-margin
        return max(80.0,pix)/max(0.5,R)

    def _update_geometry(self):
        s = self._scale()
        r1 = self.r1_draw * s
        r2 = self.r2_draw * s
        a  = self.a_draw * s
        b  = self.b_draw * s
        c  = self.e_draw * a
        self.earthOrbit.setRect(QRectF(-r1,-r1,2*r1,2*r1))
        self.neoOrbit.setRect(QRectF(-r2,-r2,2*r2,2*r2))
        self.transfer.setRect(QRectF(-a + c, -b, 2*a, 2*b))
        self.transfer.setVisible(not self.use_lambert)

        # draw full Lambert path
        if self.use_lambert:
            path = QPainterPath(QPointF(self.lambert_pts[0][0]*s, self.lambert_pts[0][1]*s))
            for (x,y) in self.lambert_pts[1:]:
                path.lineTo(QPointF(x*s, y*s))
            self.xferFull.setPath(path)
            self.xferFull.setVisible(True)
        else:
            self.xferFull.setVisible(False)

    # --------- animation tick ----------
    def _tick(self):
        s = self._scale()
        t_yrs = self.t_days / 365.25

        thetaE = (self.thetaE0 + self.n1 * t_yrs) % (2*math.pi)
        thetaN = (self.thetaN0 + self.n2 * t_yrs) % (2*math.pi)

        x_e = self.r1_draw * s * math.cos(thetaE); y_e = self.r1_draw * s * math.sin(thetaE)
        x_n = self.r2_draw * s * math.cos(thetaN); y_n = self.r2_draw * s * math.sin(thetaN)
        self.earthDot.setRect(QRectF(x_e-6,y_e-6,12,12))
        self.neoDot.setRect(QRectF(x_n-6,y_n-6,12,12))

        # spacecraft position / trail
        if self.use_lambert and self._lam_s:
            f = max(0.0, min(1.0, self.t_days / max(self.tof_days, 1e-9)))
            import bisect
            j = min(len(self._lam_s)-1, max(1, bisect.bisect_left(self._lam_s, f)))
            x_sc, y_sc = self.lambert_pts[j]
            # trail up to j
            path = QPainterPath(QPointF(self.lambert_pts[0][0]*s, self.lambert_pts[0][1]*s))
            for k in range(1, j+1):
                path.lineTo(QPointF(self.lambert_pts[k][0]*s, self.lambert_pts[k][1]*s))
            self.scTrack.setPath(path)
            self.scDot.setRect(QRectF(x_sc*s-5,y_sc*s-5,10,10))
        else:
            # simple half-ellipse guide
            a = self.a_draw * s; b = self.b_draw * s; e = self.e_draw
            M = math.pi * (self.t_days / max(self.tof_days, 1e-9))
            x_sc = a * math.cos(M) - e * a; y_sc = b * math.sin(M)
            path = QPainterPath(QPointF(a - e*a, 0))
            for i in range(1, 181):
                Ei = (i/180.0) * M
                xi = a * math.cos(Ei) - e * a; yi = b * math.sin(Ei)
                path.lineTo(QPointF(xi, yi))
            self.scTrack.setPath(path)
            self.scDot.setRect(QRectF(x_sc-5,y_sc-5,10,10))

        # labels
        self.lblEarth.setPos(x_e+12, y_e-12)
        self.lblNeo.setPos(x_n+12,  y_n-12)
        self.lblSc.setPos(self.scDot.rect().x()+14, self.scDot.rect().y()-2)

    # --------- controls ---------
    def play(self, spd=1.0): self.playing=True; self.speed=spd
    def pause(self): self.playing=False
    def reset(self): self.t_days=0.0; self._tick()

    def step(self):
        if not self.playing: return
        self.t_days = min(self.tof_days, self.t_days + (self.tof_days/900.0)*self.speed)
        self._tick()

    def resizeEvent(self, e):
        super().resizeEvent(e); self._update_geometry(); self._tick()

class MainWindow(QMainWindow):
    def __init__(self, path: Path):
        super().__init__()
        self.setWindowTitle("NEO Updater — 2-D Viewer (Lambert)")
        self.plan = load_first_plan(path)

        # UI
        cw = QWidget(); self.setCentralWidget(cw)
        v = QVBoxLayout(cw); v.setContentsMargins(6,6,6,6); v.setSpacing(6)

        self.canvas = Canvas(self.plan); v.addWidget(self.canvas, 1)

        bar = QHBoxLayout(); v.addLayout(bar)
        self.btnPlay = QPushButton("Play"); self.btnPause = QPushButton("Pause"); self.btnReset = QPushButton("Reset")
        self.speed = QSlider(Qt.Orientation.Horizontal)
        self.speed.setMinimum(1); self.speed.setMaximum(40); self.speed.setValue(10)
        self.lbl = QLabel("Speed x1.0")
        for w in (self.btnPlay, self.btnPause, self.btnReset, QLabel("  Speed:"), self.speed, self.lbl):
            bar.addWidget(w)

        self.btnPlay.clicked.connect(lambda: self.canvas.play(self.speed.value()/10.0))
        self.btnPause.clicked.connect(self.canvas.pause)
        self.btnReset.clicked.connect(self.canvas.reset)
        self.speed.valueChanged.connect(lambda v: (self.lbl.setText(f"Speed x{v/10:.1f}"), self.canvas.play(v/10.0) if self.canvas.playing else None))

        # timer
        self.timer = QTimer(self); self.timer.timeout.connect(self.canvas.step); self.timer.start(16)

def main():
    p = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_INTERCEPT
    app = QApplication(sys.argv)
    w = MainWindow(p)
    w.resize(1100, 760)
    w.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
