#!/usr/bin/env python3
"""
NEO Updater — Qt Viewer (Lambert-aware)
- Big canvas (dark theme), table on bottom
- Distinct Earth/Target orbits + gold transfer
- Kepler-timed motion for Earth/Target and transfer (Hohmann half-ellipse) in legacy mode
- **Lambert mode**: draw solver polyline, animate SC along it, mark intercept at arc end
- Intercept vs closest approach detection; trail; zoom/pan
- Visual orbit offset when r1≈r2 (clarity only; HUD note shows)

Requires:  pip install PySide6
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
from typing import List, Tuple
from PySide6.QtCore import Qt, QAbstractTableModel, QModelIndex, QTimer, QRectF, QPointF
from PySide6.QtGui import (QPalette, QColor, QBrush, QPen, QFont, QAction,
                           QPainter, QPainterPath)
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QSplitter, QTableView, QFileDialog,
    QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QSlider, QToolBar, QMessageBox,
    QGraphicsView, QGraphicsScene, QGraphicsEllipseItem, QGraphicsSimpleTextItem,
    QGraphicsPathItem, QGraphicsLineItem, QStyleFactory
)
from statistics import median as _median

# --- Qt enum compatibility shims (Pylance-friendly) ---
from typing import Any
from PySide6.QtWidgets import QGraphicsView, QAbstractItemView  # add QAbstractItemView import

def _qt_enum(ns: str, member: str, fallback: Any) -> Any:
    nsobj = getattr(Qt, ns, None)
    return getattr(nsobj, member, fallback) if nsobj is not None else fallback

QtHorizontal       = _qt_enum("Orientation",    "Horizontal",       Qt.Horizontal)
QtVertical         = _qt_enum("Orientation",    "Vertical",         Qt.Vertical)
QtDashLine         = _qt_enum("PenStyle",       "DashLine",         Qt.DashLine)
QtNoPen            = _qt_enum("PenStyle",       "NoPen",            Qt.NoPen)
QtNoBrush          = _qt_enum("BrushStyle",     "NoBrush",          Qt.NoBrush)
QtAlignCenter      = _qt_enum("AlignmentFlag",  "AlignCenter",      Qt.AlignCenter)
QtKeepAspectRatio  = _qt_enum("AspectRatioMode","KeepAspectRatio",  Qt.KeepAspectRatio)

QGVScrollHandDrag  = getattr(getattr(QGraphicsView, "DragMode", QGraphicsView), 
                             "ScrollHandDrag", QGraphicsView.ScrollHandDrag)

AU_KM = 149_597_870.7
AU_M  = 149_597_870_700.0

DEFAULT_INTERCEPT = Path("data/hazardous_neos/latest_intercept.json")

# ---------------- Table ----------------
COLUMNS = [
    ("name","Name",220),("H_mag","H (mag)",80),
    ("dia_min_m","Dia min (m)",110),("dia_max_m","Dia max (m)",110),
    ("miss_km","Miss dist (km)",130),
    ("departure_utc","Departure UTC",170),("arrival_utc","Arrival UTC",170),
    ("dv_total","Δv total (m/s)",140),("dv_LEO","Δv LEO (m/s)",130),
    ("dv_arrive","Δv arrive (m/s)",140),("rolled_forward","Rolled?",80),
]
def _safe_float(x):
    try: return float(x)
    except Exception: return None
def _extract_miss_km(neo: dict):
    ca = neo.get("close_approach")
    if isinstance(ca, dict) and "miss_distance_km" in ca: return _safe_float(ca["miss_distance_km"])
    cad = neo.get("close_approach_data")
    if isinstance(cad, list) and cad:
        e = next((e for e in cad if e.get("orbiting_body")=="Earth"), cad[0])
        md = e.get("miss_distance") or {};  return _safe_float(md.get("kilometers") or md.get("km"))
    return None
def _round_sig(x, sig=3):
    if x in (None,0): return x
    try: return round(x, sig - int(math.floor(math.log10(abs(x)))) - 1)
    except Exception: return x

class NEOTableModel(QAbstractTableModel):
    def __init__(self, rows: list[dict]): super().__init__(); self.rows=rows
    def rowCount(self, parent=QModelIndex()): return len(self.rows)
    def columnCount(self, parent=QModelIndex()): return len(COLUMNS)
    def headerData(self, s, orient, role):
        if role!=Qt.DisplayRole: return None
        return COLUMNS[s][1] if orient==Qt.Horizontal else s+1
    def data(self, idx: QModelIndex, role=Qt.DisplayRole):
        if not idx.isValid(): return None
        row=self.rows[idx.row()]; key=COLUMNS[idx.column()][0]
        if role in (Qt.DisplayRole,Qt.EditRole): return row.get(key,"")
        if role==Qt.TextAlignmentRole: return int(Qt.AlignCenter)
        return None
    def sort(self, col, order):
        key=COLUMNS[col][0]
        def kf(r):
            v=r.get(key)
            try: return (False,float(v))
            except Exception: return (True,str(v))
        self.layoutAboutToBeChanged.emit()
        self.rows.sort(key=kf, reverse=(order==Qt.DescendingOrder))
        self.layoutChanged.emit()

def load_intercepts(path: Path) -> tuple[list[dict], list[dict]]:
    blob=json.load(open(path,"r",encoding="utf-8"))
    items=blob.get("potentially_hazardous_neos") or blob.get("neos") or []
    rows=[]
    for neo in items:
        p = neo.get("intercept_plan",{}) or {}
        rows.append({
            "name": neo.get("name"),
            "H_mag": neo.get("absolute_magnitude_h"),
            "dia_min_m": _round_sig(neo.get("estimated_diameter_m_min")),
            "dia_max_m": _round_sig(neo.get("estimated_diameter_m_max")),
            "miss_km": _round_sig(_extract_miss_km(neo)),
            "departure_utc": p.get("departure_utc"),
            "arrival_utc": p.get("arrival_utc"),
            "dv_LEO": _round_sig(p.get("dv_from_LEO_m_s")),
            "dv_arrive": _round_sig(p.get("dv_arrive_heliocentric_m_s")),
            "dv_total": _round_sig(p.get("dv_total_m_s")),
            "rolled_forward": "Yes" if p.get("rolled_forward") else "No",
            "_plan": p
        })
    return rows, items

# ------------- Orbital math (AU, years) -------------
def period_years(a_AU: float) -> float: return a_AU**1.5
def mean_motion(a_AU: float) -> float: return (2*math.pi)/period_years(a_AU)
def kepler_E_from_M(M: float, e: float, tol=1e-10, maxit=32) -> float:
    M=(M+2*math.pi)%(2*math.pi); E=M if e<0.8 else math.pi
    for _ in range(maxit):
        f=E-e*math.sin(E)-M; fp=1-e*math.cos(E)
        d=-f/fp; E+=d
        if abs(d)<tol: break
    return E

# ------------- Orbit canvas -------------
class OrbitView(QGraphicsView):
    """GMAT/KSP-like canvas with proper intercept timing.
    Legacy: Hohmann half-ellipse with phasing on circular/coplanar orbits.
    Lambert mode: draw polyline from planner, animate SC on it, intercept at arc end.
    """
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setRenderHints(self.renderHints() | QPainter.Antialiasing | QPainter.SmoothPixmapTransform)
        self.setBackgroundBrush(QBrush(QColor("#0e0f13")))
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self._scene = QGraphicsScene(self)
        self.setScene(self._scene)

        # ---- PHYSICS (AU) ----
        self.r1_phys = 1.0
        self.r2_phys = 1.0
        self.tof_days = 180.0

        # transfer ellipse (PHYS)
        self.a_phys = 1.0
        self.b_phys = 1.0
        self.e_phys = 0.0

        # mean motions (PHYS)
        self.n1 = mean_motion(self.r1_phys)
        self.n2 = mean_motion(self.r2_phys)

        # initial phases
        self.thetaE0 = 0.0   # Earth anchor (Lambert dep heading)
        self.thetaN0 = 0.0   # Target start so it arrives at Lambert endpoint

        # ---- DRAWING (AU, visual only) ----
        self.r1_draw = 1.0
        self.r2_draw = 1.0
        self.a_draw = 1.0
        self.b_draw = 1.0
        self.e_draw = 0.0
        self.visual_note = ""

        # Lambert
        self.use_lambert: bool = False
        self.lambert_pts: List[Tuple[float,float]] = []  # [(x_AU,y_AU), ...]

        # ---- TIME / STATE ----
        self.t_days = 0.0
        self.speed = 1.0
        self.playing = False
        self.stop_at_intercept = True

        # intercept / closest-approach markers
        self.intercept_found = False
        self.closest_px = float("inf")
        self.closest_pos = (0.0, 0.0)
        self.closest_t = 0.0

        # scene items
        self._new_items(); self._build_static()

        # tick timer
        self.timer = QTimer(self); self.timer.timeout.connect(self._tick); self.timer.start(int(1000/60))

    # ---------- scene items ----------
    def _new_items(self):
        self.gridRings = []
        self.axes = []
        self.legendItems = []
        self.legendSwatches = []

        self.sun = QGraphicsEllipseItem()
        self.earthOrbit = QGraphicsEllipseItem()
        self.neoOrbit = QGraphicsEllipseItem()
        self.transfer = QGraphicsEllipseItem()        # legacy Hohmann
        self.xferFull = QGraphicsPathItem()   # full Lambert polyline (background)
        self.scTrack = QGraphicsPathItem()
        self.earthDot = QGraphicsEllipseItem()
        self.neoDot = QGraphicsEllipseItem()
        self.scDot = QGraphicsEllipseItem()
        self.lblEarth = QGraphicsSimpleTextItem("Earth")
        self.lblNeo = QGraphicsSimpleTextItem("Target")
        self.lblSc = QGraphicsSimpleTextItem("SC")
        self.note = QGraphicsSimpleTextItem("")

        self.interceptRing = QGraphicsEllipseItem()
        self.interceptText = QGraphicsSimpleTextItem("")
        self.closestRing = QGraphicsEllipseItem()
        self.closestText = QGraphicsSimpleTextItem("")

    def _build_static(self):
        self._scene.clear(); self._new_items()
        gridPen=QPen(QColor("#2a2f3a"),1.0,Qt.DashLine); gridPen.setCosmetic(True)
        axisPen=QPen(QColor("#3b4252"),1.25); axisPen.setCosmetic(True)
        earthPen=QPen(QColor("#3fbdf1"),2.2); earthPen.setCosmetic(True)
        neoPen  =QPen(QColor("#a4e86e"),2.2); neoPen.setCosmetic(True)
        xferPen =QPen(QColor("#f6c177"),2.2); xferPen.setCosmetic(True)
        trailPen=QPen(QColor("#f29e74")); trailPen.setWidthF(2.0); trailPen.setCosmetic(True)
        dotPen=QPen(Qt.NoPen)
        earthBrush=QBrush(QColor("#3fbdf1")); neoBrush=QBrush(QColor("#a4e86e")); scBrush=QBrush(QColor("#f29e74"))
        sunBrush=QBrush(QColor("#ffd166"))
        interceptPen=QPen(QColor("#ffd480")); interceptPen.setWidthF(2.0); interceptPen.setCosmetic(True)
        closestPen=QPen(QColor("#9aa3b2")); closestPen.setStyle(Qt.DashLine); closestPen.setCosmetic(True)
        xferFullPen = QPen(QColor("#f6c177"), 1.2)  # thin, faint
        xferFullPen.setCosmetic(True); xferFullPen.setColor(QColor(246, 193, 119, 160))

        s=self._scale_px_per_AU(); Rmax=max(self.r1_draw,self.r2_draw)
        rings=max(4,int(math.ceil(Rmax/0.5))+1)
        for i in range(1,rings+1):
            r=i*0.5*s; el=QGraphicsEllipseItem(-r,-r,2*r,2*r)
            el.setPen(gridPen); self._scene.addItem(el); self.gridRings.append(el)
        L=(Rmax+0.75)*s; xAx=QGraphicsLineItem(-L,0,L,0); yAx=QGraphicsLineItem(0,-L,0,L)
        xAx.setPen(axisPen); yAx.setPen(axisPen); self._scene.addItem(xAx); self._scene.addItem(yAx); self.axes=[xAx,yAx]

        self.sun.setRect(QRectF(-7,-7,14,14)); self.sun.setBrush(sunBrush); self.sun.setPen(QPen(Qt.NoPen)); self._scene.addItem(self.sun)

        self.earthOrbit.setPen(earthPen); self.neoOrbit.setPen(neoPen)
        self.transfer.setPen(xferPen); self.transfer.setBrush(Qt.NoBrush)
        self._scene.addItem(self.earthOrbit); self._scene.addItem(self.neoOrbit); self._scene.addItem(self.transfer)

        self.scTrack.setPen(trailPen); self.scTrack.setBrush(Qt.NoBrush); self._scene.addItem(self.scTrack)

        for dot,br in ((self.earthDot,earthBrush),(self.neoDot,neoBrush),(self.scDot,scBrush)):
            dot.setPen(dotPen); dot.setBrush(br); self._scene.addItem(dot)

        for lbl,col in ((self.lblEarth,"#9cc9ff"),(self.lblNeo,"#c3f28e"),(self.lblSc,"#ffc9a3")):
            f=QFont("Helvetica",10); lbl.setFont(f); lbl.setBrush(QBrush(QColor(col))); self._scene.addItem(lbl)

        self._legend()

        self.note.setBrush(QBrush(QColor("#e9d5ac"))); self.note.setFont(QFont("Helvetica",9)); self._scene.addItem(self.note)

        self.interceptRing.setRect(QRectF(-6,-6,12,12)); self.interceptRing.setPen(interceptPen); self.interceptRing.setBrush(Qt.NoBrush)
        self.closestRing.setRect(QRectF(-5,-5,10,10)); self.closestRing.setPen(closestPen); self.closestRing.setBrush(Qt.NoBrush)
        self.interceptText.setBrush(QBrush(QColor("#ffd480"))); self.interceptText.setFont(QFont("Helvetica",9))
        self.closestText.setBrush(QBrush(QColor("#b8c0cc"))); self.closestText.setFont(QFont("Helvetica",9))
        self._scene.addItem(self.interceptRing); self._scene.addItem(self.interceptText)
        self._scene.addItem(self.closestRing); self._scene.addItem(self.closestText)
        self._scene.addItem(self.earthOrbit); self._scene.addItem(self.neoOrbit); self._scene.addItem(self.transfer)
        self.xferFull.setPen(xferFullPen); self.xferFull.setBrush(Qt.NoBrush); self._scene.addItem(self.xferFull)

        self.interceptRing.hide(); self.interceptText.hide()
        self.closestRing.hide(); self.closestText.hide()

        self._update_geometry()
        self.fitInView(self._scene.itemsBoundingRect().adjusted(-40,-40,40,40), Qt.KeepAspectRatio)

    def _legend(self):
        entries=[("Earth orbit","#3fbdf1"),("Target orbit","#a4e86e"),("Transfer (SC)","#f6c177")]
        self.legendSwatches.clear(); self.legendItems.clear()
        for txt,col in entries:
            sw=QGraphicsEllipseItem(0,0,10,10); sw.setBrush(QBrush(QColor(col))); sw.setPen(QPen(Qt.NoPen))
            t=QGraphicsSimpleTextItem(txt); t.setBrush(QBrush(QColor("#cdd6f4"))); t.setFont(QFont("Helvetica",9))
            self._scene.addItem(sw); self._scene.addItem(t)
            self.legendSwatches.append(sw); self.legendItems.append(t)

    def _normalize_lambert_units(self, pts_xy):
        """Return polyline points in AU. Accepts AU, km, or m; auto-detects scale."""
        if not pts_xy:
            return pts_xy
        radii = []
        clean = []
        for p in pts_xy:
            try:
                x, y = float(p[0]), float(p[1])
                if math.isfinite(x) and math.isfinite(y):
                    clean.append((x, y)); radii.append(math.hypot(x, y))
            except Exception:
                pass
        if not radii or len(clean) < 2:
            return []
        med_r = _median(radii)
        r_ref_au = max(self.r1_phys, self.r2_phys, 1.0)
        ratio = med_r / r_ref_au if r_ref_au > 0 else med_r
        if 1e-3 <= ratio <= 1e3:   # already AU-ish
            return clean
        scale = AU_M if med_r > 1e9 else AU_KM
        return [(x/scale, y/scale) for (x, y) in clean]

    def _clean_and_orient_lambert_pts(self, pts_xy):
        """Drop non-finite points, ensure first point ~departure (closest to r1), resample uniformly by arc length."""
        pts = [(float(x), float(y)) for (x, y) in pts_xy if math.isfinite(x) and math.isfinite(y)]
        if len(pts) < 2:
            return []
        r1 = self.r1_phys if self.r1_phys > 0 else 1.0
        d0 = abs(math.hypot(*pts[0]) - r1)
        dN = abs(math.hypot(*pts[-1]) - r1)
        if dN < d0:
            pts.reverse()

        # remove zero-length duplicates
        dedup = [pts[0]]
        for p in pts[1:]:
            if math.hypot(p[0]-dedup[-1][0], p[1]-dedup[-1][1]) > 1e-12:
                dedup.append(p)
        if len(dedup) < 2:
            return dedup

        # cumulative arc length
        cum = [0.0]
        for i in range(1, len(dedup)):
            dx = dedup[i][0]-dedup[i-1][0]; dy = dedup[i][1]-dedup[i-1][1]
            cum.append(cum[-1] + math.hypot(dx, dy))
        total = cum[-1]
        if total <= 0:
            return dedup

        # resample to stable N
        N = max(200, min(800, len(dedup)))
        import bisect
        out = []
        for k in range(N):
            s = (k/(N-1)) * total
            j = bisect.bisect_left(cum, s)
            if j <= 0: out.append(dedup[0]); continue
            if j >= len(cum): out.append(dedup[-1]); continue
            t = (s - cum[j-1]) / (cum[j]-cum[j-1] + 1e-18)
            x = dedup[j-1][0] + t*(dedup[j][0]-dedup[j-1][0])
            y = dedup[j-1][1] + t*(dedup[j][1]-dedup[j-1][1])
            out.append((x, y))
        return out

    def _prep_lambert_path(self, pts):
        """Precompute 0..1 arclength for fast sampling & build the full background QPainterPath."""
        if not pts or len(pts) < 2:
            self._lam_pts = []
            self._lam_s = []
            return QPainterPath()
        self._lam_pts = pts
        # cumulative length
        s = [0.0]
        for i in range(1, len(pts)):
            dx = pts[i][0]-pts[i-1][0]; dy = pts[i][1]-pts[i-1][1]
            s.append(s[-1] + math.hypot(dx, dy))
        L = s[-1] if s[-1] > 0 else 1.0
        self._lam_s = [si / L for si in s]  # normalized 0..1

        # full background path
        path = QPainterPath(QPointF(pts[0][0]*self._scale_px_per_AU(), pts[0][1]*self._scale_px_per_AU()))
        sa = self._scale_px_per_AU()
        for (x, y) in pts[1:]:
            path.lineTo(QPointF(x*sa, y*sa))
        return path

    def _sample_lambert_by_f(self, f: float):
        """Return (x_au, y_au) at fraction f of arc length (0..1)."""
        if not getattr(self, "_lam_pts", None) or not getattr(self, "_lam_s", None):
            return 0.0, 0.0
        f = min(max(f, 0.0), 1.0)
        import bisect
        j = bisect.bisect_left(self._lam_s, f)
        if j <= 0: return self._lam_pts[0]
        if j >= len(self._lam_s): return self._lam_pts[-1]
        t = (f - self._lam_s[j-1]) / (self._lam_s[j]-self._lam_s[j-1] + 1e-18)
        x = self._lam_pts[j-1][0] + t*(self._lam_pts[j][0]-self._lam_pts[j-1][0])
        y = self._lam_pts[j-1][1] + t*(self._lam_pts[j][1]-self._lam_pts[j-1][1])
        return x, y

    # ---------- plan & geometry ----------
    def setPlan(self, plan: dict):
        # PHYSICS radii
        self.r1_phys = float(plan.get("r1_AU", 1.0))
        self.r2_phys = float(plan.get("r2_AU", plan.get("r2", 1.0)))
        self.tof_days = float(plan.get("tof_days", 0.0)) or (0.5 * period_years(0.5*(self.r1_phys+self.r2_phys)) * 365.25 * 2)  # fallback

        # PHYS transfer ellipse (legacy)
        self.a_phys = 0.5 * (self.r1_phys + self.r2_phys)
        self.b_phys = math.sqrt(max(1e-12, self.r1_phys * self.r2_phys))
        self.e_phys = math.sqrt(max(0.0, 1.0 - (self.b_phys*self.b_phys)/(self.a_phys*self.a_phys)))

        # PHYS mean motions
        self.n1 = mean_motion(self.r1_phys)
        self.n2 = mean_motion(self.r2_phys)

        # Legacy target phase (π for outward, 0 for inward)
        tof_yrs = self.tof_days / 365.25
        rendez_angle = math.pi if self.r2_phys >= self.r1_phys else 0.0
        self.thetaN0 = (rendez_angle - self.n2 * tof_yrs) % (2*math.pi)  # will be overridden in Lambert mode

        # DRAW radii (visual separation only if r1≈r2)
        self.r1_draw = self.r1_phys
        self.r2_draw = self.r2_phys
        self.visual_note = ""
        if abs(self.r2_phys - self.r1_phys) < 0.03:
            self.r2_draw = self.r2_phys + 0.12
            self.visual_note = "Note: target orbit visually offset by +0.12 AU because r1≈r2 (clarity only)."

        # DRAW transfer ellipse (same E parameterization, different scale)
        self.a_draw = 0.5 * (self.r1_draw + self.r2_draw)
        self.b_draw = math.sqrt(max(1e-12, self.r1_draw * self.r2_draw))
        self.e_draw = math.sqrt(max(0.0, 1.0 - (self.b_draw*self.b_draw)/(self.a_draw*self.a_draw)))

        # Lambert mode?
        lambert_raw = plan.get("lambert")
        lambert = lambert_raw if isinstance(lambert_raw, dict) else {}

        pts_raw = plan.get("lambert_polyline_xy_au") or []
        pts_norm = self._normalize_lambert_units(pts_raw)
        pts = self._clean_and_orient_lambert_pts(pts_norm)

        self.use_lambert = bool(lambert and pts)
        self.lambert_pts = pts if self.use_lambert else []

        if self.use_lambert:
            # TOF from lambert wins if present
            tof_val = lambert.get("tof_days")
            if tof_val is not None:
                try: self.tof_days = float(tof_val)
                except Exception: pass

            # Precompute arclength & draw full path (background)
            self.xferFull.setVisible(True)
            self.xferFull.setPath(self._prep_lambert_path(self.lambert_pts))

            # Anchor Earth to departure heading; phase target to arrive at end
            th_dep = math.atan2(self.lambert_pts[0][1], self.lambert_pts[0][0])
            self.thetaE0 = th_dep
            th_arr = math.atan2(self.lambert_pts[-1][1], self.lambert_pts[-1][0])
            self.thetaN0 = (th_arr - self.n2 * (self.tof_days / 365.25)) % (2 * math.pi)
        else:
            self.thetaE0 = 0.0
            self.xferFull.setVisible(False)

        # reset state/markers
        self.t_days = 0.0; self.playing = False
        self.intercept_found = False; self.closest_px = float("inf")
        self.closestRing.hide(); self.closestText.hide()
        self.interceptRing.hide(); self.interceptText.hide()

        self._build_static(); self._update_state()

    def _scale_px_per_AU(self) -> float:
        vw=max(self.viewport().width(),800); vh=max(self.viewport().height(),500); margin=80.0
        R=max(self.r1_draw,self.r2_draw); pix=0.5*min(vw,vh)-margin
        return max(80.0,pix)/max(0.5,R)

    def _update_geometry(self):
        s = self._scale_px_per_AU()
        r1 = self.r1_draw * s
        r2 = self.r2_draw * s
        a  = self.a_draw * s
        b  = self.b_draw * s
        c  = self.e_draw * a

        self.earthOrbit.setRect(QRectF(-r1,-r1,2*r1,2*r1))
        self.neoOrbit.setRect(QRectF(-r2,-r2,2*r2,2*r2))
        self.transfer.setRect(QRectF(-a + c, -b, 2*a, 2*b))
        # Hide legacy ellipse in Lambert mode
        self.transfer.setVisible(not self.use_lambert)

        # legend & note positions
        bb=self._scene.itemsBoundingRect(); lx,ly=bb.left()+14, bb.top()+14
        for i,(sw,txt) in enumerate(zip(self.legendSwatches,self.legendItems)):
            sw.setRect(QRectF(lx, ly+i*18, 10,10)); txt.setPos(lx+14, ly-2+i*18)
        self.note.setText(self.visual_note or ("Lambert mode: trajectory from planner" if self.use_lambert else ""))
        if self.visual_note or self.use_lambert: self.note.setPos(bb.left()+14, bb.bottom()-24)
        self.transfer.setVisible(not self.use_lambert)
        self.xferFull.setVisible(self.use_lambert)

    # ---------- animation ----------
    def _tick(self):
        if not self.playing: return
        step_days = self.tof_days / 900.0 * self.speed
        self.t_days = min(self.tof_days, self.t_days + step_days)
        self._update_state()

    def set_time_fraction(self, f: float):
        self.t_days = max(0.0, min(self.tof_days, f * self.tof_days))
        self._update_state()

    def play(self, spd: float):
        self.speed = max(0.1, min(4.0, spd)); self.playing = True
    def pause(self): self.playing = False
    def reset(self):
        self.playing=False; self.t_days=0.0
        self.intercept_found=False; self.closest_px=float("inf")
        self.interceptRing.hide(); self.interceptText.hide()
        self.closestRing.hide(); self.closestText.hide()
        self._update_state()

    def _sc_xy_lambert(self, s_pix_per_au: float) -> tuple[float, float]:
        """SC xy in pixels at current time using arc-length parameterization."""
        if self.tof_days <= 0:
            return 0.0, 0.0
        f = max(0.0, min(1.0, self.t_days / self.tof_days))
        x_au, y_au = self._sample_lambert_by_f(f)
        return x_au * s_pix_per_au, y_au * s_pix_per_au


    def _update_state(self):
        """Advance using PHYSICS timing; draw with DRAW geometry."""
        s = self._scale_px_per_AU()
        t_yrs = self.t_days / 365.25

        # Earth & target angles
        if self.use_lambert:
            # Earth anchored to Lambert departure heading
            thetaE = (self.thetaE0 + self.n1 * t_yrs) % (2*math.pi)
        else:
            thetaE = (self.n1 * t_yrs) % (2*math.pi)

        thetaN = (self.n2 * t_yrs + self.thetaN0) % (2*math.pi)

        # Legacy transfer anomaly for Hohmann
        M = math.pi * (self.t_days / max(self.tof_days, 1e-9))
        E = kepler_E_from_M(M, self.e_phys)

        # Compute draw coordinates
        a = self.a_draw * s; b = self.b_draw * s; e = self.e_draw

        # Spacecraft XY + trail
        if self.use_lambert:
            # position on Lambert arc (arc-length)
            x_sc, y_sc = self._sc_xy_lambert(s)

            # grow trail up to current f using precomputed s array
            path = QPainterPath()
            if self.lambert_pts:
                f = max(0.0, min(1.0, self.t_days / max(self.tof_days, 1e-9)))
                import bisect
                j = bisect.bisect_left(self._lam_s, f) if getattr(self, "_lam_s", None) else 0
                if j <= 0:
                    x0, y0 = self.lambert_pts[0]; path.moveTo(QPointF(x0*s, y0*s))
                else:
                    x0, y0 = self.lambert_pts[0]; path.moveTo(QPointF(x0*s, y0*s))
                    for k in range(1, j+1):
                        xk, yk = self.lambert_pts[k]
                        path.lineTo(QPointF(xk*s, yk*s))
            self.scTrack.setPath(path)

        else:
            # legacy half-ellipse position/trail
            x_sc = a * math.cos(E) - e * a
            y_sc = b * math.sin(E)
            path = QPainterPath(QPointF(a - e*a, 0))
            samples = max(60, int(180 * (M / math.pi)))
            for i in range(1, samples+1):
                Ei = (i / samples) * E
                xi = a * math.cos(Ei) - e * a
                yi = b * math.sin(Ei)
                path.lineTo(QPointF(xi, yi))
            self.scTrack.setPath(path)

        # Earth/Target dots (target advances to reach Lambert arrival angle at TOF)
        x_e  = self.r1_draw * s * math.cos(thetaE)
        y_e  = self.r1_draw * s * math.sin(thetaE)

        x_n  = self.r2_draw * s * math.cos(thetaN)
        y_n  = self.r2_draw * s * math.sin(thetaN)

        # dots
        self.earthDot.setRect(QRectF(x_e-6,y_e-6,12,12))
        self.neoDot.setRect(QRectF(x_n-6,y_n-6,12,12))
        self.scDot.setRect(QRectF(x_sc-5,y_sc-5,10,10))

        # labels
        self._place_label(self.lblEarth,"Earth",x_e,y_e,self.r1_draw*s)
        self._place_label(self.lblNeo,"Target",x_n,y_n,self.r2_draw*s)
        self._place_label(self.lblSc,"SC",x_sc,y_sc,max(self.r1_draw,self.r2_draw)*s,priority=True)

        # intercept / closest-approach
        dx=x_sc-x_n; dy=y_sc-y_n; dist_px=math.hypot(dx,dy)
        tol_px = max(10.0, 0.003 * s)  # 0.003 AU or 10 px
        if dist_px < self.closest_px:
            self.closest_px = dist_px; self.closest_pos = (x_sc, y_sc); self.closest_t = self.t_days

        if not self.intercept_found:
            if self.use_lambert:
                # In Lambert mode, consider the last polyline point as intercept target
                if self.t_days >= self.tof_days - 1e-9:
                    self.intercept_found = True
                    self._place_marker(self.interceptRing, self.interceptText, x_sc, y_sc,
                                       f"Intercept (Lambert): t={self.t_days:.1f} d")
                    if self.stop_at_intercept: self.pause()
            else:
                if dist_px <= tol_px:
                    self.intercept_found = True
                    self._place_marker(self.interceptRing, self.interceptText, x_sc, y_sc,
                                       f"Intercept: t={self.t_days:.1f} d, Δ≈{(dist_px/s):.4f} AU")
                    if self.stop_at_intercept: self.pause()

        if self.t_days >= self.tof_days and not self.intercept_found and math.isfinite(self.closest_px):
            x,y = self.closest_pos
            self._place_marker(self.closestRing, self.closestText, x, y,
                               f"Closest approach: t={self.closest_t:.1f} d, Δ≈{(self.closest_px/s):.4f} AU")

    def _place_marker(self, ring, textItem, x, y, label):
        ring.setRect(QRectF(x-9,y-9,18,18)); ring.show()
        textItem.setText(label); textItem.setPos(x+12, y-12); textItem.show()

    def _place_label(self, lbl, text, x, y, r_ref, priority=False):
        lbl.setText(text)
        ang=math.atan2(y,x); delta=0.22 if priority else 0.12
        ang_off=ang+(delta if x>=0 else -delta); R=r_ref+30
        lx=R*math.cos(ang_off); ly=R*math.sin(ang_off)
        lbl.setPos(lx - lbl.boundingRect().width()/2, ly - lbl.boundingRect().height()/2)

    # interactivity
    def wheelEvent(self, e):
        factor=1.15 if e.angleDelta().y()>0 else 1/1.15
        self.scale(factor,factor)
    def resizeEvent(self,e):
        super().resizeEvent(e); self._update_geometry()

# ------------- Main window (canvas top, table bottom) -------------
class MainWindow(QMainWindow):
    def __init__(self, path: Path):
        super().__init__()
        self.setWindowTitle("NEO Updater — GMAT/KSP Viewer (Lambert)")
        self.resize(1320, 860)

        self.orbit=OrbitView()
        self.hud=QLabel("Select a NEO to preview"); self.hud.setStyleSheet("color:#e5e9f0; font: 14px 'Menlo'; padding:6px;")
        self.playBtn=QPushButton("▶ Play"); self.pauseBtn=QPushButton("⏸ Pause"); self.resetBtn=QPushButton("⟲ Reset")
        self.speed=QSlider(Qt.Horizontal); self.speed.setRange(10,400); self.speed.setValue(100)
        self.time=QSlider(Qt.Horizontal); self.time.setRange(0,1000); self.time.setValue(0)
        self.timeLbl=QLabel("t = 0.0 / 0.0 days"); self.timeLbl.setStyleSheet("color:#cdd6f4")

        ctl=QWidget(); c=QHBoxLayout(ctl); c.setContentsMargins(8,4,8,4); c.setSpacing(12)
        c.addWidget(self.hud,1); c.addWidget(QLabel("Speed")); c.addWidget(self.speed)
        c.addWidget(QLabel("Time")); c.addWidget(self.time,1); c.addWidget(self.timeLbl)
        c.addWidget(self.playBtn); c.addWidget(self.pauseBtn); c.addWidget(self.resetBtn)

        top=QWidget(); tL=QVBoxLayout(top); tL.setContentsMargins(0,0,0,0); tL.addWidget(self.orbit,1); tL.addWidget(ctl,0)

        self.table=QTableView(); self.table.setSortingEnabled(True); self.table.verticalHeader().setVisible(False)
        self.table.setAlternatingRowColors(True); self.table.setSelectionBehavior(QTableView.SelectRows)
        self.table.setWordWrap(False); self.table.setMinimumHeight(220)

        split=QSplitter(Qt.Vertical); split.addWidget(top); split.addWidget(self.table)
        split.setStretchFactor(0,1); split.setStretchFactor(1,0); self.setCentralWidget(split)

        tb=QToolBar("Main"); self.addToolBar(tb)
        actOpen=QAction("Open JSON…",self); actReload=QAction("Reload",self); actQuit=QAction("Quit",self)
        tb.addAction(actOpen); tb.addAction(actReload); tb.addSeparator(); tb.addAction(actQuit)
        actOpen.triggered.connect(self.openFile); actReload.triggered.connect(self.reload); actQuit.triggered.connect(self.close)

        self.playBtn.clicked.connect(lambda: self.orbit.play(self.speed.value()/100.0))
        self.pauseBtn.clicked.connect(self.orbit.pause)
        self.resetBtn.clicked.connect(self._reset_time)
        self.speed.valueChanged.connect(lambda _: None)
        self.time.valueChanged.connect(self._scrub_time)

        self.path=path; self.rows=[]; self._apply_dark_palette(); self._load()
        self.uiTimer=QTimer(self); self.uiTimer.timeout.connect(self._ui_tick); self.uiTimer.start(100)

    def _apply_dark_palette(self):
        QApplication.setStyle(QStyleFactory.create("Fusion"))
        p=QPalette(); back=QColor(18,20,26); panel=QColor(26,28,36); text=QColor(225,229,236); mid=QColor(44,48,58); acc=QColor(70,130,180)
        p.setColor(QPalette.Window,back); p.setColor(QPalette.Base,panel); p.setColor(QPalette.AlternateBase,back)
        p.setColor(QPalette.Text,text); p.setColor(QPalette.ButtonText,text); p.setColor(QPalette.Highlight,acc)
        p.setColor(QPalette.HighlightedText,QColor(255,255,255)); p.setColor(QPalette.Button,mid); p.setColor(QPalette.WindowText,text)
        QApplication.setPalette(p)

    # file ops
    def openFile(self):
        p,_=QFileDialog.getOpenFileName(self,"Open intercept JSON",str(self.path),"JSON files (*.json)")
        if p: self.path=Path(p); self._load()
    def reload(self): self._load()

    # time controls
    def _reset_time(self): self.orbit.reset(); self.time.setValue(0)
    def _scrub_time(self,v:int): self.orbit.set_time_fraction(v/1000.0)
    def _ui_tick(self):
        tof=self.orbit.tof_days; t=self.orbit.t_days
        self.timeLbl.setText(f"t = {t:.1f} / {tof:.1f} days")
        if self.orbit.playing and tof>0:
            frac=t/tof; self.time.blockSignals(True); self.time.setValue(int(frac*1000)); self.time.blockSignals(False)

    # model & selection
    def _load(self):
        try:
            self.rows,_=load_intercepts(self.path)
            self.model=NEOTableModel(self.rows); self.table.setModel(self.model)
            for i,(_,_,w) in enumerate(COLUMNS): self.table.setColumnWidth(i,w)
            self.table.selectionModel().selectionChanged.connect(lambda *_: self._onSelect())
            if self.rows: self.table.selectRow(0); self._onSelect()
            else: self.hud.setText("No NEOs found in file.")
        except Exception as e:
            QMessageBox.critical(self,"Load error",f"Failed to load:\n{e}")

    def _onSelect(self):
        if not self.table.selectionModel(): return
        idxs=self.table.selectionModel().selectedRows()
        if not idxs: return
        r=idxs[0].row(); row=self.rows[r]; plan=row["_plan"]
        # HUD dv/TOF from plan (Lambert TOF wins if present)
        lam=plan.get("lambert") or {}
        tof = lam.get("tof_days") or plan.get("tof_days")
        dv_total = plan.get("dv_total_m_s")
        self.hud.setText(f"{row['name']}  •  Δv={_round_sig(dv_total)} m/s  •  Depart {plan.get('departure_utc')} → Arrive {plan.get('arrival_utc')}"
                         + ("  •  Lambert" if lam else ""))
        self.orbit.setPlan(plan); self.time.setValue(0)

# ------------- Entry -------------
def main():
    path=DEFAULT_INTERCEPT
    if len(sys.argv)>1: path=Path(sys.argv[1])
    app=QApplication(sys.argv); win=MainWindow(path); win.show(); sys.exit(app.exec())
if __name__=="__main__": main()
