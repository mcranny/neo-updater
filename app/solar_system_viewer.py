from __future__ import annotations

import math
from datetime import datetime, timezone
from typing import Mapping, Dict, Any, Tuple

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets  # type: ignore
import pyqtgraph as pg  # type: ignore
import pyqtgraph.opengl as gl  # type: ignore

# ----------------------------------------------------------------------------- 
# Orbital Constants
GAUSSIAN_GRAVITATIONAL_CONSTANT: float = 0.01720209895  # rad/day
JD_UNIX_EPOCH: float = 2440587.5

# Physical constants used for scaling scatter markers.
EARTH_RADIUS_KM: float = 6371.0  # km:contentReference[oaicite:0]{index=0}

# Pixel sizing constants for planets (static size on screen).
PLANET_SCALE_BASE_SIZE: float = 8.0   # Earth is ~8 px in diameter
PLANET_PIXEL_MIN_SIZE: float = 4.0    # Smallest planet marker size in pixels

# World‑unit sizing constants for the Sun (scales with zoom).
PLANET_WORLD_BASE_SIZE: float = 0.07  # Base size in AU (unused for planets here)
PLANET_WORLD_MIN_SIZE: float = 0.02   # Minimum world-unit size for markers
SUN_WORLD_SCALE: float = 0.8          # Multiplier to keep Sun within Mercury’s orbit

# Planetary data including orbital elements, colours and radii (km).
PLANETARY_DATA: Dict[str, Dict[str, Any]] = {
    "Mercury": {
        "elements": {
            "a": 0.3870983288631336,
            "e": 0.2056394563067703,
            "i": 7.003939465033488,
            "raan": 48.30843876040714,
            "argp": 29.17468769891704,
            "ma": math.degrees(1.442836675876093),  # mean anomaly in degrees
            "epoch": 2458100.340972222,
        },
        "color": (1.0, 0.0, 0.7333, 1.0),
        "radius": 2440.0,
    },
    "Venus": {
        "elements": {
            "a": 0.7233326902634808,
            "e": 0.00679866823224107,
            "i": 3.394487196632463,
            "raan": 76.62845518181165,
            "argp": 54.77727994493127,
            "ma": 113.0862278725096,
            "epoch": 2458100.340972222,
        },
        "color": (0.627, 0.439, 0.816, 1.0),
        "radius": 6052.0,
    },
    "Earth": {
        "elements": {
            "a": 1.00039544881957,
            "e": 0.01691835593284274,
            "i": 0.00280690387451405,
            "raan": 194.6496077844279,
            "argp": 266.7095883769374,
            "ma": 340.0687442843391,
            "epoch": 2458100.340972222,
        },
        "color": (0.125, 0.753, 1.0, 1.0),
        "radius": 6371.0,
    },
    "Mars": {
        "elements": {
            "a": 1.523723025279581,
            "e": 0.0934347086557555,
            "i": 1.848349492088388,
            "raan": 49.50757782581974,
            "argp": 286.6063237928008,
            "ma": 214.5554350827125,
            "epoch": 2458100.340972222,
        },
        "color": (1.0, 0.251, 0.063, 1.0),
        "radius": 3390.0,
    },
    "Jupiter": {
        "elements": {
            "a": 5.202042163311189,
            "e": 0.0489095120219927,
            "i": 1.303753431149363,
            "raan": 100.5116502623224,
            "argp": 273.7611453018681,
            "ma": 204.8067941190773,
            "epoch": 2458100.341666667,
        },
        "color": (0.96, 0.78, 0.22, 1.0),
        "radius": 69911.0,
    },
    "Saturn": {
        "elements": {
            "a": 9.574258042751591,
            "e": 0.05132315624461592,
            "i": 2.490161516580946,
            "raan": 113.5570492914928,
            "argp": 339.9959173823489,
            "ma": 175.6992517708274,
            "epoch": 2458100.342361111,
        },
        "color": (0.93, 0.84, 0.59, 1.0),
        "radius": 58232.0,
    },
    "Uranus": {
        "elements": {
            "a": 19.13080164192368,
            "e": 0.04936403218265095,
            "i": 0.7718998849235366,
            "raan": 73.98374777181797,
            "argp": 99.50085361825165,
            "ma": 216.5850480903659,
            "epoch": 2458100.343055556,
        },
        "color": (0.49, 0.97, 0.93, 1.0),
        "radius": 25362.0,
    },
    "Neptune": {
        "elements": {
            "a": 30.05366114914086,
            "e": 0.006670118099881116,
            "i": 1.768689065346207,
            "raan": 131.7490116077552,
            "argp": 268.9906654228161,
            "ma": 303.1122564055189,
            "epoch": 2458100.34375,
        },
        "color": (0.25, 0.45, 0.82, 1.0),
        "radius": 24622.0,
    },
}

# Sun data for static position but scalable size in world units.
SUN_DATA: Dict[str, Any] = {
    "radius": 695700.0,
    "color": (1.0, 0.85, 0.0, 1.0),
}

def datetime_to_jd(dt: datetime) -> float:
    """Convert a datetime (UTC) to Julian day number."""
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return JD_UNIX_EPOCH + dt.timestamp() / 86400.0

def kepler_E_from_M(M: float, e: float, tol: float = 1e-12, maxit: int = 50) -> float:
    """Solve Kepler's equation E − e·sin(E) = M for the eccentric anomaly E."""
    M = (M + math.pi) % (2.0 * math.pi) - math.pi
    E = M if e < 0.8 else math.pi
    for _ in range(maxit):
        f = E - e * math.sin(E) - M
        fp = 1.0 - e * math.cos(E)
        dE = -f / fp
        E += dE
        if abs(dE) < tol:
            break
    return E

def planet_pos_au(elements: Mapping[str, float], jd: float) -> np.ndarray:
    """Compute the planet's position in astronomical units at a given Julian date."""
    a = float(elements["a"]); e = float(elements["e"])
    i = math.radians(float(elements["i"]))
    raan = math.radians(float(elements["raan"]))
    argp = math.radians(float(elements["argp"]))
    ma0 = math.radians(float(elements["ma"]))
    epoch = float(elements["epoch"])
    n = GAUSSIAN_GRAVITATIONAL_CONSTANT / math.sqrt(a ** 3)
    dt = jd - epoch
    M = ma0 + n * dt
    E = kepler_E_from_M(M, e)
    cosE, sinE = math.cos(E), math.sin(E)
    x_p = a * (cosE - e)
    y_p = a * (math.sqrt(1.0 - e * e) * sinE)
    z_p = 0.0
    cO, sO = math.cos(raan), math.sin(raan)
    ci, si = math.cos(i), math.sin(i)
    cw, sw = math.cos(argp), math.sin(argp)
    R = np.array([
        [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
        [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
        [sw * si, cw * si, ci],
    ], dtype=float)
    r_pf = np.array([x_p, y_p, z_p])
    return (R @ r_pf).astype(np.float32)

def elliptical_orbit_points(elements: Mapping[str, float], num_points: int = 512) -> np.ndarray:
    """Generate an array of points approximating a planet's orbit."""
    a = float(elements["a"]); e = float(elements["e"])
    i = math.radians(float(elements["i"]))
    raan = math.radians(float(elements["raan"]))
    argp = math.radians(float(elements["argp"]))
    nu = np.linspace(0.0, 2.0 * math.pi, num_points, dtype=np.float64)
    p = a * (1.0 - e * e)
    r = p / (1.0 + e * np.cos(nu))
    xp = r * np.cos(nu); yp = r * np.sin(nu); zp = np.zeros_like(xp)
    rp = np.vstack((xp, yp, zp))
    cO, sO = math.cos(raan), math.sin(raan)
    ci, si = math.cos(i), math.sin(i)
    cw, sw = math.cos(argp), math.sin(argp)
    R = np.array([
        [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
        [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
        [sw * si, cw * si, ci],
    ], dtype=float)
    return (R @ rp).T.astype(np.float32)

class SolarSystemViewer(QtWidgets.QMainWindow):
    """Main window displaying the planetary orbits and animating their motion."""

    def __init__(self, day_step: float = 1.0, refresh_ms: int = 25) -> None:
        super().__init__()
        self.setWindowTitle("Solar System Viewer")
        self.view = gl.GLViewWidget()
        self.view.setMinimumSize(900, 700)
        self.view.setBackgroundColor(QtGui.QColor("#0e0f13"))
        self.setCentralWidget(self.view)
        axes = gl.GLAxisItem(); axes.setSize(2.0, 2.0, 2.0)
        self.view.addItem(axes)
        # Reference 1‑AU ring
        ring = elliptical_orbit_points({"a": 1.0, "e": 0.0, "i": 0.0,
                                        "raan": 0.0, "argp": 0.0}, 256)
        ring_item = gl.GLLinePlotItem(pos=ring, width=1.2,
                                      color=(0.2, 0.6, 1.0, 0.35),
                                      antialias=True, mode="line_strip")
        self.view.addItem(ring_item)
        self.view.setCameraPosition(distance=5.0, azimuth=45.0, elevation=22.0)
        self.sim_jd = datetime_to_jd(datetime.utcnow())
        self.day_step = float(day_step)
        self.orbit_items: Dict[str, gl.GLLinePlotItem] = {}
        self.planet_items: Dict[str, gl.GLScatterPlotItem] = {}
        self._build_scene()
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._tick)
        self.timer.setInterval(int(refresh_ms))
        self.timer.start()

    def _build_scene(self) -> None:
        """Create and add the orbital curves and scatter markers to the OpenGL view."""
        for name, data in PLANETARY_DATA.items():
            elems: Mapping[str, float] = data["elements"]
            color = data["color"]
            # Draw orbit in world coordinates
            orbit = elliptical_orbit_points({
                "a": elems["a"], "e": elems["e"], "i": elems["i"],
                "raan": elems["raan"], "argp": elems["argp"],
            }, num_points=512)
            orbit_item = gl.GLLinePlotItem(pos=orbit, width=1.5,
                                           color=color, antialias=True, mode="line_strip")
            self.view.addItem(orbit_item)
            self.orbit_items[name] = orbit_item
            # Planet markers use a fixed pixel size (pxMode=True)
            radius_km = float(data.get("radius", EARTH_RADIUS_KM))
            size_px = PLANET_SCALE_BASE_SIZE * (radius_km / EARTH_RADIUS_KM) ** (1.0/3.0)
            if size_px < PLANET_PIXEL_MIN_SIZE:
                size_px = PLANET_PIXEL_MIN_SIZE
            pos = planet_pos_au(elems, self.sim_jd)
            scatter = gl.GLScatterPlotItem(
                pos=pos.reshape(1, 3),
                size=size_px,
                color=color,
                pxMode=True,
            )
            self.view.addItem(scatter)
            self.planet_items[name] = scatter

        # Add the Sun with size in world units (pxMode=False)
        sun_radius_km = float(SUN_DATA["radius"])
        sun_color = SUN_DATA["color"]
        sun_size_world = PLANET_WORLD_BASE_SIZE * SUN_WORLD_SCALE \
                         * (sun_radius_km / EARTH_RADIUS_KM) ** (1.0/3.0)
        if sun_size_world < PLANET_WORLD_MIN_SIZE:
            sun_size_world = PLANET_WORLD_MIN_SIZE
        sun_scatter = gl.GLScatterPlotItem(
            pos=np.array([[0.0, 0.0, 0.0]]),
            size=sun_size_world,
            color=sun_color,
            pxMode=False,
        )
        self.view.addItem(sun_scatter)

    def _tick(self) -> None:
        """Advance the simulation and update planet positions."""
        self.sim_jd += self.day_step
        for name, data in PLANETARY_DATA.items():
            elems: Mapping[str, float] = data["elements"]
            pos = planet_pos_au(elems, self.sim_jd)
            scatter = self.planet_items.get(name)
            if scatter:
                scatter.setData(pos=pos.reshape(1, 3))

def main() -> int:
    app = QtWidgets.QApplication([])
    viewer = SolarSystemViewer(day_step=1.0, refresh_ms=25)
    viewer.show()
    return app.exec()

if __name__ == "__main__":
    raise SystemExit(main())
