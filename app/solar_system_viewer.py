# app/solar_system_viewer.py
"""
Solar System 3-D viewer (Qt / pyqtgraph / OpenGL).

All pure orbital math now lives in app/orbital_math.py so that scripts that
only need the math can import it without pulling in PySide6/pyqtgraph.
"""
from __future__ import annotations

from datetime import datetime
from typing import Dict, Mapping

import numpy as np
from PySide6 import QtCore, QtGui, QtWidgets  # type: ignore
import pyqtgraph as pg                          # type: ignore
import pyqtgraph.opengl as gl                   # type: ignore

# ---------------------------------------------------------------------------
# Re-export everything from orbital_math so existing callers keep working.
# ---------------------------------------------------------------------------
from app.orbital_math import (  # noqa: F401  (re-exports)
    GAUSSIAN_GRAVITATIONAL_CONSTANT,
    JD_UNIX_EPOCH,
    EARTH_RADIUS_KM,
    PLANET_SCALE_BASE_SIZE,
    PLANET_PIXEL_MIN_SIZE,
    PLANET_WORLD_BASE_SIZE,
    PLANET_WORLD_MIN_SIZE,
    SUN_WORLD_SCALE,
    PLANETARY_DATA,
    SUN_DATA,
    datetime_to_jd,
    kepler_E_from_M,
    planet_pos_au,
    elliptical_orbit_points,
)


class SolarSystemViewer(QtWidgets.QMainWindow):
    """Main window displaying planetary orbits with animated motion."""

    def __init__(self, day_step: float = 1.0, refresh_ms: int = 25) -> None:
        super().__init__()
        self.setWindowTitle("Solar System Viewer")
        self.view = gl.GLViewWidget()
        self.view.setMinimumSize(900, 700)
        self.view.setBackgroundColor(QtGui.QColor("#0e0f13"))
        self.setCentralWidget(self.view)

        axes = gl.GLAxisItem()
        axes.setSize(2.0, 2.0, 2.0)
        self.view.addItem(axes)

        # Reference 1-AU ring
        ring = elliptical_orbit_points(
            {"a": 1.0, "e": 0.0, "i": 0.0, "raan": 0.0, "argp": 0.0}, 256
        )
        ring_item = gl.GLLinePlotItem(
            pos=ring, width=1.2, color=(0.2, 0.6, 1.0, 0.35),
            antialias=True, mode="line_strip",
        )
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
        for name, data in PLANETARY_DATA.items():
            elems: Mapping[str, float] = data["elements"]
            color = data["color"]

            orbit = elliptical_orbit_points(
                {"a": elems["a"], "e": elems["e"], "i": elems["i"],
                 "raan": elems["raan"], "argp": elems["argp"]},
                num_points=512,
            )
            orbit_item = gl.GLLinePlotItem(
                pos=orbit, width=1.5, color=color, antialias=True, mode="line_strip",
            )
            self.view.addItem(orbit_item)
            self.orbit_items[name] = orbit_item

            radius_km = float(data.get("radius", EARTH_RADIUS_KM))
            size_px = PLANET_SCALE_BASE_SIZE * (radius_km / EARTH_RADIUS_KM) ** (1.0 / 3.0)
            if size_px < PLANET_PIXEL_MIN_SIZE:
                size_px = PLANET_PIXEL_MIN_SIZE

            pos = planet_pos_au(elems, self.sim_jd)
            scatter = gl.GLScatterPlotItem(
                pos=pos.reshape(1, 3), size=size_px, color=color, pxMode=True,
            )
            self.view.addItem(scatter)
            self.planet_items[name] = scatter

        sun_radius_km = float(SUN_DATA["radius"])
        sun_color = SUN_DATA["color"]
        sun_size_world = (
            PLANET_WORLD_BASE_SIZE * SUN_WORLD_SCALE
            * (sun_radius_km / EARTH_RADIUS_KM) ** (1.0 / 3.0)
        )
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
