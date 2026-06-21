"""Interactive 3D solar-system mission viewer."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pyqtgraph.opengl as gl
from PySide6 import QtCore, QtGui, QtWidgets

from app.solar_system import (
    PLANETS,
    Mission,
    asteroid_orbit_curve_au,
    asteroid_position_au,
    format_jd,
    load_missions,
    orbit_curve_au,
    planet_elements,
    planet_position_au,
)

DEFAULT_INTERCEPT = Path("data/latest_intercepts.json")


class SolarSystemView(gl.GLViewWidget):
    """Map-style OpenGL view using astronomical units for all positions."""

    def __init__(self) -> None:
        super().__init__()
        self.setBackgroundColor((4, 8, 15))
        self.opts["fov"] = 58
        self.mission: Mission | None = None
        self.planet_orbits: dict[str, gl.GLLinePlotItem] = {}
        self.planet_dots: dict[str, gl.GLScatterPlotItem] = {}
        self._build_scene()
        self.show_full_system()

    def _build_scene(self) -> None:
        grid = gl.GLGridItem()
        grid.setSize(70, 70)
        grid.setSpacing(5, 5)
        grid.setColor((48, 71, 96, 70))
        grid.setDepthValue(100)
        self.addItem(grid)

        rng = np.random.default_rng(731)
        direction = rng.normal(size=(1100, 3))
        direction /= np.linalg.norm(direction, axis=1, keepdims=True)
        radii = rng.uniform(52.0, 64.0, size=(1100, 1))
        star_positions = (direction * radii).astype(np.float32)
        star_alpha = rng.uniform(0.15, 0.65, size=(1100, 1))
        star_colors = np.hstack([np.full((1100, 3), 0.82), star_alpha]).astype(np.float32)
        self.stars = gl.GLScatterPlotItem(
            pos=star_positions, color=star_colors, size=1.2, pxMode=True
        )
        self.stars.setDepthValue(1000)
        self.addItem(self.stars)

        self.sun_halo = gl.GLScatterPlotItem(
            pos=np.zeros((1, 3), dtype=np.float32),
            color=(1.0, 0.64, 0.10, 0.18),
            size=36,
            pxMode=True,
        )
        self.sun = gl.GLScatterPlotItem(
            pos=np.zeros((1, 3), dtype=np.float32),
            color=(1.0, 0.86, 0.35, 1.0),
            size=15,
            pxMode=True,
        )
        self.addItem(self.sun_halo)
        self.addItem(self.sun)

        for planet in PLANETS:
            orbit = gl.GLLinePlotItem(
                pos=np.empty((0, 3), dtype=np.float32),
                color=(*planet.color[:3], 0.38),
                width=1.15,
                antialias=True,
                mode="line_strip",
            )
            dot = gl.GLScatterPlotItem(
                pos=np.zeros((1, 3), dtype=np.float32),
                color=planet.color,
                size=planet.marker_size,
                pxMode=True,
            )
            self.planet_orbits[planet.name] = orbit
            self.planet_dots[planet.name] = dot
            self.addItem(orbit)
            self.addItem(dot)

        self.asteroid_orbit = gl.GLLinePlotItem(
            pos=np.empty((0, 3), dtype=np.float32),
            color=(0.40, 1.0, 0.60, 0.62),
            width=1.8,
            antialias=True,
            mode="line_strip",
        )
        self.transfer_path = gl.GLLinePlotItem(
            pos=np.empty((0, 3), dtype=np.float32),
            color=(1.0, 0.74, 0.27, 0.42),
            width=2.2,
            antialias=True,
            mode="line_strip",
        )
        self.transfer_trail = gl.GLLinePlotItem(
            pos=np.empty((0, 3), dtype=np.float32),
            color=(1.0, 0.68, 0.20, 1.0),
            width=3.4,
            antialias=True,
            mode="line_strip",
        )
        self.asteroid_dot = gl.GLScatterPlotItem(
            pos=np.zeros((1, 3), dtype=np.float32),
            color=(0.40, 1.0, 0.60, 1.0),
            size=9,
            pxMode=True,
        )
        self.spacecraft_dot = gl.GLScatterPlotItem(
            pos=np.zeros((1, 3), dtype=np.float32),
            color=(1.0, 0.52, 0.12, 1.0),
            size=10,
            pxMode=True,
        )
        self.departure_dot = gl.GLScatterPlotItem(
            pos=np.zeros((1, 3), dtype=np.float32),
            color=(0.25, 0.72, 1.0, 1.0),
            size=7,
            pxMode=True,
        )
        self.arrival_dot = gl.GLScatterPlotItem(
            pos=np.zeros((1, 3), dtype=np.float32),
            color=(0.40, 1.0, 0.60, 1.0),
            size=7,
            pxMode=True,
        )
        for item in (
            self.asteroid_orbit,
            self.transfer_path,
            self.transfer_trail,
            self.departure_dot,
            self.arrival_dot,
            self.asteroid_dot,
            self.spacecraft_dot,
        ):
            self.addItem(item)

    def set_mission(self, mission: Mission) -> None:
        self.mission = mission
        self.asteroid_orbit.setData(pos=asteroid_orbit_curve_au(mission.elements))
        self.transfer_path.setData(pos=mission.polyline_au)
        self.transfer_trail.setData(pos=mission.polyline_au[:1])
        self.departure_dot.setData(pos=mission.polyline_au[:1])
        self.arrival_dot.setData(pos=mission.polyline_au[-1:])
        for planet in PLANETS:
            elements = planet_elements(planet, mission.departure_jd)
            self.planet_orbits[planet.name].setData(pos=orbit_curve_au(elements))
        self.update_time(mission.departure_jd, 0.0)

    def update_time(self, jd_tdb: float, progress: float) -> np.ndarray:
        if self.mission is None:
            return np.zeros(3)
        for planet in PLANETS:
            self.planet_dots[planet.name].setData(
                pos=planet_position_au(planet, jd_tdb).reshape(1, 3).astype(np.float32)
            )
        asteroid_position = asteroid_position_au(self.mission.elements, jd_tdb)
        self.asteroid_dot.setData(pos=asteroid_position.reshape(1, 3).astype(np.float32))

        path = self.mission.polyline_au
        position = max(0.0, min(1.0, progress)) * (len(path) - 1)
        low = int(position)
        high = min(low + 1, len(path) - 1)
        fraction = position - low
        spacecraft = path[low] * (1.0 - fraction) + path[high] * fraction
        trail = np.vstack([path[: low + 1], spacecraft]).astype(np.float32)
        self.transfer_trail.setData(pos=trail)
        self.spacecraft_dot.setData(pos=spacecraft.reshape(1, 3))
        return spacecraft

    def show_full_system(self) -> None:
        self.setCameraPosition(
            pos=QtGui.QVector3D(0.0, 0.0, 0.0), distance=68.0, elevation=62, azimuth=-38
        )

    def show_inner_system(self) -> None:
        self.setCameraPosition(
            pos=QtGui.QVector3D(0.0, 0.0, 0.0), distance=5.2, elevation=68, azimuth=-35
        )

    def focus_mission(self) -> None:
        if self.mission is None:
            return
        path = self.mission.polyline_au
        center = np.mean(path, axis=0)
        extent = float(np.max(np.linalg.norm(path - center, axis=1)))
        self.setCameraPosition(
            pos=QtGui.QVector3D(*center.astype(float)),
            distance=max(2.8, extent * 3.2),
            elevation=54,
            azimuth=-35,
        )


class MissionViewer(QtWidgets.QMainWindow):
    TIMER_INTERVAL_MS = 33

    def __init__(self, path: Path):
        super().__init__()
        self.missions = load_missions(path)
        self.mission = self.missions[0]
        self.elapsed_days = 0.0
        self.playing = True
        self._scrubbing = False
        self._resume_after_scrub = False
        self.setWindowTitle("Asteroid Intercept Planner — 3D Mission Map")
        self.setMinimumSize(1180, 760)
        self.resize(1480, 900)
        self._build_ui()
        self._set_mission(0)
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._advance)
        self.timer.start(self.TIMER_INTERVAL_MS)

    def _build_ui(self) -> None:
        splitter = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        splitter.setChildrenCollapsible(False)
        self.setCentralWidget(splitter)

        self.view = SolarSystemView()
        splitter.addWidget(self.view)

        panel = QtWidgets.QFrame()
        panel.setObjectName("controlPanel")
        panel.setMinimumWidth(350)
        panel.setMaximumWidth(410)
        panel_layout = QtWidgets.QVBoxLayout(panel)
        panel_layout.setContentsMargins(24, 24, 24, 22)
        panel_layout.setSpacing(14)

        eyebrow = QtWidgets.QLabel("MISSION CONTROL")
        eyebrow.setObjectName("eyebrow")
        title = QtWidgets.QLabel("INTERCEPT MAP")
        title.setObjectName("panelTitle")
        panel_layout.addWidget(eyebrow)
        panel_layout.addWidget(title)

        panel_layout.addWidget(self._section_label("INTERCEPTION PLAN"))
        self.mission_combo = QtWidgets.QComboBox()
        self.mission_combo.setObjectName("missionCombo")
        for mission in self.missions:
            self.mission_combo.addItem(mission.label)
        self.mission_combo.currentIndexChanged.connect(self._set_mission)
        panel_layout.addWidget(self.mission_combo)

        self.date_label = QtWidgets.QLabel()
        self.date_label.setObjectName("dateLabel")
        self.phase_label = QtWidgets.QLabel("COASTING")
        self.phase_label.setObjectName("phaseLabel")
        self.elapsed_label = QtWidgets.QLabel()
        self.elapsed_label.setObjectName("elapsedLabel")
        panel_layout.addWidget(self.date_label)
        phase_row = QtWidgets.QHBoxLayout()
        phase_row.addWidget(self.phase_label)
        phase_row.addStretch(1)
        phase_row.addWidget(self.elapsed_label)
        panel_layout.addLayout(phase_row)

        self.timeline = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.timeline.setRange(0, 1000)
        self.timeline.sliderPressed.connect(self._begin_scrub)
        self.timeline.sliderReleased.connect(self._end_scrub)
        self.timeline.valueChanged.connect(self._scrub_to)
        panel_layout.addWidget(self.timeline)

        controls = QtWidgets.QHBoxLayout()
        self.play_button = QtWidgets.QPushButton("Pause")
        self.play_button.clicked.connect(self._toggle_playback)
        self.reset_button = QtWidgets.QPushButton("Restart")
        self.reset_button.setObjectName("secondaryButton")
        self.reset_button.clicked.connect(self._restart)
        controls.addWidget(self.play_button)
        controls.addWidget(self.reset_button)
        panel_layout.addLayout(controls)

        playback_row = QtWidgets.QHBoxLayout()
        playback_row.addWidget(QtWidgets.QLabel("Simulation speed"))
        self.speed_combo = QtWidgets.QComboBox()
        for label, value in (
            ("0.5 days/s", 0.5),
            ("1 day/s", 1.0),
            ("5 days/s", 5.0),
            ("10 days/s", 10.0),
            ("30 days/s", 30.0),
            ("100 days/s", 100.0),
        ):
            self.speed_combo.addItem(label, value)
        self.speed_combo.setCurrentIndex(3)
        playback_row.addWidget(self.speed_combo, 1)
        panel_layout.addLayout(playback_row)

        self.repeat_checkbox = QtWidgets.QCheckBox("Repeat mission on arrival")
        self.repeat_checkbox.setChecked(True)
        panel_layout.addWidget(self.repeat_checkbox)

        panel_layout.addWidget(self._section_label("MISSION TELEMETRY"))
        self.telemetry: dict[str, QtWidgets.QLabel] = {}
        telemetry_grid = QtWidgets.QGridLayout()
        telemetry_grid.setHorizontalSpacing(18)
        telemetry_grid.setVerticalSpacing(10)
        for row, (key, title_text) in enumerate(
            (
                ("target", "Target"),
                ("approach", "Close approach"),
                ("distance", "Miss distance"),
                ("tof", "Time of flight"),
                ("dv", "Total Δv"),
                ("c3", "C3 estimate"),
                ("position", "Spacecraft position"),
            )
        ):
            title_label = QtWidgets.QLabel(title_text)
            title_label.setObjectName("telemetryKey")
            value_label = QtWidgets.QLabel("—")
            value_label.setObjectName("telemetryValue")
            value_label.setAlignment(QtCore.Qt.AlignmentFlag.AlignRight)
            self.telemetry[key] = value_label
            telemetry_grid.addWidget(title_label, row, 0)
            telemetry_grid.addWidget(value_label, row, 1)
        panel_layout.addLayout(telemetry_grid)

        panel_layout.addWidget(self._section_label("CAMERA"))
        camera_row = QtWidgets.QHBoxLayout()
        for label, callback in (
            ("Mission", self.view.focus_mission),
            ("Inner", self.view.show_inner_system),
            ("Full", self.view.show_full_system),
        ):
            button = QtWidgets.QPushButton(label)
            button.setObjectName("cameraButton")
            button.clicked.connect(callback)
            camera_row.addWidget(button)
        panel_layout.addLayout(camera_row)

        help_text = QtWidgets.QLabel(
            "DRAG TO ORBIT  ·  MIDDLE-DRAG TO PAN  ·  SCROLL TO ZOOM\n"
            "SPACE PLAY/PAUSE  ·  R RESTART  ·  F FOCUS MISSION"
        )
        help_text.setObjectName("helpText")
        help_text.setWordWrap(True)
        panel_layout.addWidget(help_text)
        panel_layout.addStretch(1)

        legend = QtWidgets.QLabel(
            '<span style="color:#ffe16a">● Sun</span>  '
            '<span style="color:#b7ada0">● Mercury</span>  '
            '<span style="color:#f2b757">● Venus</span>  '
            '<span style="color:#43abff">● Earth</span>  '
            '<span style="color:#f55933">● Mars</span><br>'
            '<span style="color:#e0a86e">● Jupiter</span>  '
            '<span style="color:#e8cc85">● Saturn</span>  '
            '<span style="color:#75e0e8">● Uranus</span>  '
            '<span style="color:#476df5">● Neptune</span><br>'
            '<span style="color:#66ff99">● Target</span>  '
            '<span style="color:#ff851f">● Spacecraft</span><br>'
            "Positions and orbit geometry use AU; marker sizes are exaggerated."
        )
        legend.setObjectName("legend")
        legend.setWordWrap(True)
        panel_layout.addWidget(legend)

        splitter.addWidget(panel)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([1080, 380])
        self.setStyleSheet(STYLESHEET)

    def _section_label(self, text: str) -> QtWidgets.QLabel:
        label = QtWidgets.QLabel(text)
        label.setObjectName("sectionLabel")
        return label

    @QtCore.Slot(int)
    def _set_mission(self, index: int) -> None:
        if index < 0 or index >= len(self.missions):
            return
        self.mission = self.missions[index]
        self.elapsed_days = 0.0
        self.view.set_mission(self.mission)
        self._update_telemetry_static()
        self._render_time()
        QtCore.QTimer.singleShot(0, self.view.focus_mission)

    def _update_telemetry_static(self) -> None:
        mission = self.mission
        self.telemetry["target"].setText(mission.designation)
        self.telemetry["approach"].setText(mission.approach_text)
        self.telemetry["distance"].setText(f"{mission.distance_au:.6f} AU")
        self.telemetry["tof"].setText(f"{mission.tof_days:.1f} days")
        self.telemetry["dv"].setText(f"{mission.total_dv_kms:.3f} km/s")
        self.telemetry["c3"].setText(
            f"{mission.c3_km2_s2:.3f} km²/s²" if mission.c3_km2_s2 is not None else "—"
        )

    def _render_time(self) -> None:
        mission = self.mission
        progress = self.elapsed_days / max(mission.tof_days, 1e-9)
        jd = mission.departure_jd + self.elapsed_days
        spacecraft = self.view.update_time(jd, progress)
        self.date_label.setText(format_jd(jd))
        self.elapsed_label.setText(f"T+ {self.elapsed_days:06.2f} / {mission.tof_days:.1f} d")
        self.phase_label.setText("ARRIVAL" if progress >= 0.9999 else "COASTING")
        self.telemetry["position"].setText(
            f"{spacecraft[0]:+.3f}, {spacecraft[1]:+.3f}, {spacecraft[2]:+.3f} AU"
        )
        if not self._scrubbing:
            self.timeline.blockSignals(True)
            self.timeline.setValue(round(max(0.0, min(1.0, progress)) * 1000))
            self.timeline.blockSignals(False)

    @QtCore.Slot()
    def _advance(self) -> None:
        if not self.playing or self._scrubbing:
            return
        seconds = self.TIMER_INTERVAL_MS / 1000.0
        self.elapsed_days += float(self.speed_combo.currentData()) * seconds
        if self.elapsed_days >= self.mission.tof_days:
            if self.repeat_checkbox.isChecked():
                self.elapsed_days %= self.mission.tof_days
            else:
                self.elapsed_days = self.mission.tof_days
                self.playing = False
                self.play_button.setText("Play")
        self._render_time()

    @QtCore.Slot()
    def _toggle_playback(self) -> None:
        self.playing = not self.playing
        self.play_button.setText("Pause" if self.playing else "Play")

    @QtCore.Slot()
    def _restart(self) -> None:
        self.elapsed_days = 0.0
        self.playing = True
        self.play_button.setText("Pause")
        self._render_time()

    @QtCore.Slot()
    def _begin_scrub(self) -> None:
        self._resume_after_scrub = self.playing
        self._scrubbing = True
        self.playing = False

    @QtCore.Slot()
    def _end_scrub(self) -> None:
        self._scrubbing = False
        self.playing = self._resume_after_scrub
        self.play_button.setText("Pause" if self.playing else "Play")
        self._render_time()

    @QtCore.Slot(int)
    def _scrub_to(self, value: int) -> None:
        if not self._scrubbing:
            return
        self.elapsed_days = self.mission.tof_days * value / 1000.0
        self._render_time()

    def keyPressEvent(self, event: QtGui.QKeyEvent) -> None:
        if event.key() == QtCore.Qt.Key.Key_Space:
            self._toggle_playback()
        elif event.key() == QtCore.Qt.Key.Key_R:
            self._restart()
        elif event.key() == QtCore.Qt.Key.Key_F:
            self.view.focus_mission()
        else:
            super().keyPressEvent(event)


# Compatibility name retained for existing smoke tests and integrations.
MainWindow = MissionViewer


STYLESHEET = """
QMainWindow, QWidget { background: #070b12; color: #e7f0fa; font-family: Inter, "SF Pro Text", sans-serif; }
QFrame#controlPanel { background: #0d141e; border-left: 1px solid #263446; }
QLabel#eyebrow { color: #52d6ec; font-size: 10px; font-weight: 800; letter-spacing: 2px; }
QLabel#panelTitle { color: #f0f6fc; font-size: 26px; font-weight: 800; letter-spacing: 1px; }
QLabel#sectionLabel { color: #708399; border-bottom: 1px solid #263446; padding: 10px 0 6px 0; font-size: 10px; font-weight: 800; letter-spacing: 1.5px; }
QLabel#dateLabel { color: #f1f6fc; font-size: 20px; font-weight: 750; padding-top: 4px; }
QLabel#phaseLabel { color: #071015; background: #52d6ec; border-radius: 4px; padding: 3px 8px; font-size: 10px; font-weight: 900; }
QLabel#elapsedLabel, QLabel#telemetryKey { color: #778a9f; font-size: 11px; }
QLabel#telemetryValue { color: #e0e9f2; font-size: 11px; font-weight: 650; }
QLabel#helpText, QLabel#legend { color: #63758a; font-size: 9px; line-height: 1.4; }
QComboBox { background: #111c29; border: 1px solid #2b3d52; border-radius: 7px; padding: 9px 10px; color: #e7f0fa; }
QComboBox::drop-down { border: 0; width: 24px; }
QComboBox QAbstractItemView { background: #111c29; selection-background-color: #1b4e63; color: #e7f0fa; }
QPushButton { background: #52d6ec; color: #071015; border: 0; border-radius: 7px; padding: 10px 12px; font-weight: 800; }
QPushButton:hover { background: #74e2f3; }
QPushButton#secondaryButton, QPushButton#cameraButton { background: #111c29; color: #dbe7f2; border: 1px solid #2b3d52; }
QPushButton#secondaryButton:hover, QPushButton#cameraButton:hover { border-color: #52d6ec; }
QCheckBox { color: #9fb0c2; spacing: 8px; }
QCheckBox::indicator { width: 15px; height: 15px; }
QSlider::groove:horizontal { height: 4px; background: #263446; border-radius: 2px; }
QSlider::sub-page:horizontal { background: #52d6ec; border-radius: 2px; }
QSlider::handle:horizontal { background: #f3f8fc; width: 14px; margin: -5px 0; border-radius: 7px; }
QSplitter::handle { background: #263446; width: 1px; }
"""


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_INTERCEPT
    application = QtWidgets.QApplication(sys.argv)
    application.setApplicationName("Asteroid Intercept Planner")
    window = MissionViewer(path)
    window.show()
    window.raise_()
    window.activateWindow()
    QtCore.QTimer.singleShot(0, window.raise_)
    QtCore.QTimer.singleShot(0, window.activateWindow)
    return application.exec()


if __name__ == "__main__":
    sys.exit(main())
