"""Unified PySide6 launchpad for the dashboard, ingestion, and orbit viewer."""

from __future__ import annotations

import sys
import urllib.error
import urllib.request

from PySide6 import QtCore, QtGui, QtWidgets

from app.config import ROOT_DIR, get_settings
from app.database import dashboard_stats, initialize_database, latest_ingestion
from app.export import export_viewer_json


class Launchpad(QtWidgets.QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.settings = get_settings()
        initialize_database(self.settings.database_path)
        self.dashboard_process = QtCore.QProcess(self)
        self.update_process = QtCore.QProcess(self)
        self.viewer_process = QtCore.QProcess(self)
        self.setWindowTitle("Asteroid Intercept Planner")
        self.setMinimumSize(880, 570)
        self._build_ui()
        self._connect_processes()
        self.refresh_status()

    def _build_ui(self) -> None:
        root = QtWidgets.QWidget()
        root.setObjectName("root")
        self.setCentralWidget(root)
        layout = QtWidgets.QVBoxLayout(root)
        layout.setContentsMargins(42, 36, 42, 36)
        layout.setSpacing(24)

        eyebrow = QtWidgets.QLabel("JPL DATA  /  SQLITE  /  ORBITAL MECHANICS")
        eyebrow.setObjectName("eyebrow")
        title = QtWidgets.QLabel("Asteroid Intercept Planner")
        title.setObjectName("title")
        subtitle = QtWidgets.QLabel(
            "Update the local catalog, explore close approaches, and visualize "
            "computed Lambert transfers from one control center."
        )
        subtitle.setWordWrap(True)
        subtitle.setObjectName("subtitle")
        layout.addWidget(eyebrow)
        layout.addWidget(title)
        layout.addWidget(subtitle)

        cards = QtWidgets.QHBoxLayout()
        cards.setSpacing(18)
        dashboard_card, self.dashboard_button = self._card(
            "DATA FRONTEND",
            "Explore the database",
            "Browse the SQLite catalog, inspect transfer metrics, run read-only SQL, and export CSV.",
            "Launch dashboard",
        )
        viewer_card, self.viewer_button = self._card(
            "ORBIT VIEWER",
            "See the trajectory",
            "Animate Earth-to-asteroid transfer geometry and mission timing from the latest dataset.",
            "Launch viewer",
        )
        cards.addWidget(dashboard_card)
        cards.addWidget(viewer_card)
        layout.addLayout(cards, 1)

        status_bar = QtWidgets.QFrame()
        status_bar.setObjectName("statusBar")
        status_layout = QtWidgets.QHBoxLayout(status_bar)
        status_layout.setContentsMargins(18, 14, 18, 14)
        self.status_label = QtWidgets.QLabel()
        self.status_label.setObjectName("status")
        self.update_button = QtWidgets.QPushButton("Update from JPL")
        self.update_button.setObjectName("secondaryButton")
        status_layout.addWidget(self.status_label, 1)
        status_layout.addWidget(self.update_button)
        layout.addWidget(status_bar)

        self.dashboard_button.clicked.connect(self.start_dashboard)
        self.viewer_button.clicked.connect(self.start_viewer)
        self.update_button.clicked.connect(self.start_update)
        self.setStyleSheet(STYLESHEET)

    def _card(self, kicker: str, headline: str, body: str, button: str):
        frame = QtWidgets.QFrame()
        frame.setObjectName("card")
        layout = QtWidgets.QVBoxLayout(frame)
        layout.setContentsMargins(25, 25, 25, 25)
        layout.setSpacing(12)
        label = QtWidgets.QLabel(kicker)
        label.setObjectName("cardKicker")
        title = QtWidgets.QLabel(headline)
        title.setObjectName("cardTitle")
        description = QtWidgets.QLabel(body)
        description.setWordWrap(True)
        description.setObjectName("cardBody")
        action = QtWidgets.QPushButton(button)
        layout.addWidget(label)
        layout.addWidget(title)
        layout.addWidget(description)
        layout.addStretch(1)
        layout.addWidget(action)
        return frame, action

    def _connect_processes(self) -> None:
        self.dashboard_process.readyReadStandardError.connect(
            lambda: self._capture_process_error(self.dashboard_process)
        )
        self.update_process.finished.connect(self._update_finished)
        self.update_process.readyReadStandardError.connect(
            lambda: self._capture_process_error(self.update_process)
        )

    def _process_environment(self) -> QtCore.QProcessEnvironment:
        environment = QtCore.QProcessEnvironment.systemEnvironment()
        environment.insert("PYTHONPATH", str(ROOT_DIR))
        return environment

    def _start_module(
        self, process: QtCore.QProcess, module: str, arguments: list[str] | None = None
    ) -> None:
        process.setWorkingDirectory(str(ROOT_DIR))
        process.setProcessEnvironment(self._process_environment())
        process.start(sys.executable, ["-m", module, *(arguments or [])])

    @QtCore.Slot()
    def start_dashboard(self) -> None:
        if self.dashboard_process.state() == QtCore.QProcess.ProcessState.NotRunning:
            self._start_module(self.dashboard_process, "app.web")
        self.status_label.setText("Starting the local dashboard…")
        QtCore.QTimer.singleShot(700, self._open_dashboard_when_ready)

    def _open_dashboard_when_ready(self, attempts: int = 0) -> None:
        host = "127.0.0.1" if self.settings.flask_host == "0.0.0.0" else self.settings.flask_host
        url = f"http://{host}:{self.settings.flask_port}"
        try:
            with urllib.request.urlopen(f"{url}/health", timeout=0.5):
                pass
        except (urllib.error.URLError, TimeoutError):
            if attempts < 12:
                QtCore.QTimer.singleShot(400, lambda: self._open_dashboard_when_ready(attempts + 1))
                return
            self.status_label.setText("Dashboard failed to start; inspect the process output.")
            return
        QtGui.QDesktopServices.openUrl(QtCore.QUrl(url))
        self.status_label.setText(f"Dashboard running at {url}")

    @QtCore.Slot()
    def start_viewer(self) -> None:
        try:
            export = export_viewer_json(self.settings.export_path, self.settings.database_path)
            if export.stat().st_size < 100:
                raise RuntimeError("No transfer plans are available yet.")
        except Exception as error:
            QtWidgets.QMessageBox.warning(self, "Viewer unavailable", str(error))
            return
        self._start_module(self.viewer_process, "app.neo_viewer_qt", [str(export)])
        self.status_label.setText("Orbit viewer launched.")

    @QtCore.Slot()
    def start_update(self) -> None:
        if self.update_process.state() != QtCore.QProcess.ProcessState.NotRunning:
            return
        self.update_button.setEnabled(False)
        self.update_button.setText("Updating…")
        self.status_label.setText("Fetching JPL data and computing transfer plans…")
        self._start_module(self.update_process, "app.ingest")

    def _update_finished(self, exit_code: int) -> None:
        self.update_button.setEnabled(True)
        self.update_button.setText("Update from JPL")
        if exit_code == 0:
            self.refresh_status()
        else:
            self.status_label.setText("Update failed; inspect the process output.")

    def _capture_process_error(self, process: QtCore.QProcess) -> None:
        message = bytes(process.readAllStandardError()).decode("utf-8", errors="replace").strip()
        if message:
            self.status_label.setToolTip(message[-3000:])

    def refresh_status(self) -> None:
        stats = dashboard_stats(self.settings.database_path)
        run = latest_ingestion(self.settings.database_path)
        timestamp = run.get("finished_at") if run else "never"
        self.status_label.setText(
            f"{stats.get('asteroid_count', 0)} objects  ·  "
            f"{stats.get('plan_count', 0)} transfer plans  ·  Last update: {timestamp}"
        )

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        if self.dashboard_process.state() != QtCore.QProcess.ProcessState.NotRunning:
            self.dashboard_process.terminate()
            self.dashboard_process.waitForFinished(1500)
        super().closeEvent(event)


STYLESHEET = """
QWidget#root { background: #090d14; color: #edf4ff; }
QLabel#eyebrow, QLabel#cardKicker { color: #56d7ef; font-size: 11px; font-weight: 800; letter-spacing: 2px; }
QLabel#title { color: #edf4ff; font-size: 38px; font-weight: 800; }
QLabel#subtitle { color: #8796ab; font-size: 15px; }
QFrame#card { background: #111a26; border: 1px solid #253246; border-radius: 14px; }
QLabel#cardTitle { color: #edf4ff; font-size: 22px; font-weight: 750; }
QLabel#cardBody { color: #8796ab; font-size: 14px; }
QPushButton { background: #56d7ef; color: #071015; border: 0; border-radius: 8px; padding: 12px 18px; font-weight: 800; }
QPushButton:hover { background: #78e3f4; }
QPushButton:disabled { background: #30414f; color: #788894; }
QFrame#statusBar { background: #0d141e; border: 1px solid #253246; border-radius: 11px; }
QLabel#status { color: #aebed1; }
QPushButton#secondaryButton { background: transparent; color: #edf4ff; border: 1px solid #33445b; }
"""


def main() -> int:
    application = QtWidgets.QApplication(sys.argv)
    application.setApplicationName("Asteroid Intercept Planner")
    window = Launchpad()
    window.show()
    return application.exec()


if __name__ == "__main__":
    sys.exit(main())
