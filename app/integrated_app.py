"""Single-window desktop shell for the 3D viewer and native data explorer."""

from __future__ import annotations

import csv
import sys
import threading
import urllib.parse
from pathlib import Path

from PySide6 import QtCore, QtGui, QtWidgets
from werkzeug.serving import BaseWSGIServer, make_server

from app.config import ROOT_DIR, get_settings
from app.database import (
    dashboard_stats,
    latest_ingestion,
    list_asteroids,
    list_tables,
    run_readonly_query,
)
from app.export import export_viewer_json
from app.neo_viewer_qt import MissionViewerWidget
from app.web import create_app


class LocalDashboardServer:
    """Loopback-only WSGI server retained as the browser fallback."""

    def __init__(self) -> None:
        self.server: BaseWSGIServer = make_server("127.0.0.1", 0, create_app(), threaded=True)
        self.url = f"http://127.0.0.1:{self.server.server_port}"
        self.thread = threading.Thread(
            target=self.server.serve_forever,
            name="asteroid-dashboard",
            daemon=True,
        )

    def start(self) -> None:
        self.thread.start()

    def stop(self) -> None:
        self.server.shutdown()
        self.thread.join(timeout=2.0)


class NativeDataExplorer(QtWidgets.QWidget):
    """Native catalog and read-only SQL interface backed by SQLite."""

    CATALOG_COLUMNS = (
        ("designation", "Object"),
        ("close_approach_text", "Approach (TDB)"),
        ("distance_au", "Distance (AU)"),
        ("relative_velocity_kms", "Relative km/s"),
        ("eccentricity", "Eccentricity"),
        ("inclination_deg", "Inclination"),
        ("total_dv_kms", "Total Δv km/s"),
        ("tof_days", "TOF days"),
    )

    def __init__(self, database_path: Path, browser_url: str):
        super().__init__()
        self.database_path = database_path
        self.browser_url = browser_url
        self.query_columns: list[str] = []
        self.query_rows: list[dict] = []
        self.metric_labels: dict[str, QtWidgets.QLabel] = {}
        self._build_ui()
        self.refresh()

    def _build_ui(self) -> None:
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 20)
        layout.setSpacing(14)

        heading_row = QtWidgets.QHBoxLayout()
        heading = QtWidgets.QLabel("SQLITE DATA EXPLORER")
        heading.setObjectName("dataTitle")
        subheading = QtWidgets.QLabel("Local database · read-only query tools")
        subheading.setObjectName("dataSubtitle")
        heading_stack = QtWidgets.QVBoxLayout()
        heading_stack.addWidget(heading)
        heading_stack.addWidget(subheading)
        heading_row.addLayout(heading_stack)
        heading_row.addStretch(1)
        refresh_button = QtWidgets.QPushButton("Refresh")
        browser_button = QtWidgets.QPushButton("Open browser fallback")
        refresh_button.clicked.connect(self.refresh)
        browser_button.clicked.connect(self.open_browser)
        heading_row.addWidget(refresh_button)
        heading_row.addWidget(browser_button)
        layout.addLayout(heading_row)

        metrics = QtWidgets.QHBoxLayout()
        for key, label in (
            ("asteroid_count", "TRACKED OBJECTS"),
            ("approach_count", "CLOSE APPROACHES"),
            ("plan_count", "TRANSFER PLANS"),
            ("nearest_distance_au", "NEAREST ENCOUNTER"),
            ("lowest_total_dv_kms", "LOWEST ΔV"),
        ):
            card = QtWidgets.QFrame()
            card.setObjectName("metricCard")
            card_layout = QtWidgets.QVBoxLayout(card)
            card_layout.setContentsMargins(14, 11, 14, 11)
            name = QtWidgets.QLabel(label)
            name.setObjectName("metricName")
            value = QtWidgets.QLabel("—")
            value.setObjectName("metricValue")
            self.metric_labels[key] = value
            card_layout.addWidget(name)
            card_layout.addWidget(value)
            metrics.addWidget(card, 1)
        layout.addLayout(metrics)

        self.data_tabs = QtWidgets.QTabWidget()
        self.data_tabs.setObjectName("dataTabs")
        self.data_tabs.addTab(self._build_catalog_page(), "CATALOG")
        self.data_tabs.addTab(self._build_sql_page(), "SQL QUERY")
        layout.addWidget(self.data_tabs, 1)

    def _build_catalog_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(page)
        layout.setContentsMargins(0, 12, 0, 0)
        filters = QtWidgets.QHBoxLayout()
        self.catalog_search = QtWidgets.QLineEdit()
        self.catalog_search.setPlaceholderText("Search designation or name…")
        self.catalog_search.setClearButtonEnabled(True)
        self.catalog_search.textChanged.connect(self.refresh_catalog)
        self.catalog_sort = QtWidgets.QComboBox()
        for label, value in (
            ("Approach date", "date"),
            ("Miss distance", "distance"),
            ("Relative speed", "velocity"),
            ("Transfer Δv", "dv"),
            ("Inclination", "inclination"),
            ("Designation", "designation"),
        ):
            self.catalog_sort.addItem(label, value)
        self.catalog_sort.currentIndexChanged.connect(self.refresh_catalog)
        details_button = QtWidgets.QPushButton("Open selected details")
        details_button.clicked.connect(self.open_selected_browser)
        filters.addWidget(self.catalog_search, 1)
        filters.addWidget(self.catalog_sort)
        filters.addWidget(details_button)
        layout.addLayout(filters)

        self.catalog_table = QtWidgets.QTableWidget()
        self.catalog_table.setColumnCount(len(self.CATALOG_COLUMNS))
        self.catalog_table.setHorizontalHeaderLabels(
            [title for _key, title in self.CATALOG_COLUMNS]
        )
        self._configure_table(self.catalog_table)
        self.catalog_table.doubleClicked.connect(self.open_selected_browser)
        layout.addWidget(self.catalog_table, 1)
        self.catalog_status = QtWidgets.QLabel()
        self.catalog_status.setObjectName("tableStatus")
        layout.addWidget(self.catalog_status)
        return page

    def _build_sql_page(self) -> QtWidgets.QWidget:
        page = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(page)
        layout.setContentsMargins(0, 12, 0, 0)
        tables = ", ".join(row["name"] for row in list_tables(self.database_path))
        table_hint = QtWidgets.QLabel(f"AVAILABLE TABLES / VIEWS   {tables}")
        table_hint.setObjectName("tableHint")
        layout.addWidget(table_hint)

        self.query_editor = QtWidgets.QPlainTextEdit()
        self.query_editor.setObjectName("queryEditor")
        self.query_editor.setPlainText(
            "SELECT designation, close_approach_text, distance_au, "
            "total_dv_kms\nFROM latest_asteroid_summary\n"
            "ORDER BY close_approach_jd_tdb\nLIMIT 100"
        )
        self.query_editor.setMaximumHeight(135)
        layout.addWidget(self.query_editor)

        actions = QtWidgets.QHBoxLayout()
        run_button = QtWidgets.QPushButton("Run read-only query")
        export_button = QtWidgets.QPushButton("Export results to CSV")
        run_button.clicked.connect(self.run_query)
        export_button.clicked.connect(self.export_query)
        actions.addWidget(run_button)
        actions.addWidget(export_button)
        actions.addStretch(1)
        layout.addLayout(actions)

        self.query_table = QtWidgets.QTableWidget()
        self._configure_table(self.query_table)
        layout.addWidget(self.query_table, 1)
        self.query_status = QtWidgets.QLabel("Ready")
        self.query_status.setObjectName("tableStatus")
        layout.addWidget(self.query_status)
        return page

    def _configure_table(self, table: QtWidgets.QTableWidget) -> None:
        table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        table.setAlternatingRowColors(True)
        table.setSortingEnabled(False)
        table.verticalHeader().setVisible(False)
        table.horizontalHeader().setStretchLastSection(True)

    @QtCore.Slot()
    def refresh(self) -> None:
        stats = dashboard_stats(self.database_path)
        self.metric_labels["asteroid_count"].setText(str(stats.get("asteroid_count") or 0))
        self.metric_labels["approach_count"].setText(str(stats.get("approach_count") or 0))
        self.metric_labels["plan_count"].setText(str(stats.get("plan_count") or 0))
        nearest = stats.get("nearest_distance_au")
        lowest_dv = stats.get("lowest_total_dv_kms")
        self.metric_labels["nearest_distance_au"].setText(
            f"{nearest:.5f} AU" if nearest is not None else "—"
        )
        self.metric_labels["lowest_total_dv_kms"].setText(
            f"{lowest_dv:.2f} km/s" if lowest_dv is not None else "—"
        )
        run = latest_ingestion(self.database_path)
        timestamp = run.get("finished_at") if run else "never"
        self.setToolTip(f"Latest ingestion: {timestamp}")
        self.refresh_catalog()

    @QtCore.Slot()
    @QtCore.Slot(str)
    @QtCore.Slot(int)
    def refresh_catalog(self, _value: object = None) -> None:
        rows = list_asteroids(
            search=self.catalog_search.text(),
            sort=str(self.catalog_sort.currentData()),
            path=self.database_path,
        )
        self.catalog_table.setRowCount(len(rows))
        for row_index, row in enumerate(rows):
            for column_index, (key, _title) in enumerate(self.CATALOG_COLUMNS):
                value = row.get(key)
                if isinstance(value, float):
                    precision = 6 if key == "distance_au" else 3
                    text = f"{value:.{precision}f}"
                else:
                    text = "—" if value is None else str(value)
                item = QtWidgets.QTableWidgetItem(text)
                if column_index == 0:
                    item.setData(QtCore.Qt.ItemDataRole.UserRole, row["designation"])
                self.catalog_table.setItem(row_index, column_index, item)
        self.catalog_table.resizeColumnsToContents()
        self.catalog_status.setText(f"{len(rows)} objects shown")

    @QtCore.Slot()
    def run_query(self) -> None:
        try:
            columns, rows = run_readonly_query(
                self.query_editor.toPlainText(), path=self.database_path
            )
        except Exception as error:
            QtWidgets.QMessageBox.warning(self, "Query rejected", str(error))
            self.query_status.setText(str(error))
            return
        self.query_columns = columns
        self.query_rows = rows
        self.query_table.setColumnCount(len(columns))
        self.query_table.setHorizontalHeaderLabels(columns)
        self.query_table.setRowCount(len(rows))
        for row_index, row in enumerate(rows):
            for column_index, column in enumerate(columns):
                value = row.get(column)
                self.query_table.setItem(
                    row_index,
                    column_index,
                    QtWidgets.QTableWidgetItem("—" if value is None else str(value)),
                )
        self.query_table.resizeColumnsToContents()
        self.query_status.setText(f"{len(rows)} rows returned")

    @QtCore.Slot()
    def export_query(self) -> None:
        if not self.query_rows:
            self.run_query()
        if not self.query_rows:
            return
        filename, _selected_filter = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Export query results",
            "asteroid-query.csv",
            "CSV files (*.csv)",
        )
        if not filename:
            return
        with Path(filename).open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=self.query_columns)
            writer.writeheader()
            writer.writerows(self.query_rows)
        self.query_status.setText(f"Exported {len(self.query_rows)} rows to {filename}")

    @QtCore.Slot()
    def open_browser(self) -> None:
        QtGui.QDesktopServices.openUrl(QtCore.QUrl(self.browser_url))

    @QtCore.Slot()
    @QtCore.Slot(QtCore.QModelIndex)
    def open_selected_browser(self, _index: object = None) -> None:
        selected = self.catalog_table.selectionModel().selectedRows()
        if not selected:
            return
        item = self.catalog_table.item(selected[0].row(), 0)
        designation = str(item.data(QtCore.Qt.ItemDataRole.UserRole))
        encoded = urllib.parse.quote(designation, safe="")
        QtGui.QDesktopServices.openUrl(QtCore.QUrl(f"{self.browser_url}/asteroid/{encoded}"))


class IntegratedWindow(QtWidgets.QMainWindow):
    def __init__(
        self,
        viewer_export: Path,
        dashboard_server: LocalDashboardServer,
        database_path: Path,
    ) -> None:
        super().__init__()
        self.viewer_export = viewer_export
        self.database_path = database_path
        self.dashboard_server = dashboard_server
        self.update_process = QtCore.QProcess(self)
        self.update_errors = ""
        self.setWindowTitle("Asteroid Intercept Planner")
        self.setMinimumSize(1220, 780)
        self.resize(1500, 920)

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setObjectName("mainTabs")
        self.tabs.setDocumentMode(True)
        self.tabs.setTabPosition(QtWidgets.QTabWidget.TabPosition.North)
        self.tabs.tabBar().setExpanding(True)
        self.setCentralWidget(self.tabs)

        self.viewer = MissionViewerWidget(self.viewer_export)
        self.data_explorer = NativeDataExplorer(self.database_path, self.dashboard_server.url)
        self.tabs.addTab(self.viewer, "MISSION VIEWER")
        self.tabs.addTab(self.data_explorer, "DATA EXPLORER")

        toolbar = self.addToolBar("Application")
        toolbar.setObjectName("appToolbar")
        toolbar.setMovable(False)
        update_action = toolbar.addAction("Update from JPL")
        browser_action = toolbar.addAction("Open browser fallback")
        update_action.triggered.connect(self.update_from_jpl)
        browser_action.triggered.connect(self.data_explorer.open_browser)
        self.update_action = update_action

        self.update_process.readyReadStandardError.connect(self._capture_update_error)
        self.update_process.finished.connect(self._update_finished)
        self.statusBar().showMessage("Ready — viewer and SQLite explorer are running locally.")
        self.setStyleSheet(APP_STYLESHEET)

    @QtCore.Slot()
    def update_from_jpl(self) -> None:
        if self.update_process.state() != QtCore.QProcess.ProcessState.NotRunning:
            return
        self.update_errors = ""
        self.update_action.setEnabled(False)
        self.statusBar().showMessage("Fetching JPL data and recomputing interception plans…")
        environment = QtCore.QProcessEnvironment.systemEnvironment()
        environment.insert("PYTHONPATH", str(ROOT_DIR))
        self.update_process.setProcessEnvironment(environment)
        self.update_process.setWorkingDirectory(str(ROOT_DIR))
        self.update_process.start(sys.executable, ["-m", "app.ingest"])

    @QtCore.Slot()
    def _capture_update_error(self) -> None:
        self.update_errors += bytes(self.update_process.readAllStandardError()).decode(
            "utf-8", errors="replace"
        )

    @QtCore.Slot(int, QtCore.QProcess.ExitStatus)
    def _update_finished(self, exit_code: int, exit_status: QtCore.QProcess.ExitStatus) -> None:
        self.update_action.setEnabled(True)
        if exit_code != 0 or exit_status == QtCore.QProcess.ExitStatus.CrashExit:
            QtWidgets.QMessageBox.critical(
                self,
                "JPL update failed",
                self.update_errors[-4000:] or self.update_process.errorString(),
            )
            self.statusBar().showMessage("JPL update failed.")
            return
        export_viewer_json(self.viewer_export, self.database_path)
        replacement = MissionViewerWidget(self.viewer_export)
        old_viewer = self.viewer
        self.viewer = replacement
        self.tabs.removeTab(0)
        self.tabs.insertTab(0, replacement, "MISSION VIEWER")
        old_viewer.deleteLater()
        self.data_explorer.refresh()
        self.tabs.setCurrentIndex(0)
        self.statusBar().showMessage("JPL data and mission plans updated.", 8000)

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        if self.update_process.state() != QtCore.QProcess.ProcessState.NotRunning:
            self.update_process.terminate()
            self.update_process.waitForFinished(1500)
        self.dashboard_server.stop()
        super().closeEvent(event)


APP_STYLESHEET = """
QMainWindow, QWidget { background: #070b12; color: #e7f0fa; }
QTabWidget#mainTabs::pane, QTabWidget#dataTabs::pane { border: 0; background: #070b12; }
QTabBar::tab { background: #0d141e; color: #74879c; border: 0; border-right: 1px solid #253446; padding: 12px 26px; font-size: 11px; font-weight: 800; letter-spacing: 1px; }
QTabBar::tab:selected { color: #071015; background: #52d6ec; }
QTabBar::tab:hover:!selected { color: #e7f0fa; background: #14202d; }
QToolBar#appToolbar { background: #0a1018; border-bottom: 1px solid #253446; spacing: 5px; padding: 5px 10px; }
QToolBar#appToolbar QToolButton { color: #a9bacb; background: #111c29; border: 1px solid #2b3d52; border-radius: 5px; padding: 6px 11px; }
QToolBar#appToolbar QToolButton:hover { color: #ffffff; border-color: #52d6ec; }
QStatusBar { color: #778a9f; background: #090f16; border-top: 1px solid #253446; }
QLabel#dataTitle { color: #f0f6fc; font-size: 24px; font-weight: 800; }
QLabel#dataSubtitle, QLabel#tableStatus, QLabel#tableHint { color: #718499; font-size: 10px; }
QFrame#metricCard { background: #101925; border: 1px solid #263446; border-radius: 8px; }
QLabel#metricName { color: #718499; font-size: 9px; font-weight: 800; letter-spacing: 1px; }
QLabel#metricValue { color: #f0f6fc; font-size: 19px; font-weight: 800; }
QPushButton { background: #52d6ec; color: #071015; border: 0; border-radius: 6px; padding: 8px 12px; font-weight: 800; }
QPushButton:hover { background: #78e3f4; }
QLineEdit, QComboBox, QPlainTextEdit { background: #0a1018; border: 1px solid #2b3d52; border-radius: 6px; padding: 7px 9px; color: #e7f0fa; }
QPlainTextEdit#queryEditor { font-family: "SF Mono", monospace; font-size: 11px; }
QTableWidget { background: #090f17; alternate-background-color: #0d1621; color: #c4d2df; border: 1px solid #263446; gridline-color: #1d2a3a; selection-background-color: #17475c; selection-color: #ffffff; }
QHeaderView::section { background: #101925; color: #7f92a7; border: 0; border-right: 1px solid #263446; border-bottom: 1px solid #263446; padding: 7px; font-size: 9px; font-weight: 800; }
"""


def main() -> int:
    settings = get_settings()
    export = export_viewer_json(settings.export_path, settings.database_path)
    server = LocalDashboardServer()
    server.start()
    application = QtWidgets.QApplication(sys.argv)
    application.setApplicationName("Asteroid Intercept Planner")
    window = IntegratedWindow(export, server, settings.database_path)
    window.show()
    window.raise_()
    window.activateWindow()
    return application.exec()


if __name__ == "__main__":
    sys.exit(main())
