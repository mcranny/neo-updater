"""Flask dashboard for browsing asteroid and transfer-plan data."""

from __future__ import annotations

import csv
import io

from flask import Flask, Response, jsonify, redirect, render_template, request, url_for

from app.config import get_settings
from app.database import (
    asteroid_approaches,
    dashboard_stats,
    get_asteroid,
    initialize_database,
    latest_ingestion,
    list_asteroids,
    list_tables,
    run_readonly_query,
)


def create_app() -> Flask:
    settings = get_settings()
    initialize_database(settings.database_path)
    application = Flask(__name__)
    application.config.update(SECRET_KEY=settings.secret_key)

    @application.template_filter("number")
    def format_number(value: object, digits: int = 3) -> str:
        if value is None:
            return "—"
        try:
            return f"{float(value):,.{digits}f}"
        except (TypeError, ValueError):
            return str(value)

    @application.route("/")
    def index():
        search = request.args.get("q", "").strip()
        sort = request.args.get("sort", "date")
        order = request.args.get("order", "asc")
        return render_template(
            "index.html",
            rows=list_asteroids(search, sort, order, path=settings.database_path),
            stats=dashboard_stats(settings.database_path),
            latest_run=latest_ingestion(settings.database_path),
            search=search,
            sort=sort,
            order=order,
        )

    @application.route("/asteroids")
    def asteroids_redirect():
        return redirect(url_for("index"))

    @application.route("/asteroid/<designation>")
    def asteroid_detail(designation: str):
        asteroid = get_asteroid(designation, settings.database_path)
        if asteroid is None:
            return render_template("not_found.html", designation=designation), 404
        return render_template(
            "detail.html",
            asteroid=asteroid,
            approaches=asteroid_approaches(designation, settings.database_path),
        )

    @application.route("/data", methods=["GET", "POST"])
    def data_explorer():
        query = request.form.get(
            "query",
            "SELECT * FROM latest_asteroid_summary ORDER BY close_approach_jd_tdb LIMIT 100",
        )
        columns: list[str] = []
        rows: list[dict] = []
        error = None
        if request.method == "POST":
            try:
                columns, rows = run_readonly_query(query, path=settings.database_path)
            except Exception as exception:
                error = str(exception)
        return render_template(
            "data.html",
            tables=list_tables(settings.database_path),
            query=query,
            columns=columns,
            rows=rows,
            error=error,
        )

    @application.post("/data/export")
    def export_query():
        query = request.form.get("query", "")
        try:
            columns, rows = run_readonly_query(query, path=settings.database_path)
        except ValueError as exception:
            return str(exception), 400
        stream = io.StringIO()
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
        return Response(
            stream.getvalue(),
            mimetype="text/csv",
            headers={"Content-Disposition": "attachment; filename=asteroid-query.csv"},
        )

    @application.get("/api/asteroids")
    def asteroid_api():
        return jsonify(
            list_asteroids(
                request.args.get("q", ""),
                request.args.get("sort", "date"),
                request.args.get("order", "asc"),
                path=settings.database_path,
            )
        )

    @application.get("/health")
    def health():
        return {"status": "ok", "database": str(settings.database_path)}

    return application


app = create_app()


def main() -> None:
    settings = get_settings()
    app.run(
        host=settings.flask_host,
        port=settings.flask_port,
        debug=settings.flask_debug,
        use_reloader=False,
    )


if __name__ == "__main__":
    main()
