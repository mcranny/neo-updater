from __future__ import annotations

from app.database import upsert_payload


def test_dashboard_and_api(monkeypatch, tmp_path, sample_payload):
    database = tmp_path / "asteroids.db"
    monkeypatch.setenv("ASTEROID_DATABASE_PATH", str(database))
    upsert_payload(sample_payload, database)

    from app.web import create_app

    application = create_app()
    application.config.update(TESTING=True)
    client = application.test_client()

    dashboard = client.get("/")
    assert dashboard.status_code == 200
    assert b"TEST-1" in dashboard.data
    assert b"5.30" in dashboard.data

    detail = client.get("/asteroid/TEST-1")
    assert detail.status_code == 200
    assert b"Synthetic test asteroid" in detail.data

    api = client.get("/api/asteroids")
    assert api.status_code == 200
    assert api.get_json()[0]["designation"] == "TEST-1"

    health = client.get("/health")
    assert health.get_json()["status"] == "ok"
