import os, json, base64, datetime, pathlib
from scripts.gh_client import GitHubClient
from scripts.neo_fetch import fetch_hazardous_neos

def _decode_existing_local(path: str):
    p = pathlib.Path(path)
    if not p.exists():
        return {}
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}

def _merge(existing: dict, *, date_iso: str, now_utc: datetime.datetime, pha: dict) -> dict:
    merged = dict(existing) if isinstance(existing, dict) else {}
    merged.update({
        "date_utc": date_iso,
        "snapshot_utc": now_utc.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "count": pha.get("count"),
        "potentially_hazardous_neos": pha.get("potentially_hazardous_neos", []),
    })
    return merged

def main():
    out_path   = os.environ.get("OUT_FILE", "data/hazardous_neos/latest.json")
    nasa_key   = os.environ.get("NASA_API_KEY", "DEMO_KEY")

    now_utc = datetime.datetime.now(datetime.timezone.utc)
    date_override = os.environ.get("DATE_ISO")
    date_iso = date_override or now_utc.date().isoformat()
    pha = fetch_hazardous_neos(date_iso, nasa_key)

    existing = _decode_existing_local(out_path)
    merged = _merge(existing, date_iso=date_iso, now_utc=now_utc, pha=pha)

    pathlib.Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(merged, f, ensure_ascii=False, indent=2)

    print(f"[LOCAL RUN] wrote {out_path}  date={date_iso} count={pha.get('count')}")

if __name__ == "__main__":
    main()
