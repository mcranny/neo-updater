import json, urllib.request, urllib.parse
API_NASA = "https://api.nasa.gov/neo/rest/v1/feed"

def fetch_hazardous_neos(date_iso: str, api_key: str, user_agent: str = "neo-snapshot/1.0"):
    params = {"start_date": date_iso, "end_date": date_iso, "api_key": api_key}
    url = f"{API_NASA}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"User-Agent": user_agent}, method="GET")
    with urllib.request.urlopen(req, timeout=30) as r:
        if r.status != 200:
            raise RuntimeError(f"NeoWs {r.status}")
        data = json.loads(r.read().decode("utf-8"))

    neos = data.get("near_earth_objects", {}).get(date_iso, [])
    out = []
    for obj in neos:
        if not obj.get("is_potentially_hazardous_asteroid"):
            continue
        diam = obj.get("estimated_diameter", {}).get("meters", {})
        dmin, dmax = diam.get("estimated_diameter_min"), diam.get("estimated_diameter_max")
        ca = next((c for c in obj.get("close_approach_data", []) if c.get("close_approach_date") == date_iso), None)
        rec = {
            "name": obj.get("name"),
            "neo_reference_id": obj.get("neo_reference_id"),
            "absolute_magnitude_h": obj.get("absolute_magnitude_h"),
            "estimated_diameter_m_min": dmin,
            "estimated_diameter_m_max": dmax,
        }
        if ca:
            try:    miss_km = float(ca["miss_distance"]["kilometers"])
            except: miss_km = None
            try:    vel_kms = float(ca["relative_velocity"]["kilometers_per_second"])
            except: vel_kms = None
            rec["close_approach"] = {
                "date_full": ca.get("close_approach_date_full"),
                "orbiting_body": ca.get("orbiting_body"),
                "miss_distance_km": miss_km,
                "relative_velocity_km_s": vel_kms
            }
        out.append(rec)
    return {"count": len(out), "potentially_hazardous_neos": out}
