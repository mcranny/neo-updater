#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, sys
from pathlib import Path

from scripts.sbdb_client import discover_earth_close_approaches
from scripts.hud_metrics import compute_hud

def main():
    ap = argparse.ArgumentParser(description="Update Earth close-approach feed (CAD-first).")
    ap.add_argument("--out", default="data/hazardous_neos/latest.json")
    ap.add_argument("--date-min", default="now")
    ap.add_argument("--date-max", default="+60")
    ap.add_argument("--dist-max", default="0.05")  # AU or e.g. 10LD
    ap.add_argument("--no-elements", action="store_true")
    ap.add_argument("--limit", type=int, default=2000)
    ap.add_argument("--with-intercepts", action="store_true")
    ap.add_argument("--leo-alt-m", type=float, default=500e3)
    args = ap.parse_args()

    # 1) Pull CAD + SBDB
    records = discover_earth_close_approaches(
        date_min=args.date_min,
        date_max=args.date_max,
        dist_max=args.dist_max,
        limit=args.limit,
        attach_elements=(not args.no_elements),
    )

    # 2) Attach HUD metrics per object (at the next CA time)
    for obj in records:
        el = obj.get("elements")
        ca = obj.get("next_ca")
        hud = None
        if el and ca and ca.get("jd_tdb") is not None:
            try:
                hud = compute_hud(el, float(ca["jd_tdb"]))
                hud["timestamp_cd_tdb"] = ca.get("cd_tdb")  # for display
            except Exception:
                hud = None
        obj["hud"] = hud

    # 3) Package payload
    payload = {
        "schema": "close_approach_feed/1.1",
        "source": "JPL SSD CAD + SBDB",
        "query": {"date_min": args.date_min, "date_max": args.date_max, "dist_max": args.dist_max},
        "count": len(records),
        "objects": records,
    }

    # 4) (Optional) Attach Lambert intercepts
    if args.with_intercepts:
        # relative import so it works when run as a module
        from .plan_intercepts import attach_intercepts
        payload = attach_intercepts(payload, leo_alt_m=args.leo_alt_m)

    # 5) Write file
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2))
    print(f"Wrote {out} with {len(records)} objects"
          f"{' + intercepts' if args.with_intercepts else ''}.")

if __name__ == "__main__":
    sys.exit(main())
