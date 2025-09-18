#!/usr/bin/env bash
set -euo pipefail

# make sure we have fresh PHA data
python -m scripts.update_hazardous_neos

# try a small grid; expand if needed
for D in 50 60 70 80 90 100 110 120; do
  for T in 80 100 120 140 160 180 200 220; do
    OUT="/tmp/lambert_${D}_${T}.json"
    echo "Trying depart=${D}d TOF=${T}d"

    # run planner (fallback disabled in your file, so it will either succeed or cleanly fail)
    python -m scripts.neo_intercept_planner \
      --inp data/hazardous_neos/latest.json \
      --out "$OUT" \
      --with-sbdb \
      --polyline-n 600 \
      --depart-days "$D" \
      --tof-days "$T" \
      --debug || true

    # check if a Lambert polyline was produced
    ok="$(python - "$OUT" <<'PY'
import json, sys
try:
    j = json.load(open(sys.argv[1]))
    for r in j.get("potentially_hazardous_neos", []):
        ip = (r or {}).get("intercept_plan") or {}
        poly = ip.get("lambert_polyline_xyz_au")
        if poly and len(poly) >= 2:
            print(1)
            break
    else:
        print(0)
except Exception:
    print(0)
PY
)"
    if [ "$ok" = "1" ]; then
      echo "SUCCESS at depart=${D}d TOF=${T}d -> $OUT"
      cp "$OUT" data/hazardous_neos/latest_intercept.json
      exit 0
    fi
  done
done

echo "No valid Lambert arc found in the sweep."
exit 1
