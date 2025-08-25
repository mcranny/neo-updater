# scripts/sbdb_client.py
from __future__ import annotations
import math, requests
from typing import Dict, Any, Optional

_SBDB = "https://ssd-api.jpl.nasa.gov/sbdb.api"

def fetch_elements(sstr: str) -> Optional[Dict[str, Any]]:
    """
    Return SBDB osculating elements for a small body.
    Keys: a_AU,e,i_deg,raan_deg,argp_deg,M_deg,epoch_jd,source
    """
    r = requests.get(_SBDB, params={"sstr": sstr, "phys-par": "0"}, timeout=30)
    r.raise_for_status()
    j = r.json()

    # SBDB packs elements under orbit.elements as list of {name,value,units,...}
    # Field names follow: a,e,i,om,w,ma + epoch (JD) .  See SBDB docs.  :contentReference[oaicite:2]{index=2}
    orbit = (j.get("orbit") or {})
    elems = {e.get("name"): float(e.get("value"))
             for e in (orbit.get("elements") or []) if "value" in e}
    if not elems:
        return None

    epoch = orbit.get("epoch")
    if isinstance(epoch, dict):
        epoch = epoch.get("value")

    out = {
        "a_AU": elems.get("a"),
        "e": elems.get("e"),
        "i_deg": elems.get("i"),
        "raan_deg": elems.get("om"),
        "argp_deg": elems.get("w"),
        "M_deg": elems.get("ma"),
        "epoch_jd": float(epoch) if epoch is not None else None,
        "source": "JPL SBDB",
    }
    # sanity
    if any(out[k] is None for k in ("a_AU","e","i_deg","raan_deg","argp_deg","M_deg")):
        return None
    return out
