#!/usr/bin/env python3
"""
NEO Intercept Planner (JSON → JSON)

Reads a hazardous-NEO snapshot JSON (default: OUT_FILE or data/hazardous_neos/latest.json)
and writes a separate planning file (default: sibling *_intercept.json).

For each NEO, computes:
- Circular/coplanar Hohmann transfer from r1=1 AU → r2 (configurable / inferred)
- Time of flight and v_inf at Earth
- LEO Δv using v_inf and Oberth-aware escape from circular LEO
- Simple arrival Δv (heliocentric insertion to target circular speed)
- Heuristic on-asteroid orbit suggestion (>=3R and >=1 km altitude)
- Launch/arrival timestamps; roll-forward of past arrivals using synodic period

CLI:
  python -m scripts.neo_intercept_planner \
      --input data/hazardous_neos/latest.json \
      --output data/hazardous_neos/latest_intercept.json \
      --roll-past \
      --profile neosurveyor

Notes:
- If target orbital radius r2 is unavailable, defaults to 1.0 AU (Earth-like).
- If an arrival timestamp is in the past and --roll-past is set, we add k * P_syn
  (synodic period between Earth at 1 AU and target at r2) until arrival > now.
- Uses only circular/coplanar approximations; good for quick sizing, not mission design.
"""

from __future__ import annotations
import argparse
import json
import math
import os
import sys
from datetime import datetime, timezone, timedelta
from typing import Any, Dict, Optional, Tuple, List

# --- Constants (SI) ---
AU = 1.495978707e11           # m
MU_SUN = 1.32712440018e20     # m^3/s^2
MU_EARTH = 3.986004418e14     # m^3/s^2
R_EARTH = 6378137.0           # m
G0 = 9.80665                  # m/s^2

# Default LEO altitude
DEFAULT_LEO_ALT_M = 500e3     # 500 km

# Heuristic asteroid density (if unknown)
RHO_ASTEROID = 2000.0         # kg/m^3 (2 g/cc)

SCHEMA_VERSION = "1.1.0"


# -------------------- Utilities --------------------

def iso_z(dt: datetime) -> str:
    return dt.astimezone(timezone.utc).replace(tzinfo=timezone.utc).isoformat().replace("+00:00", "Z")


def parse_any_datetime(s: str) -> Optional[datetime]:
    """
    Accepts formats like:
      - "2025-Aug-24 17:55"
      - "2025-08-24T17:55:00Z" or offsets
      - "2025-08-24"
    Returns an aware datetime in UTC.
    """
    if not s:
        return None
    s = s.strip()
    fmts = [
        "%Y-%b-%d %H:%M",      # "2025-Aug-24 17:55"
        "%Y-%m-%dT%H:%M:%SZ",  # ISO Z
        "%Y-%m-%dT%H:%M:%S.%fZ",
        "%Y-%m-%dT%H:%M:%S%z",
        "%Y-%m-%d %H:%M:%S",
        "%Y-%m-%d",
    ]
    for f in fmts:
        try:
            dt = datetime.strptime(s, f)
            if dt.tzinfo is None:
                dt = dt.replace(tzinfo=timezone.utc)
            else:
                dt = dt.astimezone(timezone.utc)
            return dt
        except Exception:
            pass
    try:
        dt = datetime.fromisoformat(s.replace("Z", "+00:00"))
        return dt.astimezone(timezone.utc)
    except Exception:
        return None


def hohmann_delta_v_and_tof(r1_m: float, r2_m: float, mu: float = MU_SUN) -> Tuple[float, float, float, float]:
    """
    Hohmann transfer between two circular, coplanar orbits.
    Returns (dv1, dv2, tof_sec, v_inf_depart_approx)

    dv1: departure heliocentric Δv at r1 (m/s)
    dv2: arrival heliocentric Δv at r2 (m/s)
    tof_sec: π * sqrt(a_t^3 / μ)
    v_inf_depart_approx: magnitude of departure Δv relative to planetary circular speed (≈ v_inf)
    """
    v1 = math.sqrt(mu / r1_m)
    v2 = math.sqrt(mu / r2_m)
    a_t = 0.5 * (r1_m + r2_m)
    v_p = math.sqrt(mu * (2.0 / r1_m - 1.0 / a_t))  # on transfer at perihelion (r1)
    v_a = math.sqrt(mu * (2.0 / r2_m - 1.0 / a_t))  # on transfer at aphelion (r2)
    dv1 = abs(v_p - v1)
    dv2 = abs(v2 - v_a)
    tof = math.pi * math.sqrt(a_t ** 3 / mu)
    v_inf = dv1  # patched-conic approx: dv1 ≈ v_inf wrt Earth
    return dv1, dv2, tof, v_inf


def leo_escape_dv(v_inf: float, leo_alt_m: float = DEFAULT_LEO_ALT_M) -> float:
    """
    Δv from circular LEO to achieve hyperbolic excess v_inf leaving Earth.
    Uses: Δv_LEO = sqrt(v_inf^2 + v_esc^2) - v_circ, where v_esc = sqrt(2)*v_circ.
    """
    r_park = R_EARTH + leo_alt_m
    v_circ = math.sqrt(MU_EARTH / r_park)
    v_esc = math.sqrt(2.0) * v_circ
    return math.sqrt(v_inf**2 + v_esc**2) - v_circ


def synodic_period_days(r1_m: float, r2_m: float, mu: float = MU_SUN) -> float:
    """
    Synodic period between two circular orbits: P_syn = 2π / |n1 - n2|
    where n = sqrt(μ / r^3). Returns days.
    """
    n1 = math.sqrt(mu / (r1_m ** 3))
    n2 = math.sqrt(mu / (r2_m ** 3))
    if abs(n1 - n2) < 1e-15:
        return 1e9  # effectively infinite; avoid divide-by-zero
    psyn = 2.0 * math.pi / abs(n1 - n2)  # seconds
    return psyn / 86400.0


def roll_forward(arrival_utc: datetime, psyn_days: float, now_utc: datetime) -> Tuple[datetime, bool, int]:
    """
    If arrival_utc <= now_utc, add k * P_syn (in days) until it's future.
    Returns (rolled_arrival, rolled?, k)
    """
    if arrival_utc > now_utc:
        return arrival_utc, False, 0
    if psyn_days > 1e8:  # guard for effectively-infinite synodic
        k = math.ceil((now_utc - arrival_utc).total_seconds() / 86400.0)
        return arrival_utc + timedelta(days=k), True, k
    delta_days = (now_utc - arrival_utc).total_seconds() / 86400.0
    k = math.ceil(delta_days / psyn_days)
    rolled = arrival_utc + timedelta(days=k * psyn_days)
    return rolled, True, k


def flatten_neows(blob: Dict[str, Any]) -> List[Dict[str, Any]]:
    """
    Accepts multiple snapshot shapes and returns a flat list of NEOs:
      - {'potentially_hazardous_neos': [...]}
      - {'neos': [...]}
      - {'near_earth_objects': {'YYYY-MM-DD': [...]}}  # feed
      - {'near_earth_objects': [...]}                  # browse
      - the blob itself might already be a list
    """
    # Preferred pipeline keys
    for k in ("potentially_hazardous_neos", "neos"):
        v = blob.get(k)
        if isinstance(v, list):
            return v

    neo = blob.get("near_earth_objects")
    if isinstance(neo, list):
        return neo
    if isinstance(neo, dict):
        flat: List[Dict[str, Any]] = []
        for day_list in neo.values():
            if isinstance(day_list, list):
                flat.extend(day_list)
        return flat

    # Sometimes the file IS the list
    if isinstance(blob, list):
        return blob

    # Loose fallbacks some feeds use
    for k in ("data", "objects"):
        v = blob.get(k)
        if isinstance(v, list):
            return v

    raise SystemExit("Could not find NEO list in input JSON.")


def get_estimated_diameter_m_bounds(neo: Dict[str, Any]) -> tuple[Optional[float], Optional[float]]:
    """
    Return (d_min_m, d_max_m) from either flattened fields or NeoWs nested units.
    """
    # Flattened (your prior shape)
    dmin = neo.get("estimated_diameter_m_min")
    dmax = neo.get("estimated_diameter_m_max")
    if dmin is not None and dmax is not None:
        try:
            return float(dmin), float(dmax)
        except Exception:
            pass

    # NeoWs nested shape
    ed = neo.get("estimated_diameter")
    if isinstance(ed, dict):
        meters = ed.get("meters") or {}
        dmin2 = meters.get("estimated_diameter_min")
        dmax2 = meters.get("estimated_diameter_max")
        try:
            if dmin2 is not None and dmax2 is not None:
                return float(dmin2), float(dmax2)
        except Exception:
            pass

    return None, None


def pick_target_r2_au(neo: Dict[str, Any], default_r2_au: float) -> float:
    """
    Try to infer a plausible r2 (AU). If snapshot includes orbital data with 'semi_major_axis', use it.
    Else fallback to default_r2_au.
    """
    for key in ("orbital_data", "orbit", "elements"):
        d = neo.get(key)
        if isinstance(d, dict):
            a = d.get("semi_major_axis") or d.get("a") or d.get("a_AU")
            try:
                if a:
                    return float(a)
            except Exception:
                pass
    return float(default_r2_au)


def estimate_asteroid_mass_and_mu(neo: Dict[str, Any]) -> Tuple[Optional[float], Optional[float], Optional[float]]:
    """
    Estimate radius (m), mass (kg), and mu=G*M (m^3/s^2) from diameter bounds and density.
    Returns (R, M, mu) or (None, None, None) if insufficient data.
    """
    dmin, dmax = get_estimated_diameter_m_bounds(neo)
    if dmin is None or dmax is None:
        return None, None, None
    d_mean = 0.5 * (float(dmin) + float(dmax))  # meters
    r = 0.5 * d_mean
    volume = (4.0 / 3.0) * math.pi * r**3
    M = RHO_ASTEROID * volume
    G = 6.67430e-11
    mu = G * M
    return r, M, mu


# -------------------- Core Planning --------------------

def plan_for_neo(
    neo: Dict[str, Any],
    r1_au: float,
    default_r2_au: float,
    leo_alt_m: float,
    profile: Dict[str, Any],
    arrival_override: Optional[datetime],
    roll_past: bool,
    now_utc: datetime
) -> Dict[str, Any]:
    # r2 selection
    r2_au = pick_target_r2_au(neo, default_r2_au)
    r1_m = r1_au * AU
    r2_m = r2_au * AU

    # Hohmann
    dv1, dv2, tof_sec, v_inf = hohmann_delta_v_and_tof(r1_m, r2_m, MU_SUN)
    dv_leo = leo_escape_dv(v_inf, leo_alt_m)
    total_dv = dv_leo + dv2  # not including asteroid capture

    # Dates: try to find arrival in the NEO blob if not overridden
    arrival_dt = None
    if arrival_override:
        arrival_dt = arrival_override
    else:
        ca = neo.get("close_approach") or neo.get("close_approach_data")
        if isinstance(ca, dict):
            arrival_dt = parse_any_datetime(ca.get("date_full") or ca.get("close_approach_date_full") or ca.get("date"))
        elif isinstance(ca, list) and ca:
            earth_first = None
            for entry in ca:
                if entry.get("orbiting_body") == "Earth":
                    earth_first = entry
                    break
            cand = earth_first or ca[0]
            arrival_dt = parse_any_datetime(
                cand.get("date_full") or cand.get("close_approach_date_full") or cand.get("close_approach_date")
            )
    if arrival_dt is None:
        # Default to "now + TOF", i.e., launch now
        arrival_dt = now_utc + timedelta(seconds=tof_sec)

    # Roll forward if requested
    psyn_days = synodic_period_days(r1_m, r2_m, MU_SUN)
    rolled = False
    k_roll = 0
    if roll_past:
        arrival_dt, rolled, k_roll = roll_forward(arrival_dt, psyn_days, now_utc)

    # Launch = arrival - TOF
    depart_dt = arrival_dt - timedelta(seconds=tof_sec)

    # Simple on-asteroid orbit heuristic
    R_ast, M_ast, mu_ast = estimate_asteroid_mass_and_mu(neo)
    if R_ast is not None and mu_ast is not None and mu_ast > 0:
        alt = max(3.0 * R_ast, 1000.0)  # meters
        r_orb = R_ast + alt
        v_circ = math.sqrt(mu_ast / r_orb)
        T = (2.0 * math.pi * r_orb / v_circ) if v_circ > 0 else None
        suggested_orbit = {
            "radius_m": r_orb,
            "altitude_m": alt,
            "circular_speed_m_s": v_circ,
            "period_s": T,
            "note": "alt >= max(3R, 1 km); density 2 g/cc assumed",
        }
    else:
        suggested_orbit = {
            "radius_m": None,
            "altitude_m": 1000.0,
            "circular_speed_m_s": None,
            "period_s": None,
            "note": "Diameter missing; defaulting to >=1 km altitude heuristic.",
        }

    # Rocket equation prop estimate for total_dv (very rough)
    m0 = profile["m0_kg"]
    isp = profile["Isp_s"]
    if isp > 0 and total_dv > 0:
        mr = math.exp(total_dv / (isp * G0))  # mass ratio m0/mf
        mf = m0 / mr
        prop_used = max(m0 - mf, 0.0)
    else:
        prop_used = 0.0

    plan = {
        "schema_version": SCHEMA_VERSION,
        "profile": profile["name"],
        "r1_AU": r1_au,
        "r2_AU": r2_au,
        "tof_days": tof_sec / 86400.0,
        "synodic_days": psyn_days,
        "arrival_utc": iso_z(arrival_dt),
        "departure_utc": iso_z(depart_dt),
        "rolled_forward": rolled,
        "roll_periods_added": k_roll,
        "dv_depart_heliocentric_m_s": dv1,
        "dv_arrive_heliocentric_m_s": dv2,
        "v_inf_m_s": v_inf,
        "dv_from_LEO_m_s": dv_leo,
        "dv_total_m_s": total_dv,
        "leo_altitude_m": leo_alt_m,
        "spacecraft_mass_kg": m0,
        "Isp_s": isp,
        "propellant_estimate_kg": prop_used,
        "suggested_orbit": suggested_orbit,
        "notes": "Circular/coplanar Hohmann; Earth SOI patched-conic; arrival date rolled if in past and --roll-past set.",
    }
    return plan


# -------------------- IO & CLI --------------------

def load_json(path: str) -> Dict[str, Any]:
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def write_json(obj: Dict[str, Any], path: str) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, sort_keys=False)


def derive_default_output(in_path: str) -> str:
    d, base = os.path.split(in_path)
    name, _ = os.path.splitext(base)
    out_name = f"{name}_intercept.json"
    return os.path.join(d, out_name)


def build_profile(name: str, m0_kg: float, isp_s: float) -> Dict[str, Any]:
    return {"name": name, "m0_kg": float(m0_kg), "Isp_s": float(isp_s)}


def main():
    parser = argparse.ArgumentParser(description="Compute simple intercept plans for hazardous NEOs.")
    parser.add_argument("--input", default=os.getenv("OUT_FILE", "data/hazardous_neos/latest.json"))
    parser.add_argument("--output", default=None, help="Explicit output path (JSON). If omitted, writes *_intercept.json next to input.")
    parser.add_argument("--r1-au", type=float, default=1.0, help="Departure circular orbit radius around Sun (AU). Default 1.0 AU.")
    parser.add_argument("--r2-au", type=float, default=None, help="Override target circular orbit radius (AU) for ALL NEOs.")
    parser.add_argument("--leo-km", type=float, default=500.0, help="Parking orbit altitude in km. Default 500.")
    parser.add_argument("--arrive-utc", type=str, default=None, help="Force arrival time (ISO/Z or 'YYYY-Mon-DD HH:MM').")
    parser.add_argument("--roll-past", action="store_true", help="If arrival is in the past, roll it forward by whole synodic periods.")
    parser.add_argument("--profile", choices=["neosurveyor", "custom"], default="neosurveyor")
    parser.add_argument("--m0-kg", type=float, default=1300.0, help="Only for --profile custom.")
    parser.add_argument("--Isp-s", type=float, default=230.0, help="Only for --profile custom.")

    args = parser.parse_args()
    now_utc = datetime.now(timezone.utc)

    if args.profile == "neosurveyor":
        prof = build_profile("neosurveyor", 1300.0, 230.0)
    else:
        prof = build_profile("custom", args.m0_kg, args.Isp_s)

    in_path = args.input
    out_path = args.output or derive_default_output(in_path)
    r1_au = float(args.r1_au)
    default_r2_au = float(args.r2_au) if args.r2_au else 1.0
    leo_alt_m = float(args.leo_km) * 1e3
    arrival_override = parse_any_datetime(args.arrive_utc) if args.arrive_utc else None

    data = load_json(in_path)
    try:
        neos = flatten_neows(data)
    except SystemExit as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(2)

    plans: List[Dict[str, Any]] = []
    for neo in neos:
        try:
            plan = plan_for_neo(
                neo=neo,
                r1_au=r1_au,
                default_r2_au=default_r2_au,
                leo_alt_m=leo_alt_m,
                profile=prof,
                arrival_override=arrival_override,
                roll_past=args.roll_past,
                now_utc=now_utc
            )
            neo_out = dict(neo)
            neo_out["intercept_plan"] = plan
            plans.append(neo_out)
        except Exception as e:
            plans.append({
                "name": neo.get("name"),
                "neo_reference_id": neo.get("neo_reference_id") or neo.get("id"),
                "error": f"{type(e).__name__}: {e}"
            })

    out_obj: Dict[str, Any] = {
        "date_utc": now_utc.date().isoformat(),
        "snapshot_utc": iso_z(now_utc),
        "count": len(plans),
        "schema": {
            "version": SCHEMA_VERSION,
            "source": "neo_intercept_planner.py",
            "assumptions": {
                "circular_coplanar": True,
                "density_kg_m3": RHO_ASTEROID,
                "leo_alt_km": leo_alt_m / 1e3
            }
        },
        "potentially_hazardous_neos": plans
    }
    write_json(out_obj, out_path)
    print(f"Wrote intercept plans → {out_path}")


if __name__ == "__main__":
    main()
