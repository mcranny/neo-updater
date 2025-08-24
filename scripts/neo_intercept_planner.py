#!/usr/bin/env python3
"""
NEO Intercept Planner (JSON → JSON)
- Reads your NEO snapshot JSON, computes a 500 km LEO → heliocentric Hohmann transfer
  using a NEO Surveyor–like spacecraft by default, and appends an "intercept_plan"
  block to each NEO entry. Writes a new JSON with the same shape as input.

Usage examples:
  python3 -m scripts.neo_intercept_planner
  python3 -m scripts.neo_intercept_planner data/hazardous_neos/latest.json --output data/hazardous_neos/latest_intercept.json
  python3 -m scripts.neo_intercept_planner --r2-au 1.30 --output data/hazardous_neos/latest_intercept.json
"""

import os
import json
import math
import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
from datetime import datetime, timezone, timedelta

# --- Constants (SI) ---
AU = 1.495978707e11
MU_SUN = 1.32712440018e20
MU_EARTH = 3.986004418e14
R_EARTH = 6378.1363e3
G = 6.67430e-11
g0 = 9.80665

# --- Craft profile presets ---
PROFILES = {
    "neosurveyor": {
        "m0_kg": 1300.0,
        "Isp_s": 230.0,
        "h_leo_km": 500.0,
        "note": "NEO Surveyor–like mass; Isp assumed hydrazine monoprop (~230 s)."
    }
}

@dataclass
class MissionInputs:
    r2_AU: float = 1.0
    h_leo_km: float = 500.0
    Isp_s: float = 230.0
    m0_kg: float = 1300.0
    density_kg_m3: float = 2200.0
    rendezvous: bool = True
    aerocapture: bool = False
    profile: Optional[str] = "neosurveyor"
    profile_note: Optional[str] = PROFILES["neosurveyor"]["note"]

@dataclass
class HohmannResult:
    dv_helio_dep: float
    dv_helio_arr: float
    tof_days: float
    phase_angle_deg: float
    v_inf_dep: float
    dv_leo_escape: float
    dv_total_to_rendezvous: float

@dataclass
class AsteroidScenario:
    diameter_m: float
    radius_m: float
    density_kg_m3: float
    mass_kg: float
    mu_m3s2: float
    hill_m: float
    surface_g_m_s2: float
    chosen_orbit_radius_m: float
    altitude_over_surface_m: float
    v_circ_m_s: float
    period_hours: float
    dv_insertion_m_s: float
    dv_total_including_insert_m_s: float
    m_final_after_all_burns_kg: float
    propellant_used_kg: float

# --- Orbital math ---
def hohmann_1AU_to_r2(r2_AU: float) -> Tuple[float, float, float, float, float, float, float]:
    r1 = AU
    r2 = r2_AU * AU
    v1 = math.sqrt(MU_SUN / r1)
    v2 = math.sqrt(MU_SUN / r2)
    a_t = 0.5 * (r1 + r2)
    v_t1 = math.sqrt(MU_SUN * (2/r1 - 1/a_t))
    v_t2 = math.sqrt(MU_SUN * (2/r2 - 1/a_t))
    dv_dep = abs(v_t1 - v1)
    dv_arr = abs(v2 - v_t2)
    tof = math.pi * math.sqrt(a_t**3 / MU_SUN)
    if r2 > r1:
        phi = math.pi * (1 - math.sqrt((2*r2)/(r1 + r2)))
    else:
        phi = math.pi * (math.sqrt((2*r1)/(r1 + r2)) - 1)
    return dv_dep, dv_arr, tof, math.degrees(phi), v1, v2, v_t1

def leo_escape_dv_from_vinf(h_leo_km: float, v_inf: float) -> float:
    r_leo = R_EARTH + h_leo_km*1e3
    v_circ = math.sqrt(MU_EARTH / r_leo)
    v_esc = math.sqrt(2) * v_circ
    dv = math.sqrt(v_inf**2 + v_esc**2) - v_circ
    return dv

def mass_after_burn(m0: float, dv: float, Isp: float) -> float:
    if dv <= 0:
        return m0
    return m0 * math.exp(-dv/(g0*Isp))

def asteroid_params(diameter_m: float, density_kg_m3: float, a_AU: float) -> Tuple[float, float, float, float, float]:
    R = 0.5 * diameter_m
    volume = (4.0/3.0) * math.pi * R**3
    mass = density_kg_m3 * volume
    mu = G * mass
    a = a_AU * AU
    r_hill = a * (mu/(3.0*MU_SUN))**(1.0/3.0)
    g_surface = mu / (R**2) if R > 0 else 0.0
    return R, mass, mu, r_hill, g_surface

def choose_stable_orbit_radius(R_body: float, r_hill: float) -> float:
    # Heuristic: >= 3R and >= 1 km; < ~30% Hill
    r_lower = max(3.0*R_body, 1000.0)
    r_upper = 0.30 * r_hill
    if r_upper < r_lower:
        return r_lower
    return r_lower

def circular_orbit_stats(mu: float, r_orbit: float) -> Tuple[float, float]:
    if mu <= 0 or r_orbit <= 0:
        return 0.0, 0.0
    v = math.sqrt(mu / r_orbit)
    T = 2.0*math.pi * math.sqrt(r_orbit**3 / mu)
    return v, T

def plan_chain(params: MissionInputs) -> HohmannResult:
    dv_dep_helio, dv_arr_helio, tof, phi_deg, v1, v2, v_t1 = hohmann_1AU_to_r2(params.r2_AU)
    v_inf_dep = abs(v_t1 - v1)
    dv_leo = leo_escape_dv_from_vinf(params.h_leo_km, v_inf_dep)
    dv_arr_final = 0.0 if (params.aerocapture or not params.rendezvous) else dv_arr_helio
    dv_total = dv_leo + dv_dep_helio + dv_arr_final
    return HohmannResult(
        dv_helio_dep=dv_dep_helio,
        dv_helio_arr=dv_arr_final,
        tof_days=tof/86400.0,
        phase_angle_deg=phi_deg,
        v_inf_dep=v_inf_dep,
        dv_leo_escape=dv_leo,
        dv_total_to_rendezvous=dv_total,
    )

def scenario_for_diameter(d_m: float, params: MissionInputs, chain: HohmannResult) -> AsteroidScenario:
    R, m_ast, mu_ast, hill, g_s = asteroid_params(d_m, params.density_kg_m3, params.r2_AU)
    r_orbit = choose_stable_orbit_radius(R, hill)
    v_circ, T = circular_orbit_stats(mu_ast, r_orbit)
    dv_insert = v_circ if params.rendezvous and not params.aerocapture else 0.0

    # Mass flow across burns
    m1 = mass_after_burn(params.m0_kg, chain.dv_leo_escape, params.Isp_s)
    m2 = mass_after_burn(m1, chain.dv_helio_dep, params.Isp_s)
    m3 = mass_after_burn(m2, chain.dv_helio_arr, params.Isp_s)
    m_final = mass_after_burn(m3, dv_insert, params.Isp_s)

    return AsteroidScenario(
        diameter_m=d_m,
        radius_m=R,
        density_kg_m3=params.density_kg_m3,
        mass_kg=m_ast,
        mu_m3s2=mu_ast,
        hill_m=hill,
        surface_g_m_s2=g_s,
        chosen_orbit_radius_m=r_orbit,
        altitude_over_surface_m=r_orbit - R,
        v_circ_m_s=v_circ,
        period_hours=(T/3600.0 if T > 0 else 0.0),
        dv_insertion_m_s=dv_insert,
        dv_total_including_insert_m_s=chain.dv_total_to_rendezvous + dv_insert,
        m_final_after_all_burns_kg=m_final,
        propellant_used_kg=params.m0_kg - m_final
    )

# --- Time & I/O helpers ---
def load_neo_items(path: Path) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    p = Path(path)
    if not p.exists():
        sys.stderr.write(f"ERROR: Input file not found: {p}\n"
                         f"Tip: use an absolute path, or cd into the folder with the JSON.\n")
        sys.exit(2)
    try:
        data = json.loads(p.read_text())
    except json.JSONDecodeError as e:
        sys.stderr.write(f"ERROR: Failed to parse JSON in {p}:\n{e}\n")
        sys.exit(3)
    items = data.get("potentially_hazardous_neos") or data.get("neos") or []
    if isinstance(items, dict):
        items = [items]
    return data, items

def parse_close_approach_utc(item: Dict[str, Any]) -> Optional[datetime]:
    ca = item.get("close_approach") or {}
    s = ca.get("date_full") or ca.get("date") or None
    if not s:
        return None
    for fmt in ("%Y-%b-%d %H:%M", "%Y-%m-%d", "%Y-%m-%d %H:%M", "%Y-%b-%d %H:%M:%S"):
        try:
            return datetime.strptime(s, fmt).replace(tzinfo=timezone.utc)
        except Exception:
            pass
    return None

def synodic_days(r2_AU: float) -> float:
    """Approx synodic period (days) for Earth (1 AU) vs a circular target at r2_AU."""
    T_earth = 365.25
    T2 = 365.25 * (r2_AU ** 1.5)
    denom = abs((1.0/T_earth) - (1.0/T2))
    return 99999.0 if denom == 0 else (1.0 / denom)

def roll_forward_arrival(arrive_dt: datetime, now_dt: datetime, r2_AU: float) -> Tuple[datetime, int, float]:
    Tsyn = synodic_days(r2_AU)
    delta_days = (now_dt - arrive_dt).total_seconds() / 86400.0
    k = max(1, math.ceil(delta_days / Tsyn))
    return arrive_dt + timedelta(days=k*Tsyn), k, Tsyn

# --- Main ---
def main():
    ap = argparse.ArgumentParser(description="NEO Intercept Planner → JSON (Hohmann + size-based stable orbit)")
    ap.add_argument("input", nargs="?", help="Path to NEO JSON (defaults to $OUT_FILE or data/hazardous_neos/latest.json)")
    ap.add_argument("--output", type=str, help="Output JSON path (default: <input>_intercept.json)")
    ap.add_argument("--profile", type=str, choices=["neosurveyor", "custom"], default="neosurveyor",
                    help="Spacecraft profile (default: neosurveyor)")
    ap.add_argument("--r2-au", type=float, default=1.0, help="Target semi-major axis (AU). Use asteroid 'a' for accuracy.")
    ap.add_argument("--leo-alt", type=float, default=None, help="LEO altitude (km). Overrides profile default if set.")
    ap.add_argument("--isp", type=float, default=None, help="Specific impulse (s). Overrides profile default if set.")
    ap.add_argument("--m0", type=float, default=None, help="Initial mass in LEO (kg). Overrides profile default if set.")
    ap.add_argument("--density", type=float, default=2200.0, help="Asteroid bulk density (kg/m^3).")
    ap.add_argument("--no-rendezvous", action="store_true", help="Do not match asteroid heliocentric velocity (flyby).")
    ap.add_argument("--aerocapture", action="store_true", help="Skip arrival Δv (placeholder).")
    ap.add_argument("--arrive-utc", type=str, default=None,
                    help="Arrival UTC 'YYYY-mm-dd[ HH:MM]'. Overrides close_approach.date_full.")
    ap.add_argument("--window-policy", choices=["roll", "error", "ignore"], default="roll",
                    help="If arrival is in the past: roll (add synodic periods), error, or ignore. Default: roll.")
    args = ap.parse_args()

    # Resolve input/output paths (defaults to your snapshot and a separate intercept file)
    default_in = os.environ.get("OUT_FILE", "data/hazardous_neos/latest.json")
    input_path = Path(args.input or default_in)
    output_path = Path(args.output) if args.output else input_path.with_name(input_path.stem + "_intercept.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Input JSON
    base_data, items = load_neo_items(input_path)
    if not items:
        # Write a passthrough file with a generation timestamp so CI still commits an artifact
        out_data = dict(base_data)
        out_data["solutions_generated_utc"] = datetime.now(timezone.utc).isoformat()
        output_path.write_text(json.dumps(out_data, indent=2))
        print(f"No NEOs found; wrote passthrough file: {output_path}")
        return

    # Build mission inputs (profile first)
    if args.profile == "neosurveyor":
        preset = PROFILES["neosurveyor"]
        h_leo_km = preset["h_leo_km"] if args.leo_alt is None else args.leo_alt
        isp_s = preset["Isp_s"] if args.isp is None else args.isp
        m0_kg = preset["m0_kg"] if args.m0 is None else args.m0
        prof_note = preset["note"]
    else:
        h_leo_km = 500.0 if args.leo_alt is None else args.leo_alt
        isp_s = 230.0 if args.isp is None else args.isp
        m0_kg = 1300.0 if args.m0 is None else args.m0
        prof_note = "Custom profile (defaults align with NEOSurveyor-like unless overridden)."

    params = MissionInputs(
        r2_AU=args.r2_au,
        h_leo_km=h_leo_km,
        Isp_s=isp_s,
        m0_kg=m0_kg,
        density_kg_m3=args.density,
        rendezvous=(not args.no_rendezvous),
        aerocapture=args.aerocapture,
        profile=args.profile,
        profile_note=prof_note
    )

    # Arrival time (user override or from the JSON)
    arrive_dt_user: Optional[datetime] = None
    if args.arrive_utc:
        for fmt in ("%Y-%m-%d %H:%M", "%Y-%m-%d"):
            try:
                arrive_dt_user = datetime.strptime(args.arrive_utc, fmt).replace(tzinfo=timezone.utc)
                break
            except ValueError:
                continue
        if arrive_dt_user is None:
            print("WARNING: --arrive-utc not parseable; ignoring.", file=sys.stderr)

    # Shared Hohmann chain (depends on r2_AU & parking orbit)
    chain = plan_chain(params)

    # Build output by cloning input and annotating each NEO with intercept_plan
    out_data = dict(base_data)
    key = "potentially_hazardous_neos" if "potentially_hazardous_neos" in base_data else "neos"
    annotated_items = []
    now_utc = datetime.now(timezone.utc)

    for it in items:
        neo = dict(it)

        # Arrival/departure timestamps (if known)
        timing_meta: Dict[str, Any] = {}
        arrive_dt = arrive_dt_user or parse_close_approach_utc(neo)
        if arrive_dt and arrive_dt <= now_utc:
            if args.window_policy == "error":
                raise SystemExit(f"Arrival {arrive_dt.isoformat()} is in the past. "
                                 f"Use --arrive-utc (future) or --window-policy roll.")
            elif args.window_policy == "roll":
                new_arrive, k, Tsyn = roll_forward_arrival(arrive_dt, now_utc, params.r2_AU)
                timing_meta = {"rolled_forward_windows": k, "synodic_days": Tsyn}
                arrive_dt = new_arrive

        depart_dt = (arrive_dt - timedelta(days=chain.tof_days)) if arrive_dt else None

        # Diameter candidates
        dmin = neo.get("estimated_diameter_m_min") \
               or neo.get("estimated_diameter", {}).get("meters", {}).get("estimated_diameter_min")
        dmax = neo.get("estimated_diameter_m_max") \
               or neo.get("estimated_diameter", {}).get("meters", {}).get("estimated_diameter_max")

        scenarios_block: Dict[str, Any] = {}
        if dmin is not None:
            scn_min = scenario_for_diameter(float(dmin), params, chain)
            scenarios_block["min_diameter"] = {
                "diameter_m": scn_min.diameter_m,
                "asteroid_mu_m3s2": scn_min.mu_m3s2,
                "asteroid_hill_radius_m": scn_min.hill_m,
                "surface_g_m_s2": scn_min.surface_g_m_s2,
                "chosen_orbit_radius_m": scn_min.chosen_orbit_radius_m,
                "altitude_over_surface_m": scn_min.altitude_over_surface_m,
                "v_circ_m_s": scn_min.v_circ_m_s,
                "period_hours": scn_min.period_hours,
                "dv_insertion_m_s": scn_min.dv_insertion_m_s,
                "dv_total_including_insert_m_s": scn_min.dv_total_including_insert_m_s,
                "final_mass_kg": scn_min.m_final_after_all_burns_kg,
                "propellant_used_kg": scn_min.propellant_used_kg
            }
        if dmax is not None and (dmin is None or float(dmax) != float(dmin)):
            scn_max = scenario_for_diameter(float(dmax), params, chain)
            scenarios_block["max_diameter"] = {
                "diameter_m": scn_max.diameter_m,
                "asteroid_mu_m3s2": scn_max.mu_m3s2,
                "asteroid_hill_radius_m": scn_max.hill_m,
                "surface_g_m_s2": scn_max.surface_g_m_s2,
                "chosen_orbit_radius_m": scn_max.chosen_orbit_radius_m,
                "altitude_over_surface_m": scn_max.altitude_over_surface_m,
                "v_circ_m_s": scn_max.v_circ_m_s,
                "period_hours": scn_max.period_hours,
                "dv_insertion_m_s": scn_max.dv_insertion_m_s,
                "dv_total_including_insert_m_s": scn_max.dv_total_including_insert_m_s,
                "final_mass_kg": scn_max.m_final_after_all_burns_kg,
                "propellant_used_kg": scn_max.propellant_used_kg
            }

        neo["intercept_plan"] = {
            "profile": params.profile,
            "profile_note": params.profile_note,
            "assumptions": {
                "r2_AU": params.r2_AU,
                "leo_alt_km": params.h_leo_km,
                "Isp_s": params.Isp_s,
                "m0_kg": params.m0_kg,
                "density_kg_m3": params.density_kg_m3,
                "rendezvous": params.rendezvous,
                "aerocapture": params.aerocapture
            },
            "timing": {
                "arrive_utc": (arrive_dt.isoformat() if arrive_dt else None),
                "depart_utc": (depart_dt.isoformat() if depart_dt else None),
                "transfer_time_days": chain.tof_days,
                "phase_angle_deg": chain.phase_angle_deg,
                "meta": timing_meta or None
            },
            "heliocentric_chain": {
                "dv_helio_dep_m_s": chain.dv_helio_dep,
                "dv_helio_arr_m_s": chain.dv_helio_arr,
                "dv_leo_escape_m_s": chain.dv_leo_escape,
                "dv_total_to_rendezvous_m_s": chain.dv_total_to_rendezvous,
                "v_inf_dep_m_s": chain.v_inf_dep
            },
            "stable_orbit_scenarios": scenarios_block
        }

        annotated_items.append(neo)

    out_data[key] = annotated_items
    out_data["solutions_generated_utc"] = datetime.now(timezone.utc).isoformat()
    output_path.write_text(json.dumps(out_data, indent=2))

    print(f"Wrote: {output_path}")
    print(f"Assumed profile: {params.profile} | r2_AU={params.r2_AU} | LEO={params.h_leo_km} km | Isp={params.Isp_s} s | m0={params.m0_kg} kg")
    print(f"Chain: Δv_dep={chain.dv_helio_dep:.1f} m/s, Δv_arr={chain.dv_helio_arr:.1f} m/s, "
          f"LEO→esc={chain.dv_leo_escape:.1f} m/s, v∞={chain.v_inf_dep:.1f} m/s, TOF={chain.tof_days:.2f} d")

if __name__ == "__main__":
    main()
