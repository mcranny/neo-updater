# scripts/plan_intercepts.py
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from scripts.ephem_earth import earth_rv
from scripts.orbital import (
    AU_M,
    MU_SUN,
    deg2rad,
    kepler_E_from_M,
    lambert_universal,  # r1_m, r2_m, tof_s, mu -> (v1_mps, v2_mps)
    sample_transfer_polyline,  # r1_m, v1_mps, tof_s, N -> list[[xAU,yAU],...]
    wrap_2pi,
)

DAY_S = 86400.0


def _asteroid_M_at_time(el: dict[str, Any], t_jd: float) -> float | None:
    a = float(el["a_AU"]) * AU_M
    n = math.sqrt(MU_SUN / (a**3))  # rad/s
    tp = el.get("tp_jd_tdb")
    if tp is not None:
        return wrap_2pi(n * ((t_jd - float(tp)) * DAY_S))
    ma_deg = el.get("ma_deg")
    epoch = el.get("epoch_jd_tdb")
    if ma_deg is None or epoch is None:
        return None
    M0 = deg2rad(float(ma_deg))
    dt = (t_jd - float(epoch)) * DAY_S
    return wrap_2pi(M0 + n * dt)


def asteroid_rv_from_elements(el: dict[str, Any], jd_tdb: float) -> tuple[np.ndarray, np.ndarray]:
    a = float(el["a_AU"]) * AU_M
    e = float(el["e"])
    inc = deg2rad(float(el["i_deg"]))
    raan = deg2rad(float(el["raan_deg"]))
    argp = deg2rad(float(el["argp_deg"]))

    M = _asteroid_M_at_time(el, jd_tdb)
    if M is None:
        raise ValueError("insufficient elements to compute mean anomaly (need tp or ma@epoch)")

    E = kepler_E_from_M(M, e)
    cE, sE = math.cos(E), math.sin(E)
    r_peri = np.array([a * (cE - e), a * (math.sqrt(1.0 - e * e) * sE), 0.0])
    r_mag = a * (1.0 - e * cE)
    fac = math.sqrt(MU_SUN * a) / r_mag
    v_peri = np.array([-fac * sE, fac * math.sqrt(1.0 - e * e) * cE, 0.0])

    cO, sO = math.cos(raan), math.sin(raan)
    ci, si = math.cos(inc), math.sin(inc)
    cw, sw = math.cos(argp), math.sin(argp)
    R = np.array(
        [
            [cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
            [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
            [sw * si, cw * si, ci],
        ]
    )
    return R @ r_peri, R @ v_peri


def leo_departure_dv(vinf_kms: float, leo_alt_m: float = 500e3) -> float:
    MU_EARTH = 3.986004418e14
    R_EARTH = 6378137.0
    r = R_EARTH + leo_alt_m
    v_circ = math.sqrt(MU_EARTH / r)
    v_esc = math.sqrt(2 * MU_EARTH / r)
    v_hyp = math.sqrt((vinf_kms * 1000.0) ** 2 + v_esc * v_esc)
    return (v_hyp - v_circ) / 1000.0


def plan_intercept_for_object(
    obj: dict[str, Any],
    tof_days_grid: list[int] | None = None,
    arrive_offset_hours: list[int] | None = None,
    leo_alt_m: float = 500e3,
) -> dict[str, Any] | None:
    if tof_days_grid is None:
        tof_days_grid = list(range(30, 181, 10))
    if arrive_offset_hours is None:
        arrive_offset_hours = [-12, -6, 0, 6, 12]
    el = obj.get("elements")
    ca = obj.get("next_ca")
    if not el or not ca:
        return None

    jd_ca = float(ca["jd_tdb"])
    best: dict[str, Any] | None = None

    for d_off in arrive_offset_hours:
        jd_arr = jd_ca + d_off / 24.0
        r_ast_arr, v_ast_arr = asteroid_rv_from_elements(el, jd_arr)

        for tof_d in tof_days_grid:
            jd_dep = jd_arr - tof_d
            r_e_dep, v_e_dep = earth_rv(jd_dep)
            tof_s = tof_d * DAY_S

            solution = None
            for prograde in (True, False):
                try:
                    solution = lambert_universal(
                        r_e_dep, r_ast_arr, tof_s, MU_SUN, prograde=prograde
                    )
                    break
                except (ValueError, RuntimeError, OverflowError, ArithmeticError):
                    continue
            if solution is None:
                continue
            v1, v2 = solution

            dv_dep = float(np.linalg.norm(v1 - v_e_dep))
            dv_arr = float(np.linalg.norm(v_ast_arr - v2))
            dv_sum = dv_dep + dv_arr

            if best is None or dv_sum < best["dv_sum_mps"]:
                try:
                    poly = sample_transfer_polyline(r_e_dep, v1, tof_s, 180)
                except Exception:
                    poly = []
                best = {
                    "depart_jd_tdb": jd_dep,
                    "arrive_jd_tdb": jd_arr,
                    "tof_days": tof_d,
                    "dv_depart_kms": dv_dep / 1000.0,
                    "dv_arrive_kms": dv_arr / 1000.0,
                    "dv_sum_mps": dv_sum,
                    "c3_km2_s2": (dv_dep / 1000.0) ** 2,
                    "vinf_kms": dv_dep / 1000.0,
                    "leo_dv_kms": leo_departure_dv(dv_dep / 1000.0, leo_alt_m=leo_alt_m),
                    "lambert_polyline_xy_au": poly,
                }
    return best


def attach_intercepts(
    payload: dict[str, Any],
    tof_days_grid: list[int] | None = None,
    arrive_offset_hours: list[int] | None = None,
    leo_alt_m: float = 500e3,
) -> dict[str, Any]:
    if tof_days_grid is None:
        tof_days_grid = list(range(30, 181, 10))
    if arrive_offset_hours is None:
        arrive_offset_hours = [-12, -6, 0, 6, 12]

    for o in payload.get("objects", []):
        o["intercept"] = plan_intercept_for_object(
            o,
            tof_days_grid=tof_days_grid,
            arrive_offset_hours=arrive_offset_hours,
            leo_alt_m=leo_alt_m,
        )

    payload["intercept_note"] = {
        "frame": "heliocentric ecliptic (Sun μ)",
        "earth_model": "DE ephemeris via Skyfield (fallback: Kepler)",
        "cost": "minimize |v1−vE| + |vast−v2|",
        "dv_units": "km/s",
        "time_scale": "TDB JD",
    }
    return payload


def main():
    ap = argparse.ArgumentParser(description="Attach Lambert intercepts to CAD-first feed.")
    ap.add_argument("--inp", default="data/hazardous_neos/latest.json")
    ap.add_argument("--out", default="data/hazardous_neos/latest_intercepts.json")
    ap.add_argument("--leo-alt-m", type=float, default=500e3)
    ap.add_argument("--tof-min", type=int, default=30)
    ap.add_argument("--tof-max", type=int, default=180)
    ap.add_argument("--tof-step", type=int, default=10)
    ap.add_argument("--arrive-hours", default="-12,-6,0,6,12")
    args = ap.parse_args()

    src = json.loads(Path(args.inp).read_text())
    tof = list(range(args.tof_min, args.tof_max + 1, args.tof_step))
    offs = [int(s) for s in args.arrive_hours.split(",") if s]
    out = attach_intercepts(
        src, tof_days_grid=tof, arrive_offset_hours=offs, leo_alt_m=args.leo_alt_m
    )
    Path(args.out).write_text(json.dumps(out, indent=2))
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
