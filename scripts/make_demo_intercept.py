# scripts/make_demo_intercept.py
#!/usr/bin/env python3
from __future__ import annotations
import json, math, time
from pathlib import Path
import numpy as np

AU_M = 149_597_870_700.0
MU_SUN = 1.32712440018e20
SID_Y = 365.256363004
DEG = math.pi/180.0

def jd_now():
    return time.time()/86400.0 + 2440587.5

def earth_r_au(jd):
    th = 2*math.pi*((jd-2451545.0)/SID_Y)
    return np.array([math.cos(th), math.sin(th), 0.0], float)

def kepler_E(M,e,tol=1e-12,maxit=50):
    M=(M+math.pi)%(2*math.pi)-math.pi
    E=M if e<0.8 else math.pi
    for _ in range(maxit):
        f=E-e*math.sin(E)-M; fp=1-e*math.cos(E); d=-f/fp; E+=d
        if abs(d)<tol: break
    return E

def target_r_au(a_AU,e,i_deg,raan_deg,argp_deg,ma0_deg,epoch0_jd,t_jd):
    a=a_AU*AU_M; i=i_deg*DEG; Om=raan_deg*DEG; w=argp_deg*DEG; M0=ma0_deg*DEG
    n=math.sqrt(MU_SUN/a**3); dt=(t_jd-epoch0_jd)*86400.0; M=M0+n*dt
    E=kepler_E(M,e)
    x_p=a*(math.cos(E)-e); y_p=a*(math.sqrt(1-e*e)*math.sin(E))
    cO,sO,ci,si,cw,sw = math.cos(Om),math.sin(Om),math.cos(i),math.sin(i),math.cos(w),math.sin(w)
    R=np.array([[cO*cw - sO*sw*ci, -cO*sw - sO*cw*ci,  sO*si],
                [sO*cw + cO*sw*ci, -sO*sw + cO*cw*ci, -cO*si],
                [           sw*si,             cw*si,      ci]], float)
    r = (R @ np.array([x_p,y_p,0.0]))/AU_M
    return r

def demo_json(n=500, depart_in_days=40.0, tof_days=220.0):
    # Target elements
    elems = {"a_AU": 1.6, "e": 0.22, "i_deg": 15.0, "raan_deg": 120.0, "argp_deg": 40.0}
    elems["epoch_jd"] = jd_now()
    elems["ma_deg"] = 10.0  # mean anomaly at epoch

    depart_jd = elems["epoch_jd"] + depart_in_days
    arrive_jd = depart_jd + tof_days

    r1 = earth_r_au(depart_jd)
    r2 = target_r_au(elems["a_AU"], elems["e"], elems["i_deg"], elems["raan_deg"], elems["argp_deg"],
                     elems["ma_deg"], elems["epoch_jd"], arrive_jd)

    # Smooth curved path from r1 to r2 (bezier-ish; just for demo)
    t  = np.linspace(0, 1, n)[:,None]
    ctrl = (r1*1.10 + r2*0.30) + np.array([[0.0, 0.0, 0.12]])  # gentle bump
    p   = (1-t)*((1-t)*r1 + t*ctrl) + t*((1-t)*ctrl + t*r2)

    plan = {
        "elements": elems,
        "depart_epoch_jd": float(depart_jd),
        "arrive_epoch_jd": float(arrive_jd),
        "tof_days": float(tof_days),
        "r1_au": r1.tolist(),
        "r2_au": r2.tolist(),
        "lambert_polyline_xyz_au": p.reshape(-1,3).tolist()
    }
    return {
      "date_utc": "demo",
      "snapshot_utc": "demo",
      "count": 1,
      "potentially_hazardous_neos": [
        { "name": "DEMO", "neo_reference_id": "0", "intercept_plan": plan }
      ]
    }

if __name__ == "__main__":
    out = Path("data/hazardous_neos/latest_intercept_poly.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(demo_json(), indent=2))
    print("Wrote", out)
