```markdown
# NEO Updater & Intercept Planner

Pipeline to:
- Pull **Potentially Hazardous Asteroids (PHAs)** from NASA **NeoWs**.
- Write a daily snapshot → `data/hazardous_neos/latest.json`.
- Plan Earth→NEO transfers:
  - **Baseline:** circular/coplanar **Hohmann** sizing.
  - **High-fidelity (opt-in):** **SBDB elements + Lambert** rendezvous.
- Visualize trajectories with a Qt viewer (local tool; not used in CI).

---

## Architecture

```

NeoWs (daily)  ──▶  data/hazardous\_neos/latest.json
│
├─▶  Hohmann (default CI) ──▶  latest\_intercept.json
│
└─▶  +SBDB +Lambert (opt-in) ─▶ latest\_intercept.json (adds lambert/\*)

````

- **Updater**: `scripts.update_hazardous_neos` (requires `NASA_API_KEY`).
- **Planner**: `scripts.neo_intercept_planner`  
  - Accepts NeoWs feed/browse or flattened lists.  
  - Emits per-NEO `intercept_plan` (Hohmann fields always; Lambert fields when enabled).
- **Viewer**: `app/neo_viewer_qt.py`  
  - Hohmann mode: half-ellipse + phasing.  
  - Lambert mode: draws solver polyline; SC rides the arc; intercept at arrival.

---

## Install (condensed)

- Python **3.10+**; `pip` in a venv recommended.
- Dependencies: `requests`, `python-dateutil`, `python-dotenv`, `numpy`, `PySide6` (viewer only; skipped in CI).

```bash
python -m venv .venv && source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
````

`.env` (for updater):

```
NASA_API_KEY=YOUR_REAL_KEY
OUT_FILE=data/hazardous_neos/latest.json
```

---

## Usage

### Hohmann baseline (default/CI)

```bash
python -m scripts.update_hazardous_neos
python -m scripts.neo_intercept_planner \
  --input  data/hazardous_neos/latest.json \
  --output data/hazardous_neos/latest_intercept.json \
  --roll-past
python app/neo_viewer_qt.py
```

### Elements + Lambert (true rendezvous)

```bash
python -m scripts.neo_intercept_planner \
  --input  data/hazardous_neos/latest.json \
  --output data/hazardous_neos/latest_intercept.json \
  --roll-past \
  --elements \
  --lambert \
  --depart-utc 2026-01-05T00:00:00Z \
  --tof-days 180
python app/neo_viewer_qt.py data/hazardous_neos/latest_intercept.json
```

> `--elements` pulls osculating elements from JPL **SBDB** (no key).
> `--lambert` solves $\mathbf{r}_1(t_1) \to \mathbf{r}_2(t_2)$ with specified $t_1$ and TOF.

---

## Planner CLI (selected)

```
--r1-au <float>           Departure heliocentric radius (AU), default 1.0.
--r2-au <float>           Override target radius for all NEOs (AU).
--leo-km <float>          Parking altitude (km) for LEO escape model.
--roll-past               If arrival is in the past, roll forward by synodic periods.
--elements                Attach SBDB elements (a,e,i,Ω,ω,M at epoch).
--lambert                 Compute a single prograde, single-rev Lambert.
--depart-utc <ISO>        Departure epoch for Lambert.
--tof-days <float>        Time of flight for Lambert.
```

---

## Math & Modeling

### Units & frame

* **Canonical solar units** for Hohmann/propagation reasoning: AU, year, $\mu_\odot = 4\pi^2$ → $P^2 = a^3$, $n = 2\pi/a^{3/2}$.
* SI for Δv and LEO injection. Viewer draws the **ecliptic XY** projection.

### Circular Hohmann (baseline sizing)

Given $r_1, r_2$ (AU):

* Transfer ellipse:
  $a_t = \tfrac{1}{2}(r_1 + r_2), \quad b_t = \sqrt{r_1 r_2}, \quad e_t = \sqrt{1 - (b_t/a_t)^2}$
* Speeds (Sun-centric):
  $v_1 = \sqrt{\mu/r_1}, \; v_2 = \sqrt{\mu/r_2}$
  On transfer: $v_p = \sqrt{\mu(2/r_1 - 1/a_t)}, \; v_a = \sqrt{\mu(2/r_2 - 1/a_t)}$
* Δv (Sun-centric): $\Delta v_1 = |v_p - v_1|, \; \Delta v_2 = |v_2 - v_a|$
* Time of flight: $\text{TOF} = \pi \sqrt{a_t^3/\mu}$

**Phasing** (viewer legacy mode): for outward transfers we offset the target start angle by
$\theta_T(0) = \pi - n_2 \cdot \text{TOF}$ so the target reaches the end of the half-ellipse at TOF.

**LEO escape** (patched conic): with LEO radius $r_L$,

$$
\Delta v_{\text{LEO}} = \sqrt{v_\infty^2 + \left(\sqrt{2}\,v_L\right)^2} - v_L, \quad v_L = \sqrt{\mu_\oplus / r_L}
$$

where $v_\infty \approx \Delta v_1$ (Sun-centric to Earth-relative mapping).

**Synodic roll-forward**: $P_{\text{syn}} = \frac{2\pi}{|n_1 - n_2|}$ to move past stale arrival epochs.

### Lambert rendezvous (opt-in)

We solve the two-point boundary value problem in **universal variables** with Stumpff functions $C(z), S(z)$:

* Inputs: $\mathbf{r}_1(t_1), \mathbf{r}_2(t_2), \Delta t$.
* Output: $\mathbf{v}_1, \mathbf{v}_2$ (short-way, **prograde**, single-rev by default).
* Δv accounting:
  $\Delta v_{\text{depart}} = \|\mathbf{v}_1 - \mathbf{v}_{\oplus}(t_1)\|, \;
     \Delta v_{\text{arrive}} = \|\mathbf{v}_2 - \mathbf{v}_{\text{NEO}}(t_2)\|$.
  Report **C3** via $C_3 = v_\infty^2$ if/when exposed.

**Propagation** uses universal-variable Kepler to sample the arc into a polyline (XY in AU) for the viewer.
**Elements** are ingested as $a,e,i,\Omega,\omega,M_0$ at epoch; target states are rotated 3-D vectors.

---

## JSON Contract

Baseline (always):

```jsonc
"intercept_plan": {
  "schema_version": "1.3.0",
  "r1_AU": 1.0,
  "r2_AU": 1.23,
  "tof_days": 184.2,
  "synodic_days": 398.9,
  "departure_utc": "YYYY-MM-DDTHH:MM:SSZ",
  "arrival_utc":   "YYYY-MM-DDTHH:MM:SSZ",
  "rolled_forward": true,
  "dv_depart_heliocentric_m_s": 3215.4,
  "dv_arrive_heliocentric_m_s": 1712.8,
  "dv_from_LEO_m_s": 3940.1,
  "dv_total_m_s": 5652.9,
  "leo_altitude_m": 500000,
  "spacecraft_mass_kg": 1300.0,
  "Isp_s": 230.0,
  "propellant_estimate_kg": 421.3
}
```

When `--elements`:

```jsonc
"elements": {
  "a_AU": 1.234, "e": 0.12,
  "i_deg": 5.8, "raan_deg": 87.3, "argp_deg": 210.4,
  "M_deg": 14.2, "epoch_utc": "YYYY-MM-DDTHH:MM:SSZ", "source": "JPL SBDB"
}
```

When `--lambert`:

```jsonc
"lambert": {
  "t_depart_utc": "YYYY-MM-DDTHH:MM:SSZ",
  "t_arrive_utc": "YYYY-MM-DDTHH:MM:SSZ",
  "tof_days": 180.0,
  "dv_depart_m_s": 3120.2,
  "dv_arrive_m_s": 1650.8
},
"lambert_polyline_xy_au": [[x,y], ...]   // XY ecliptic polyline for the viewer
```

---

## CI

* Workflow: `.github/workflows/pha-4x-daily.yml` scheduled at **00/06/12/18 UTC** plus `workflow_dispatch`.
* Default job runs **Hohmann** only (back-compat).
  Opt-in **elements + Lambert** via dispatch inputs or by setting `env.LAMBERT="true"` in the workflow.
* Viewer deps (Qt) are **skipped** on CI; it’s a local tool.

---

## Limitations

* Two-body dynamics; no J2/third-body/finite burns.
* **Lambert**: prograde, single-rev by default (no multi-rev search yet).
* Earth is modeled as circular 1 AU for departure $\mathbf{v}_\oplus(t_1)$; mapping to launch vehicle C3 is approximate.
* Arrival treats heliocentric matching; capture in the asteroid SOI is not modeled.

---

## Validation & Tests

* Unit tests for Kepler solver (E(M,e)) and frame rotations (PQW→IJK).
* Numerical checks for Lambert residuals $|\Delta t_{\text{prop}} - \Delta t|$.
* JSON schema snapshots for planner output.

---
