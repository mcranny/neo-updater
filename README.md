# NEO Updater & Intercept Planner

A compact pipeline to:
- Pull **Potentially Hazardous Asteroids (PHAs)** from NASA **NeoWs**.
- Write a daily snapshot → `data/hazardous_neos/latest.json`.
- Plan Earth→NEO transfers:
  - **Baseline**: circular/coplanar **Hohmann** sizing (default in CI).
  - **High-fidelity (opt-in)**: **SBDB elements + Lambert** rendezvous.
- Visualize trajectories locally with a Qt viewer (CI skips Qt).

For full derivations and algorithms, see **[docs/math.md](docs/math.md)**.

---

## Architecture

    NeoWs (daily) ─┬─▶ data/hazardous_neos/latest.json
                   │
                   ├─▶ Hohmann (default CI) ─▶ latest_intercept.json
                   │
                   └─▶ +SBDB +Lambert (opt-in) ─▶ latest_intercept.json (adds lambert/*)

- **Updater** — `scripts/update_hazardous_neos.py`  
  Uses `NASA_API_KEY` (single-day window; `count=0` is normal).
- **Planner** — `scripts/neo_intercept_planner.py`  
  Parses NeoWs feed/browse/flat → emits per-NEO `intercept_plan`.  
  Adds `elements/*`, `lambert/*`, and a `lambert_polyline_xy_au` when enabled.
- **Viewer** — `app/neo_viewer_qt.py`  
  Hohmann mode (half-ellipse + phasing) or Lambert mode (polyline; SC rides the arc).

---

## Install (condensed)

- Python **3.10+**
- Virtual environment recommended
- Deps: `requests`, `python-dateutil`, `python-dotenv`, `numpy`, `PySide6` (viewer only)

~~~bash
python -m venv .venv && source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
~~~

Optional `.env` for the updater:
~~~text
NASA_API_KEY=YOUR_REAL_KEY
OUT_FILE=data/hazardous_neos/latest.json
~~~

---

## Usage

### Baseline Hohmann (default / CI)
~~~bash
python -m scripts.update_hazardous_neos

python -m scripts.neo_intercept_planner \
  --input  data/hazardous_neos/latest.json \
  --output data/hazardous_neos/latest_intercept.json \
  --roll-past

python app/neo_viewer_qt.py
~~~

### Elements + Lambert (true rendezvous)
~~~bash
python -m scripts.neo_intercept_planner \
  --input  data/hazardous_neos/latest.json \
  --output data/hazardous_neos/latest_intercept.json \
  --roll-past \
  --elements \
  --lambert \
  --depart-utc 2026-01-05T00:00:00Z \
  --tof-days 180

python app/neo_viewer_qt.py data/hazardous_neos/latest_intercept.json
~~~

**Selected flags**
~~~text
--r1-au <float>      Departure heliocentric radius (AU), default 1.0
--r2-au <float>      Override target radius for all NEOs (AU)
--leo-km <float>     Parking altitude for LEO escape model (km)
--roll-past          If arrival is in the past, add synodic periods until future
--elements           Attach SBDB elements (a,e,i,Ω,ω,M at epoch)
--lambert            Single-rev, prograde Lambert between r1(t1) and r2(t2)
--depart-utc <ISO>   Departure epoch for Lambert
--tof-days <float>   Time of flight for Lambert
~~~

---

## Math (essentials)

### Units & frame
- Canonical solar units for reasoning: AU, year, with $\mu_\odot = 4\pi^2$. Then $P^2 = a^3$, $n = 2\pi/a^{3/2}$.
- SI for $\Delta v$ and LEO injection. Viewer draws the ecliptic-plane XY projection.

### Hohmann (circular→circular, coplanar)
Given heliocentric radii $r_1, r_2$ (AU), transfer ellipse:

$$a_t=\tfrac{1}{2}(r_1+r_2),\quad b_t=\sqrt{r_1 r_2},\quad e_t=\sqrt{1-(b_t/a_t)^2}$$

Speeds:

$$v_1=\sqrt{\mu/r_1},\quad v_2=\sqrt{\mu/r_2},\quad v_p=\sqrt{\mu(2/r_1-1/a_t)},\quad v_a=\sqrt{\mu(2/r_2-1/a_t)}$$

Δv (Sun-centric):  

$$\Delta v_1=\lvert v_p-v_1\rvert,\quad \Delta v_2=\lvert v_2-v_a\rvert$$  

Time of flight:  

$$\mathrm{TOF}=\pi\sqrt{a_t^3/\mu}$$

**LEO escape (patched conic)** with LEO radius $r_L$, $v_L=\sqrt{\mu_\oplus/r_L}$:

$$\Delta v_{\text{LEO}}=\sqrt{v_\infty^2+(\sqrt{2}\,v_L)^2}-v_L,\qquad v_\infty \approx \Delta v_1$$

**Phasing (viewer legacy):** pick target start angle so it reaches the end of the half-ellipse at $t=\mathrm{TOF}$.
Outward example:  

$$\theta_T(0)=\pi - n_2 \cdot \mathrm{TOF}\pmod{2\pi}$$

## Lambert rendezvous (opt-in)

Propagate Earth $\mathbf{r}_\oplus(t_1)$ and NEO $\mathbf{r}_A(t_2)$ from elements.

Solve universal-variables **Lambert** (prograde, single-rev) to get $\mathbf{v}_1,\ \mathbf{v}_2$.

Report:

$$
\Delta v_{\mathrm{depart}}=\left\|\mathbf{v}_1-\mathbf{v}_\oplus(t_1)\right\|,\quad
\Delta v_{\mathrm{arrive}}=\left\|\mathbf{v}_2-\mathbf{v}_A(t_2)\right\|
$$

Optionally expose $C_3=v_\infty^2$ (future).

The viewer draws the Lambert polyline (XY in AU); the SC rides the arc and intercepts at arrival.

For deeper coverage (Kepler solvers, Stumpff functions, residual checks), see **[docs/math.md](docs/math.md)**.

---

## Output schema (subset)

Per NEO (always present):
~~~json
{
  "intercept_plan": {
    "schema_version": "1.3.0",
    "r1_AU": 1.0,
    "r2_AU": 1.23,
    "tof_days": 184.2,
    "synodic_days": 398.9,
    "departure_utc": "YYYY-MM-DDTHH:MM:SSZ",
    "arrival_utc": "YYYY-MM-DDTHH:MM:SSZ",
    "rolled_forward": true,
    "dv_depart_heliocentric_m_s": 3215.4,
    "dv_arrive_heliocentric_m_s": 1712.8,
    "dv_from_LEO_m_s": 3940.1,
    "dv_total_m_s": 5652.9
  }
}
~~~

When `--elements`:
~~~json
{
  "elements": {
    "a_AU": 1.234,
    "e": 0.12,
    "i_deg": 5.8,
    "raan_deg": 87.3,
    "argp_deg": 210.4,
    "M_deg": 14.2,
    "epoch_utc": "YYYY-MM-DDTHH:MM:SSZ",
    "source": "JPL SBDB"
  }
}
~~~

When `--lambert`:
~~~json
{
  "lambert": {
    "t_depart_utc": "YYYY-MM-DDTHH:MM:SSZ",
    "t_arrive_utc": "YYYY-MM-DDTHH:MM:SSZ",
    "tof_days": 180.0,
    "dv_depart_m_s": 3120.2,
    "dv_arrive_m_s": 1650.8
  },
  "lambert_polyline_xy_au": [[x_au, y_au], ...]
}
~~~

---

## CI

- Workflow: `.github/workflows/pha-4x-daily.yml`
  - Scheduled at **00/06/12/18 UTC** + `workflow_dispatch`.
  - Default job = Hohmann (back-compat).
  - Enable **elements + Lambert** via dispatch inputs or set `env.LAMBERT="true"`.

---

## Limits

- Two-body dynamics; no third-body/J2/finite-burn effects.
- Lambert search: prograde, single-rev (no multi-rev/window scan yet).
- Earth departure uses circular 1 AU for sizing; launcher C3 mapping is approximate.
- Arrival models heliocentric matching; micro-SOI capture not included.

---

## Validation

- Kepler \(E(M,e)\) and frame rotations (PQW→IJK).
- Universal-variable propagation residuals.
- Lambert residual \( \|\mathbf{r}(t_1+\Delta t)-\mathbf{r}_2\| \).
- JSON schema snapshots.

---

## References

Prussing & Conway; Battin; Vallado; Izzo (“Revisiting Lambert’s Problem”).
