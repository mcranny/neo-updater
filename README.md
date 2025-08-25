
# NEO Updater & Intercept Planner

A small pipeline that:
- Fetches **Potentially Hazardous Asteroids (PHAs)** from NASA NeoWs.
- Writes a daily snapshot to `data/hazardous_neos/latest.json`.
- Plans a **Hohmann**-style Earth→NEO transfer for each object and writes `latest_intercept.json`.
- Includes a **Qt viewer** that animates Earth, target, and transfer with intercept/closest-approach detection.
- Runs **4×/day on GitHub Actions** and can be run locally.

**Repo:** `mcranny/neo-updater` ─ **Branch:** `main`

---

## Quick Start

### Requirements
- Python **3.10+** (macOS/Windows/Linux)
- `pip` (in a **virtual environment** recommended)
- NASA API key from https://api.nasa.gov (free)

### Install
```bash
git clone https://github.com/mcranny/neo-updater.git
cd neo-updater
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
python -m pip install -r requirements.txt   # requests, python-dateutil, PySide6
```

### Configure environment
Create `.env` (kept out of git by `.gitignore`):
```bash
cat > .env << 'EOF'
NASA_API_KEY=YOUR_REAL_KEY
OUT_FILE=data/hazardous_neos/latest.json
EOF
```

### Run locally
```bash
# 1) Update snapshot from NeoWs
python -m scripts.update_hazardous_neos

# 2) Plan Hohmann transfers (adds *_intercept.json)
python -m scripts.neo_intercept_planner --output data/hazardous_neos/latest_intercept.json --roll-past

# 3) Visualize
python app/neo_viewer_qt.py
```

**Outputs:**  
`data/hazardous_neos/latest.json` and `data/hazardous_neos/latest_intercept.json`

---

## GitHub Actions (CI)

Workflow: `.github/workflows/pha-4x-daily.yml`  
Triggers: `0 0,6,12,18 * * *` **UTC** (≈ **17:00, 23:00, 05:00, 11:00** US‑Pacific).

Setup:
1. Add repo **Secret** `NASA_API_KEY`.
2. Ensure workflow has `contents: write` permissions.
3. (Optional) Set author to your **noreply** email so commits count on your graph.

The job runs updater → planner, and commits changes if any.

---

## Minimal Math Overview (What the Viewer Shows)

We use a **coplanar, circular** Sun‑centered model in **canonical solar units**:

- Length **AU**, time **year**, Sun μ: `μ☉ = 4π² AU³/yr²` ⇒ `P²=a³`, `n=2π/a^(3/2)`.
- Earth radius: `r₁ = 1 AU`. Target radius from planner: `r₂` (AU).

**Hohmann transfer (circular→circular):**
- `a_t = (r₁ + r₂)/2`,  `b_t = √(r₁ r₂)`,  `e_t = √(1 - (b_t²/a_t²))`
- Speeds: `v₁=√(μ☉/r₁)`, `v₂=√(μ☉/r₂)`, `v_p=√( μ☉(2/r₁ - 1/a_t) )`, `v_a=√( μ☉(2/r₂ - 1/a_t) )`
- Δv (heliocentric): `Δv₁=|v_p - v₁|`, `Δv₂=|v₂ - v_a|`, `Δv_total≈Δv₁+Δv₂`
- **Time of Flight:** `TOF = π √(a_t³/μ☉)`

**Phasing (so target is there at `t=TOF`):**
- Choose target start angle `θ_T(0)` so that `θ_T(TOF)` aligns with transfer end. For outward transfers:  
  `θ_T(0) = π − n₂·TOF` (in radians, modulo `2π`).

**Spacecraft motion on the ellipse:**
- Mean anomaly: `M(t) = π·(t/TOF)`; solve `M = E − e_t sin E` by Newton’s method.
- Position: `x = a_t (cosE − e_t)`, `y = b_t sinE` (Sun at a focus, perihelion on +x).

**LEO injection (patched conics):**
- Hyperbolic excess `v∞ = Δv₁` (convert AU/yr→m/s).
- For a circular LEO of radius `r_LEO = (R⊕ + 500 km)`:
  - `v_circ = √(μ⊕/r_LEO)`, `v_peri = √(v∞² + 2μ⊕/r_LEO)`
  - `Δv_LEO = v_peri − v_circ`, report `Δv_total ≈ Δv_LEO + Δv₂`.

**Intercept logic in viewer:** mark **Intercept** when screen‑space separation ≤ `max(10 px, 0.003 AU × scale)`, otherwise mark **Closest approach** at `t=TOF`. When `r₂≈r₁`, only the **drawing** target orbit is offset by `+0.12 AU` (physics unchanged).

---

## JSON Contract (Planner → Viewer)

Per NEO (subset):
```jsonc
\"intercept_plan\": {
  \"r1_AU\": 1.0,
  \"r2_AU\": <float>,
  \"tof_days\": <float>,
  \"dv_from_LEO_m_s\": <float>,
  \"dv_arrive_heliocentric_m_s\": <float>,
  \"dv_total_m_s\": <float>,
  \"departure_utc\": \"YYYY-MM-DDTHH:MM:SSZ\",
  \"arrival_utc\": \"YYYY-MM-DDTHH:MM:SSZ\",
  \"rolled_forward\": true|false
}
```

---

## Troubleshooting

- **`ModuleNotFoundError: PySide6`** → activate venv, `pip install -r requirements.txt`.
- **macOS PEP 668** → always use a venv; avoid Homebrew’s system Python for `pip` installs.
- **Qt crash on reload** → already handled by recreating scene items (just upgrade to this version).

---

## Roadmap (Short)
- Pull **full Kepler elements** from JPL SBDB/Horizons; propagate true 3D orbits.
- Add a **Lambert solver** for epoch‑accurate windows and Δv; replace simple phasing.
- Visualize plane changes, node markers, and Δv vectors; export PNG/MP4.

---

## License
Classical astrodynamics formulas (Hohmann, Kepler). Add algorithm citations (e.g., Izzo) if/when a Lambert solver is introduced.
