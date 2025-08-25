# Math & Methods

This document expands the math behind the planner and viewer.

---

## 1. Units, Frames, Conventions

- **Canonical solar units**: AU, year, \( \mu_\odot = 4\pi^2 \). Then
  \[
  P^2=a^3,\qquad n=\frac{2\pi}{a^{3/2}}.
  \]
- **SI** for Δv and LEO injection reporting.
- States are heliocentric in the ecliptic frame; the viewer projects to **XY**.

Angles are in radians internally; JSON stores degrees where indicated.

---

## 2. Circular Hohmann (Baseline Sizing)

Given departure radius \( r_1 \) and target radius \( r_2 \) (AU):

**Transfer ellipse**
\[
a_t=\frac{r_1+r_2}{2},\qquad
b_t=\sqrt{r_1 r_2},\qquad
e_t=\sqrt{1-\left(\frac{b_t}{a_t}\right)^2}.
\]

**Speeds**
\[
v_1=\sqrt{\mu/r_1},\quad
v_2=\sqrt{\mu/r_2},\quad
v_p=\sqrt{\mu\left(2/r_1 - 1/a_t\right)},\quad
v_a=\sqrt{\mu\left(2/r_2 - 1/a_t\right)}.
\]

**Δv (Sun-centric)**
\[
\Delta v_1=|v_p-v_1|,\qquad
\Delta v_2=|v_2-v_a|.
\]

**Time of flight**
\[
\mathrm{TOF}=\pi\sqrt{a_t^3/\mu}.
\]

**Viewer phasing (legacy):** For outward transfers,
\[
\theta_T(0)=\pi - n_2 \cdot \mathrm{TOF} \pmod{2\pi},
\]
so the target reaches the end of the half-ellipse at TOF.

---

## 3. Patched-Conic LEO Escape

Let \( r_L \) be LEO radius and \( v_L=\sqrt{\mu_\oplus/r_L} \). With target \( v_\infty \) relative to Earth:
\[
\Delta v_{\text{LEO}}=\sqrt{v_\infty^2+\left(\sqrt{2}\,v_L\right)^2}-v_L.
\]
In baseline Hohmann, \( v_\infty\approx\Delta v_1 \). This is a sizing-level approximation; true injection depends on parking orbit geometry, steering losses, and burn strategy.

---

## 4. Elements → State Vectors

Osculating elements at epoch \( t_0 \):
\[
(a,e,i,\Omega,\omega,M_0)_{t_0}.
\]

**Propagate mean anomaly:**
\[
M(t)=M_0 + n(t-t_0),\qquad n=\sqrt{\frac{\mu}{a^3}}.
\]

**Solve Kepler:**
\[
M=E-e\sin E
\]
(Newton: \(E_0=M\) if \(e<0.8\), else \(E_0=\pi\)).

**True anomaly and PQW state:**
\[
\nu = 2\arctan\!\Big(\sqrt{\tfrac{1+e}{1-e}}\tan\tfrac{E}{2}\Big),\qquad
p=a(1-e^2).
\]
\[
\mathbf{r}_{\text{pf}}=
\begin{bmatrix}
p\cos\nu/(1+e\cos\nu)\\
p\sin\nu/(1+e\cos\nu)\\
0
\end{bmatrix},\quad
\mathbf{v}_{\text{pf}}=
\sqrt{\frac{\mu}{p}}
\begin{bmatrix}
-\sin\nu\\
e+\cos\nu\\
0
\end{bmatrix}.
\]

**Rotate PQW→IJK (ecliptic):**
\[
\mathbf{r}=\mathbf{R}_3(\Omega)\mathbf{R}_1(i)\mathbf{R}_3(\omega)\,\mathbf{r}_{\text{pf}},\qquad
\mathbf{v}=\mathbf{R}_3(\Omega)\mathbf{R}_1(i)\mathbf{R}_3(\omega)\,\mathbf{v}_{\text{pf}}.
\]

---

## 5. Lambert (Universal Variables)

Given \( \mathbf{r}_1 \) at \( t_1 \), \( \mathbf{r}_2 \) at \( t_2=t_1+\Delta t \), find \( \mathbf{v}_1,\mathbf{v}_2 \) such that propagation from \( (\mathbf{r}_1,\mathbf{v}_1) \) reaches \( \mathbf{r}_2 \) in \( \Delta t \).

Define
\[
\Delta\theta=\arccos\left(\frac{\mathbf{r}_1\cdot\mathbf{r}_2}{\|\mathbf{r}_1\|\|\mathbf{r}_2\|}\right),
\quad
A=\sin\Delta\theta\sqrt{\frac{r_1 r_2}{1-\cos\Delta\theta}}.
\]

Use Stumpff functions
\[
C(z)=
\begin{cases}
\frac{1-\cos\sqrt{z}}{z}, & z>0\\[4pt]
\frac{1}{2}, & z=0\\[4pt]
\frac{\cosh\sqrt{-z}-1}{-z}, & z<0
\end{cases}\!,
\quad
S(z)=
\begin{cases}
\frac{\sqrt{z}-\sin\sqrt{z}}{z^{3/2}}, & z>0\\[4pt]
\frac{1}{6}, & z=0\\[4pt]
\frac{\sinh\sqrt{-z}-\sqrt{-z}}{(-z)^{3/2}}, & z<0
\end{cases}\!.
\]

Then
\[
y(z)=r_1+r_2 + A\frac{zS(z)-1}{\sqrt{C(z)}},\qquad
F(z)=\left(\frac{y}{C}\right)^{3/2}S + A\sqrt{y} - \sqrt{\mu}\,\Delta t.
\]

Solve \(F(z)=0\) (secant/Brent-like) for the chosen branch (short-way **prograde** by default). With \(y(z)\), compute Lagrange coefficients:
\[
f=1-\frac{y}{r_1},\quad
g=\frac{A\sqrt{y}}{\sqrt{\mu}},\quad
\dot g=1-\frac{y}{r_2}.
\]
Then
\[
\mathbf{v}_1=\frac{\mathbf{r}_2-f\,\mathbf{r}_1}{g},\qquad
\mathbf{v}_2=\frac{\dot g\,\mathbf{r}_2-\mathbf{r}_1}{g}.
\]

**Residual check:** Propagate \((\mathbf{r}_1,\mathbf{v}_1)\) by \(\Delta t\) (Section 6) and verify it matches \(\mathbf{r}_2\) within tolerance.

---

## 6. Universal-Variable Kepler Propagation

Given \((\mathbf{r}_0,\mathbf{v}_0)\) and \(\Delta t\), define
\[
\alpha=\frac{2}{\|\mathbf{r}_0\|}-\frac{\|\mathbf{v}_0\|^2}{\mu}.
\]
Iterate on the universal anomaly \(\chi\) solving the time-of-flight equation with \(C(z),S(z)\) where \(z=\alpha\chi^2\). Then
\[
f=1-\frac{\chi^2 C}{\|\mathbf{r}_0\|},\quad
g=\Delta t - \frac{\chi^3 S}{\sqrt{\mu}}.
\]
\[
\mathbf{r}=f\,\mathbf{r}_0+g\,\mathbf{v}_0,\qquad
\mathbf{v}=\dot f\,\mathbf{r}_0+\dot g\,\mathbf{v}_0,
\]
with
\[
\dot f=\frac{\sqrt{\mu}}{\|\mathbf{r}\|\|\mathbf{r}_0\|}(\!z S-1)\chi,\qquad
\dot g=1-\frac{\chi^2 C}{\|\mathbf{r}\|}.
\]

Works for elliptic, parabolic (limit), and hyperbolic cases.

---

## 7. Δv Accounting

With Earth \( \mathbf{v}_\oplus(t_1) \) and asteroid \( \mathbf{v}_A(t_2) \):
\[
\Delta v_{\text{depart}}=\|\mathbf{v}_1-\mathbf{v}_\oplus(t_1)\|,\qquad
\Delta v_{\text{arrive}}=\|\mathbf{v}_2-\mathbf{v}_A(t_2)\|.
\]
Optionally expose launcher **C3** via \( C_3=v_\infty^2=\Delta v_{\text{depart}}^2 \).

Total mission Δv for sizing may combine LEO injection and arrival:
\[
\Delta v_{\text{mission}}\approx \Delta v_{\text{LEO}} + \Delta v_{\text{arrive}}.
\]

---

## 8. Plane-Change Sizing (Future)

For inclination differences:
- **At node**: \( \Delta v_\perp \approx 2v\sin(\Delta i/2) \) at speed \(v\).
- **Split-burn**: divide the plane change between departure and arrival to reduce peak cost.
- **Coupled solves**: include out-of-plane components in Lambert or via DSMs.

---

## 9. Window Search (Future)

- Sweep \(t_1\) (e.g., ±2 years) and TOF grid (e.g., 60–360 days).
- Solve Lambert per sample; record \( \Delta v \), \(C_3\), constraints.
- Pareto-filter (Δv, TOF); report top candidates with margins.

---

## 10. Validation

- Kepler \(M(E,e)=E-e\sin E\), round-trip \(E\leftrightarrow\nu\).
- PQW→IJK rotation orthonormality and inverse consistency.
- Propagation residual \( \|\mathbf{r}(t_1+\Delta t)-\mathbf{r}_2\| \).
- JSON schema snapshots and backward-compat loading in the viewer.

---

## 11. References

- Prussing & Conway, *Orbital Mechanics*.  
- Battin, *An Introduction to the Mathematics and Methods of Astrodynamics*.  
- Vallado, *Fundamentals of Astrodynamics and Applications*.  
- Izzo, “Revisiting Lambert’s Problem”.
