# Math & Methods
---

## 1. Units, Frames, Conventions

- **Canonical solar units** for reasoning: AU, year, with $\mu_\odot = 4\pi^2$. Then $P^2 = a^3$ and $n = \dfrac{2\pi}{a^{3/2}}$.
- **SI** for $\Delta v$ and LEO injection reporting.
- States are heliocentric in the ecliptic frame; the viewer projects to **XY**.

Angles are in radians internally; JSON stores degrees where indicated.

---

## 2. Circular Hohmann (Baseline Sizing)

Given departure radius $r_1$ and target radius $r_2$ (AU):

**Transfer ellipse**

$$
a_t = \frac{r_1 + r_2}{2}, \qquad
b_t = \sqrt{r_1 r_2}, \qquad
e_t = \sqrt{1 - \left(\frac{b_t}{a_t}\right)^2}
$$

**Speeds**

$$
v_1 = \sqrt{\mu/r_1}, \qquad
v_2 = \sqrt{\mu/r_2}, \qquad
v_p = \sqrt{\mu\!\left(2/r_1 - 1/a_t\right)}, \qquad
v_a = \sqrt{\mu\!\left(2/r_2 - 1/a_t\right)}
$$

**Sun-centric $\Delta v$**

$$
\Delta v_1 = \lvert v_p - v_1 \rvert, \qquad
\Delta v_2 = \lvert v_2 - v_a \rvert
$$

**Time of flight**

$$
\mathrm{TOF} = \pi \sqrt{\frac{a_t^3}{\mu}}
$$

**Viewer phasing (legacy)**

For outward transfers:

$$
\theta_T(0) = \pi - n_2 \cdot \mathrm{TOF} \pmod{2\pi}
$$

so the target reaches the end of the half-ellipse at TOF.

---

## 3. Patched-Conic LEO Escape

Let $r_L$ be LEO radius and $v_L = \sqrt{\mu_\oplus/r_L}$. With target $v_\infty$ relative to Earth:

$$
\Delta v_{\text{LEO}} = \sqrt{\,v_\infty^2 + \left(\sqrt{2}\,v_L\right)^2\,} - v_L
$$

In baseline Hohmann, $v_\infty \approx \Delta v_1$. This is a sizing-level approximation; real injection depends on parking-orbit geometry, steering losses, and burn strategy.

---

## 4. Elements → State Vectors

Osculating elements at epoch $t_0$: $(a,e,i,\Omega,\omega,M_0)_{t_0}$.

**Mean anomaly at time $t$**

$$
M(t) = M_0 + n\,(t - t_0), \qquad n = \sqrt{\frac{\mu}{a^3}}
$$

**Kepler’s equation (elliptic)**

$$
M = E - e\sin E
$$

(Newton: $E_0 = M$ if $e < 0.8$, else $E_0 = \pi$.)

**True anomaly and PQW state**

$$
\nu = 2\arctan\!\Big(\sqrt{\tfrac{1+e}{1-e}}\tan\tfrac{E}{2}\Big), \qquad
p = a(1-e^2)
$$

$$
\mathbf{r}_{\mathrm{pf}} =
\begin{bmatrix}
\dfrac{p\cos\nu}{1+e\cos\nu}\\
\dfrac{p\sin\nu}{1+e\cos\nu}\\
0
\end{bmatrix}
$$

$$
\mathbf{v}_{\mathrm{pf}} = \sqrt{\dfrac{\mu}{p}}\,
\begin{bmatrix}
-\sin\nu\\
e+\cos\nu\\
0
\end{bmatrix}
$$

**Rotate PQW→IJK (ecliptic)**

$$
\mathbf{r} = R_3(\Omega)\,R_1(i)\,R_3(\omega)\,\mathbf{r}_{\mathrm{pf}}
$$

$$
\mathbf{v} = R_3(\Omega)\,R_1(i)\,R_3(\omega)\,\mathbf{v}_{\mathrm{pf}}
$$

---

## 5. Lambert (Universal Variables)

Given $\mathbf{r}_1$ at $t_1$ and $\mathbf{r}_2$ at $t_2=t_1+\Delta t$, find $\mathbf{v}_1,\mathbf{v}_2$ such that propagation from $(\mathbf{r}_1,\mathbf{v}_1)$ reaches $\mathbf{r}_2$ in $\Delta t$ (short-way **prograde** by default).

Geometry:

$$
\Delta\theta=\arccos\!\left(\frac{\mathbf{r}_1\!\cdot\!\mathbf{r}_2}{\lVert\mathbf{r}_1\rVert\,\lVert\mathbf{r}_2\rVert}\right),\qquad
A=\sin\Delta\theta\,\sqrt{\frac{r_1 r_2}{1-\cos\Delta\theta}}
$$

Stumpff functions:

$$
C(z)=
\begin{cases}
\dfrac{1-\cos\sqrt{z}}{z}, & \text{if } z>0,\\
\dfrac{1}{2}, & \text{if } z=0,\\
\dfrac{\cosh\sqrt{-z}-1}{-z}, & \text{if } z<0
\end{cases}
\qquad
S(z)=
\begin{cases}
\dfrac{\sqrt{z}-\sin\sqrt{z}}{z^{3/2}}, & \text{if } z>0,\\
\dfrac{1}{6}, & \text{if } z=0,\\
\dfrac{\sinh\sqrt{-z}-\sqrt{-z}}{(-z)^{3/2}}, & \text{if } z<0
\end{cases}
$$

Define

$$
y(z)=r_1+r_2 + A\,\frac{z\,S(z)-1}{\sqrt{C(z)}}
$$

and solve the scalar equation

$$
F(z)=\left(\frac{y}{C}\right)^{3/2}S + A\sqrt{y} - \sqrt{\mu}\,\Delta t = 0
$$

With the root $z^\*$,

$$
f=1-\frac{y}{r_1},\qquad
g=\frac{A\sqrt{y}}{\sqrt{\mu}},\qquad
\dot g=1-\frac{y}{r_2}
$$

and

$$
\mathbf{v}_1=\frac{\mathbf{r}_2 - f\,\mathbf{r}_1}{g},\qquad
\mathbf{v}_2=\frac{\dot g\,\mathbf{r}_2 - \mathbf{r}_1}{g}
$$

**Residual check**: propagate $(\mathbf{r}_1,\mathbf{v}_1)$ by $\Delta t$ (Section 6) and verify $\lVert\mathbf{r}(t_1+\Delta t)-\mathbf{r}_2\rVert$ is within tolerance.

---

## 6. Universal-Variable Kepler Propagation

Given $(\mathbf{r}_0,\mathbf{v}_0)$ and $\Delta t$, define

$$
\alpha=\frac{2}{\lVert\mathbf{r}_0\rVert}-\frac{\lVert\mathbf{v}_0\rVert^2}{\mu},\qquad z=\alpha\,\chi^2
$$

Iterate on the universal anomaly $\chi$ using $C(z),S(z)$ to satisfy the time-of-flight equation. Then

$$
f=1-\frac{\chi^2 C}{\lVert\mathbf{r}_0\rVert},\qquad
g=\Delta t - \frac{\chi^3 S}{\sqrt{\mu}}
$$

$$
\mathbf{r}=f\,\mathbf{r}_0+g\,\mathbf{v}_0,\qquad
\mathbf{v}=\dot{f}\,\mathbf{r}_0+\dot{g}\,\mathbf{v}_0
$$

with

$$
\dot f=\frac{\sqrt{\mu}}{\lVert\mathbf{r}\rVert\,\lVert\mathbf{r}_0\rVert}\,(\,zS-1\,)\,\chi,\qquad
\dot g=1-\frac{\chi^2 C}{\lVert\mathbf{r}\rVert}
$$

This works for elliptic, parabolic (limit), and hyperbolic cases.

---

## 7. $\Delta v$ Accounting

Let $\mathbf v_\oplus(t_1)$ be Earth’s heliocentric velocity at departure, and $\mathbf v_A(t_2)$ the asteroid’s at arrival.

$$
\Delta v_{\mathrm{depart}}
= \sqrt{\,(\mathbf v_1 - \mathbf v_\oplus(t_1))\cdot(\mathbf v_1 - \mathbf v_\oplus(t_1))\,}
$$

$$
\Delta v_{\mathrm{arrive}}
= \sqrt{\,(\mathbf v_2 - \mathbf v_A(t_2))\cdot(\mathbf v_2 - \mathbf v_A(t_2))\,}
$$

Optionally expose launcher $C_3$ via $C_3 = v_\infty^{2} = (\Delta v_{\mathrm{depart}})^{2}$.

Total mission sizing (simple):

$$
\Delta v_{\mathrm{mission}} \approx \Delta v_{\mathrm{LEO}} + \Delta v_{\mathrm{arrive}}
$$

---

## 8. Plane-Change Sizing (Future)

For inclination differences $\Delta i$:

- **At node**: $\Delta v_\perp \approx 2v\sin(\Delta i/2)$ at speed $v$.
- **Split-burn**: divide plane change between departure/arrival to reduce peak cost.
- **Coupled solves**: include out-of-plane components in Lambert or via DSMs.

---

## 9. Window Search (Future)

- Sweep $t_1$ (e.g., $\pm 2$ years) and TOF grid (e.g., 60–360 days).
- Solve Lambert per sample; record $\Delta v$, $C_3$, constraints.
- Pareto-filter (Δv, TOF); report top candidates with margins.

---

## 10. Validation

- Kepler $M(E,e)=E-e\sin E$; round-trip $E\leftrightarrow\nu$.
- PQW→IJK rotation orthonormality and inverse consistency.
- Propagation residual $\lVert\mathbf{r}(t_1+\Delta t)-\mathbf{r}_2\rVert$.
- JSON schema snapshots and backward-compat loading in the viewer.

---

## 11. References

- Prussing & Conway, *Orbital Mechanics*  
- Battin, *An Introduction to the Mathematics and Methods of Astrodynamics*  
- Vallado, *Fundamentals of Astrodynamics and Applications*  
- Izzo, “Revisiting Lambert’s Problem”
