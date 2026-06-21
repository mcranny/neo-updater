# Mathematical model and validation

This document defines the frames, equations, search procedure, validation checks, and limitations used by Asteroid Intercept Planner.

## 1. Units and reference frame

Internal dynamics use SI units:

- distance: metres;
- velocity: metres per second;
- time: seconds;
- solar gravitational parameter: $\mu = 1.32712440018 \times 10^{20}\ \mathrm{m^3/s^2}$.

Times from JPL are interpreted as TDB Julian dates. Asteroid elements and computed state vectors use the ecliptic J2000 frame. Skyfield's ICRF state is rotated about the x-axis by the J2000 mean obliquity before use.

## 2. Asteroid state propagation

For semi-major axis $a$, eccentricity $e$, and solar parameter $\mu$, the mean motion is

$$n = \sqrt{\mu / a^3}.$$

If SBDB provides a mean anomaly $M_0$ at epoch $t_0$, then

$$M(t) = \big(M_0 + n(t - t_0)\big) \bmod 2\pi.$$

If a perihelion epoch $t_p$ is available, the equivalent expression is

$$M(t) = \big(n(t - t_p)\big) \bmod 2\pi.$$

For an elliptical orbit, eccentric anomaly $E$ is found with Newton iteration on Kepler's equation:

$$E - e \sin(E) = M.$$

The perifocal position and velocity are then rotated by $R_3(\Omega)\, R_1(i)\, R_3(\omega)$ into ecliptic J2000.

## 3. Earth state

The default `analytic` provider uses standard low-precision J2000 osculating elements and Kepler propagation. This mode is deterministic, small, and suitable for repeatable transfer sizing.

With `EPHEMERIS_MODE=skyfield`, the application instead computes Earth's heliocentric state from a JPL Development Ephemeris kernel. `JPL_EPHEMERIS_PATH` can point to a local kernel; otherwise Skyfield obtains `de440s.bsp` through its normal loader.

The ephemeris choice materially affects precise departure velocity. Results should always report the selected mode.

## 4. Lambert formulation

Given departure position $\mathbf{r}_1$, arrival position $\mathbf{r}_2$, transfer angle $\Delta\theta$, and time of flight $\Delta t$, define

$$A = \sin(\Delta\theta)\sqrt{\dfrac{|\mathbf{r}_1||\mathbf{r}_2|}{1 - \cos(\Delta\theta)}}.$$

The universal-variable solution uses Stumpff functions $C(z)$ and $S(z)$ and

$$y(z) = |\mathbf{r}_1| + |\mathbf{r}_2| + A\,\frac{zS(z) - 1}{\sqrt{C(z)}},$$

$$F(z) = \left(\frac{y(z)}{C(z)}\right)^{3/2} S(z) + A\sqrt{y(z)} - \sqrt{\mu}\,\Delta t.$$

The code scans the finite single-revolution domain for a sign-changing bracket and then applies bisection. Invalid regions where $y(z) \le 0$ are skipped; they are not treated as valid infinite-time samples.

After obtaining the root,

$$f = 1 - \frac{y}{|\mathbf{r}_1|}, \qquad g = A\sqrt{\frac{y}{\mu}}, \qquad \dot{g} = 1 - \frac{y}{|\mathbf{r}_2|},$$

$$\mathbf{v}_1 = \frac{\mathbf{r}_2 - f\,\mathbf{r}_1}{g}, \qquad \mathbf{v}_2 = \frac{\dot{g}\,\mathbf{r}_2 - \mathbf{r}_1}{g}.$$

Both prograde and retrograde branches are attempted. The first branch that solves and passes validation is eligible for the search.

## 5. Independent endpoint validation

Each Lambert departure state is propagated over $\Delta t$ with the universal Kepler equation:

$$F(\chi) = \frac{\mathbf{r}_0 \cdot \mathbf{v}_0}{\sqrt{\mu}}\,\chi^2\, C(\alpha\chi^2) + \big(1 - \alpha|\mathbf{r}_0|\big)\,\chi^3\, S(\alpha\chi^2) + |\mathbf{r}_0|\chi - \sqrt{\mu}\,\Delta t,$$

$$\alpha = \frac{2}{|\mathbf{r}_0|} - \frac{|\mathbf{v}_0|^2}{\mu}.$$

The propagated endpoint must agree with the requested $\mathbf{r}_2$ within

$$\max\left(10\ \mathrm{km},\ 10^{-7}|\mathbf{r}_2|\right).$$

This check is intentionally separate from the Lambert time equation. It caught and prevents two subtle implementation errors that can otherwise generate smooth-looking but incorrect transfer curves:

1. dividing the Lambert $y(z)$ expression by $C(z)$ instead of $\sqrt{C(z)}$; and
2. using the radial-distance expression as the universal Kepler residual.

Accepted transfers are sampled at uniform time intervals in all three ecliptic J2000 coordinates. The viewer therefore interpolates simulation time along the validated conic rather than moving at a visually convenient constant arc-length speed.

## 6. Search objective

For each close approach, the application evaluates a configured grid of arrival offsets and times of flight. For a candidate solution,

$$\Delta v_{\text{depart}} = |\mathbf{v}_1 - \mathbf{v}_{\text{Earth}}(t_1)|, \qquad \Delta v_{\text{arrive}} = |\mathbf{v}_{\text{asteroid}}(t_2) - \mathbf{v}_2|,$$

$$\Delta v_{\text{total}} = \Delta v_{\text{depart}} + \Delta v_{\text{arrive}}.$$

The stored plan minimizes $\Delta v_{\text{total}}$ over the sampled grid. This is a rendezvous-style velocity-matching metric, not merely a geometric flyby intercept.

## 7. C3 and parking-orbit estimate

The heliocentric departure mismatch is treated as an approximate Earth-relative hyperbolic excess speed:

$$C_3 \approx v_\infty^2 \approx \Delta v_{\text{depart}}^2.$$

For parking-orbit radius $r = R_{\text{Earth}} + h$,

$$v_{\text{circular}} = \sqrt{\mu_{\text{Earth}} / r}, \qquad v_{\text{escape}} = \sqrt{2\mu_{\text{Earth}} / r},$$

$$\Delta v_{\text{LEO}} = \sqrt{v_\infty^2 + v_{\text{escape}}^2} - v_{\text{circular}}.$$

This patched-conic mapping is useful for first-order comparison only. A rigorous Earth-departure design must transform the heliocentric asymptote into the correct geocentric frame and include launch-site and finite-burn constraints.

## 8. Automated checks

The test suite verifies:

- Kepler-equation residuals;
- one-period closure of a circular Earth orbit;
- Lambert propagation to a known Earth-to-Mars-like endpoint;
- finite, physically plausible states and 3D orbit curves for all eight planets;
- mission selection data and legacy-path compatibility;
- rejection of nonpositive transfer time;
- idempotent database upserts;
- foreign-key-compatible viewer export;
- read-only SQL enforcement; and
- dashboard, detail, API, and health routes.

## 9. Limitations

The current model omits:

- n-body and relativistic perturbations;
- covariance propagation and encounter uncertainty;
- multi-revolution Lambert branches;
- continuous launch-window optimization;
- planetary sphere-of-influence transitions;
- launch vehicle performance and finite burns;
- terminal guidance, capture, and proximity operations.

The results are appropriate for preliminary comparison, visualization, and software demonstration—not flight decisions.

## References

- Bate, Mueller, and White, *Fundamentals of Astrodynamics*.
- Battin, *An Introduction to the Mathematics and Methods of Astrodynamics*.
- Vallado, *Fundamentals of Astrodynamics and Applications*.
- JPL Solar System Dynamics API documentation.
