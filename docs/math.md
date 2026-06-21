# Mathematical model and validation

This document defines the frames, equations, search procedure, validation checks, and limitations used by Asteroid Intercept Planner.

## 1. Units and reference frame

Internal dynamics use SI units:

- distance: metres;
- velocity: metres per second;
- time: seconds;
- solar gravitational parameter: `1.32712440018e20 m³/s²`.

Times from JPL are interpreted as TDB Julian dates. Asteroid elements and computed state vectors use the ecliptic J2000 frame. Skyfield's ICRF state is rotated about the x-axis by the J2000 mean obliquity before use.

## 2. Asteroid state propagation

For semi-major axis `a`, eccentricity `e`, and solar parameter `μ`, the mean motion is

```text
n = sqrt(μ / a³).
```

If SBDB provides a mean anomaly `M₀` at epoch `t₀`, then

```text
M(t) = wrap₂π(M₀ + n(t - t₀)).
```

If a perihelion epoch `tₚ` is available, the equivalent expression is

```text
M(t) = wrap₂π(n(t - tₚ)).
```

For an elliptical orbit, eccentric anomaly `E` is found with Newton iteration on Kepler's equation:

```text
E - e sin(E) = M.
```

The perifocal position and velocity are then rotated by `R₃(Ω) R₁(i) R₃(ω)` into ecliptic J2000.

## 3. Earth state

The default `analytic` provider uses standard low-precision J2000 osculating elements and Kepler propagation. This mode is deterministic, small, and suitable for repeatable transfer sizing.

With `EPHEMERIS_MODE=skyfield`, the application instead computes Earth's heliocentric state from a JPL Development Ephemeris kernel. `JPL_EPHEMERIS_PATH` can point to a local kernel; otherwise Skyfield obtains `de440s.bsp` through its normal loader.

The ephemeris choice materially affects precise departure velocity. Results should always report the selected mode.

## 4. Lambert formulation

Given departure position **r₁**, arrival position **r₂**, transfer angle `Δθ`, and time of flight `Δt`, define

```text
A = sin(Δθ) sqrt(|r₁||r₂| / (1 - cos(Δθ))).
```

The universal-variable solution uses Stumpff functions `C(z)` and `S(z)` and

```text
y(z) = |r₁| + |r₂| + A (z S(z) - 1) / sqrt(C(z)),

F(z) = (y(z) / C(z))^(3/2) S(z) + A sqrt(y(z)) - sqrt(μ) Δt.
```

The code scans the finite single-revolution domain for a sign-changing bracket and then applies bisection. Invalid regions where `y(z) <= 0` are skipped; they are not treated as valid infinite-time samples.

After obtaining the root,

```text
f    = 1 - y / |r₁|,
g    = A sqrt(y / μ),
gdot = 1 - y / |r₂|,

v₁ = (r₂ - f r₁) / g,
v₂ = (gdot r₂ - r₁) / g.
```

Both prograde and retrograde branches are attempted. The first branch that solves and passes validation is eligible for the search.

## 5. Independent endpoint validation

Each Lambert departure state is propagated over `Δt` with the universal Kepler equation:

```text
F(χ) = (r₀·v₀ / sqrt(μ)) χ² C(αχ²)
     + (1 - α|r₀|) χ³ S(αχ²)
     + |r₀|χ - sqrt(μ)Δt,

α = 2/|r₀| - |v₀|²/μ.
```

The propagated endpoint must agree with the requested **r₂** within

```text
max(10 km, 10⁻⁷ |r₂|).
```

This check is intentionally separate from the Lambert time equation. It caught and prevents two subtle implementation errors that can otherwise generate smooth-looking but incorrect transfer curves:

1. dividing the Lambert `y(z)` expression by `C(z)` instead of `sqrt(C(z))`; and
2. using the radial-distance expression as the universal Kepler residual.

## 6. Search objective

For each close approach, the application evaluates a configured grid of arrival offsets and times of flight. For a candidate solution,

```text
Δv_depart = |v₁ - v_Earth(t₁)|,
Δv_arrive = |v_asteroid(t₂) - v₂|,
Δv_total  = Δv_depart + Δv_arrive.
```

The stored plan minimizes `Δv_total` over the sampled grid. This is a rendezvous-style velocity-matching metric, not merely a geometric flyby intercept.

## 7. C3 and parking-orbit estimate

The heliocentric departure mismatch is treated as an approximate Earth-relative hyperbolic excess speed:

```text
C3 ≈ v∞² ≈ Δv_depart².
```

For parking-orbit radius `r = R_Earth + h`,

```text
v_circular = sqrt(μ_Earth / r),
v_escape   = sqrt(2 μ_Earth / r),
Δv_LEO     = sqrt(v∞² + v_escape²) - v_circular.
```

This patched-conic mapping is useful for first-order comparison only. A rigorous Earth-departure design must transform the heliocentric asymptote into the correct geocentric frame and include launch-site and finite-burn constraints.

## 8. Automated checks

The test suite verifies:

- Kepler-equation residuals;
- one-period closure of a circular Earth orbit;
- Lambert propagation to a known Earth-to-Mars-like endpoint;
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
