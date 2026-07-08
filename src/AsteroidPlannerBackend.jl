module AsteroidPlannerBackend

using JSON
using LinearAlgebra

const AU_M = 149_597_870_700.0
const DAY_S = 86_400.0
const MU_SUN = 1.32712440018e20
const MU_EARTH = 3.986004418e14
const R_EARTH_M = 6_378_137.0
const G_M3_KG_S2 = 6.67430e-11
const JD_J2000 = 2_451_545.0
const DEFAULT_DENSITY_KG_M3 = 2_000.0
const DEFAULT_ALBEDO = 0.14

deg2rad(x) = Float64(x) * pi / 180.0
wrap_2pi(x) = mod(Float64(x), 2.0 * pi)

function finite_positive(x)
    return x !== nothing && isfinite(Float64(x)) && Float64(x) > 0.0
end

function maybe_float(value)
    value === nothing && return nothing
    try
        return Float64(value)
    catch
        return nothing
    end
end

function get_float(dict::AbstractDict, key::AbstractString)
    return haskey(dict, key) ? maybe_float(dict[key]) : nothing
end

function kepler_E_from_M(M, e; tol=1e-12, itmax=60)
    M = wrap_2pi(M)
    e = Float64(e)
    e < 1e-12 && return M
    E = e < 0.8 ? M : pi
    for _ in 1:itmax
        f = E - e * sin(E) - M
        fp = 1.0 - e * cos(E)
        dE = -f / fp
        E += dE
        abs(dE) < tol && break
    end
    return E
end

function oe_to_rv(a_m, e, i, raan, argp, nu; mu=MU_SUN)
    a_m = Float64(a_m)
    e = Float64(e)
    p = a_m * (1.0 - e * e)
    p > 0.0 || throw(ArgumentError("orbital elements require positive semilatus rectum"))
    cnu, snu = cos(nu), sin(nu)
    r_pf = [p * cnu / (1.0 + e * cnu), p * snu / (1.0 + e * cnu), 0.0]
    v_pf = [-sqrt(mu / p) * snu, sqrt(mu / p) * (e + cnu), 0.0]

    ci, si = cos(i), sin(i)
    cO, sO = cos(raan), sin(raan)
    co, so = cos(argp), sin(argp)
    R = [
        cO * co - sO * so * ci  -cO * so - sO * co * ci   sO * si
        sO * co + cO * so * ci  -sO * so + cO * co * ci  -cO * si
        so * si                  co * si                   ci
    ]
    return R * r_pf, R * v_pf
end

function stumpff_C(z)
    z = Float64(z)
    abs(z) < 1e-8 && return 0.5 - z / 24.0 + z * z / 720.0
    if z > 0.0
        s = sqrt(z)
        return (1.0 - cos(s)) / z
    end
    s = sqrt(-z)
    s <= 700.0 || throw(OverflowError("Stumpff C overflow"))
    return (cosh(s) - 1.0) / (-z)
end

function stumpff_S(z)
    z = Float64(z)
    abs(z) < 1e-8 && return (1.0 / 6.0) - z / 120.0 + z * z / 5040.0
    if z > 0.0
        s = sqrt(z)
        return (s - sin(s)) / s^3
    end
    s = sqrt(-z)
    s <= 700.0 || throw(OverflowError("Stumpff S overflow"))
    return (sinh(s) - s) / s^3
end

function kepler_universal_propagate(r0, v0, dt; mu=MU_SUN, itmax=80, tol=1e-6)
    r0 = Float64.(collect(r0))
    v0 = Float64.(collect(v0))
    dt = Float64(dt)
    r0n = norm(r0)
    r0n > 0.0 || throw(ArgumentError("kepler_universal_propagate: |r0| is zero"))
    r0_dot_v0 = dot(r0, v0)
    alpha = 2.0 / r0n - dot(v0, v0) / mu
    sqrt_mu = sqrt(mu)

    if abs(alpha) > 1e-12
        chi = sqrt_mu * abs(alpha) * dt
    else
        h = norm(cross(r0, v0))
        p = h * h / mu
        p > 0.0 || throw(ArgumentError("kepler_universal_propagate: degenerate angular momentum"))
        s = 0.5 * (pi / 2.0 - atan(3.0 * sqrt(mu / p^3) * dt))
        w = atan(tan(s)^(1.0 / 3.0))
        chi = sqrt(p) * 2.0 / tan(2.0 * w)
    end

    converged = false
    for _ in 1:itmax
        z = alpha * chi * chi
        C = stumpff_C(z)
        S = stumpff_S(z)
        F = (r0_dot_v0 / sqrt_mu) * chi^2 * C +
            (1.0 - alpha * r0n) * chi^3 * S +
            r0n * chi - sqrt_mu * dt
        dF = (r0_dot_v0 / sqrt_mu) * chi * (1.0 - z * S) +
            (1.0 - alpha * r0n) * chi^2 * C + r0n
        if abs(F) < tol * sqrt_mu
            converged = true
            break
        end
        abs(dF) > 1e-16 || throw(ErrorException("Kepler propagation derivative vanished"))
        chi -= F / dF
    end
    converged || throw(ErrorException("Kepler universal-variable solver did not converge"))

    z = alpha * chi * chi
    C = stumpff_C(z)
    S = stumpff_S(z)
    f = 1.0 - chi^2 * C / r0n
    g = dt - chi^3 * S / sqrt_mu
    r = f .* r0 .+ g .* v0
    rn = norm(r)
    fdot = (sqrt_mu / (rn * r0n)) * (z * S - 1.0) * chi
    gdot = 1.0 - chi^2 * C / rn
    v = fdot .* r0 .+ gdot .* v0
    return r, v
end

function lambert_universal(r1, r2, tof_s; mu=MU_SUN, prograde=true, max_iter=80, tol=1e-8)
    r1 = Float64.(collect(r1))
    r2 = Float64.(collect(r2))
    tof_s = Float64(tof_s)
    r1n = norm(r1)
    r2n = norm(r2)
    r1n > 0.0 && r2n > 0.0 || throw(ArgumentError("Lambert: |r1| or |r2| is zero"))
    tof_s > 0.0 || throw(ArgumentError("Lambert: time of flight must be positive"))
    cosd = clamp(dot(r1, r2) / (r1n * r2n + 1e-16), -1.0, 1.0)
    dtheta = acos(cosd)
    cz = cross(r1, r2)[3]
    if prograde && cz < 0.0
        dtheta = 2.0 * pi - dtheta
    elseif !prograde && cz > 0.0
        dtheta = 2.0 * pi - dtheta
    end

    denom = 1.0 - cosd
    abs(denom) >= 1e-16 || throw(ArgumentError("Lambert geometry singular"))
    A = sin(dtheta) * sqrt(max(0.0, r1n * r2n / denom))
    abs(A) >= 1e-14 || throw(ArgumentError("Lambert geometry singular: A≈0"))
    sqrt_mu = sqrt(mu)

    function y_of_z(z)
        C = stumpff_C(z)
        S = stumpff_S(z)
        Ceff = abs(C) > 1e-16 ? C : 1e-16
        return r1n + r2n + A * ((z * S - 1.0) / sqrt(Ceff))
    end

    function tof_of_z(z)
        C = stumpff_C(z)
        S = stumpff_S(z)
        y = y_of_z(z)
        y > 0.0 || return Inf
        Ceff = abs(C) > 1e-16 ? C : 1e-16
        x = sqrt(max(0.0, y / Ceff))
        return (x^3 * S + A * sqrt(y)) / sqrt_mu
    end

    z_limit = 4.0 * pi^2 - 1e-6
    previous = nothing
    bracket = nothing
    for candidate in range(-z_limit, z_limit; length=1200)
        local value
        try
            value = tof_of_z(candidate) - tof_s
        catch
            continue
        end
        isfinite(value) || continue
        if abs(value) <= tol
            bracket = (Float64(candidate), Float64(candidate))
            break
        end
        if previous !== nothing && previous[2] * value < 0.0
            bracket = (previous[1], Float64(candidate))
            break
        end
        previous = (Float64(candidate), value)
    end
    bracket !== nothing || throw(ErrorException("Lambert: no single-revolution time-of-flight root"))

    z_lo, z_hi = bracket
    z = z_lo
    converged = z_lo == z_hi
    for _ in 1:max_iter
        converged && break
        z = 0.5 * (z_lo + z_hi)
        f_mid = tof_of_z(z) - tof_s
        if abs(f_mid) <= tol
            converged = true
            break
        end
        f_lo = tof_of_z(z_lo) - tof_s
        if f_lo * f_mid <= 0.0
            z_hi = z
        else
            z_lo = z
        end
    end
    if !converged
        converged = abs(z_hi - z_lo) < 1e-12
        z = 0.5 * (z_lo + z_hi)
    end
    converged || throw(ErrorException("Lambert: time-of-flight root did not converge"))

    y = y_of_z(z)
    y > 0.0 || throw(ErrorException("Lambert: y(z*) <= 0 at solution"))
    f = 1.0 - y / r1n
    g = A * sqrt(y / mu)
    gdot = 1.0 - y / r2n
    if abs(g) < 1e-14
        g = copysign(1e-14, g)
    end
    v1 = (r2 .- f .* r1) ./ g
    v2 = (gdot .* r2 .- r1) ./ g
    propagated, _ = kepler_universal_propagate(r1, v1, tof_s; mu=mu)
    residual_m = norm(propagated .- r2)
    tolerance_m = max(10_000.0, 1e-7 * r2n)
    if !isfinite(residual_m) || residual_m > tolerance_m
        throw(ErrorException("Lambert endpoint residual $(residual_m) m exceeds $(tolerance_m) m"))
    end
    return v1, v2
end

function sample_transfer_polyline(r1, v1, tof_s, n=180; mu=MU_SUN)
    points = Vector{Vector{Float64}}()
    for k in 0:n
        t = (k / n) * Float64(tof_s)
        r, _ = kepler_universal_propagate(r1, v1, t; mu=mu)
        push!(points, [r[1] / AU_M, r[2] / AU_M, r[3] / AU_M])
    end
    return points
end

function earth_rv_analytic(jd_tdb)
    semi_major_axis = 1.00000011 * AU_M
    eccentricity = 0.01671022
    inclination = deg2rad(0.00005)
    longitude_ascending_node = deg2rad(-11.26064)
    longitude_perihelion = deg2rad(102.94719)
    argument_perihelion = longitude_perihelion - longitude_ascending_node
    mean_longitude_j2000 = deg2rad(100.46435)
    mean_anomaly_j2000 = mean_longitude_j2000 - longitude_perihelion
    mean_motion = sqrt(MU_SUN / semi_major_axis^3)
    mean_anomaly = wrap_2pi(mean_anomaly_j2000 + mean_motion * (Float64(jd_tdb) - JD_J2000) * DAY_S)
    eccentric_anomaly = kepler_E_from_M(mean_anomaly, eccentricity)
    true_anomaly = 2.0 * atan(
        sqrt(1.0 + eccentricity) * sin(eccentric_anomaly / 2.0),
        sqrt(1.0 - eccentricity) * cos(eccentric_anomaly / 2.0),
    )
    return oe_to_rv(
        semi_major_axis,
        eccentricity,
        inclination,
        longitude_ascending_node,
        argument_perihelion,
        true_anomaly,
    )
end

function asteroid_M_at_time(elements::AbstractDict, jd_tdb)
    a = Float64(elements["a_AU"]) * AU_M
    n = sqrt(MU_SUN / a^3)
    tp = get_float(elements, "tp_jd_tdb")
    if tp !== nothing
        return wrap_2pi(n * ((Float64(jd_tdb) - tp) * DAY_S))
    end
    ma_deg = get_float(elements, "ma_deg")
    epoch = get_float(elements, "epoch_jd_tdb")
    if ma_deg === nothing || epoch === nothing
        return nothing
    end
    return wrap_2pi(deg2rad(ma_deg) + n * (Float64(jd_tdb) - epoch) * DAY_S)
end

function asteroid_rv_from_elements(elements::AbstractDict, jd_tdb)
    a = Float64(elements["a_AU"]) * AU_M
    e = Float64(elements["e"])
    inc = deg2rad(elements["i_deg"])
    raan = deg2rad(elements["raan_deg"])
    argp = deg2rad(elements["argp_deg"])
    M = asteroid_M_at_time(elements, jd_tdb)
    M !== nothing || throw(ArgumentError("insufficient elements to compute mean anomaly"))
    E = kepler_E_from_M(M, e)
    cE, sE = cos(E), sin(E)
    r_peri = [a * (cE - e), a * (sqrt(1.0 - e * e) * sE), 0.0]
    r_mag = a * (1.0 - e * cE)
    fac = sqrt(MU_SUN * a) / r_mag
    v_peri = [-fac * sE, fac * sqrt(1.0 - e * e) * cE, 0.0]

    cO, sO = cos(raan), sin(raan)
    ci, si = cos(inc), sin(inc)
    cw, sw = cos(argp), sin(argp)
    R = [
        cO * cw - sO * sw * ci  -cO * sw - sO * cw * ci   sO * si
        sO * cw + cO * sw * ci  -sO * sw + cO * cw * ci  -cO * si
        sw * si                  cw * si                   ci
    ]
    return R * r_peri, R * v_peri
end

function leo_departure_dv(vinf_kms, leo_alt_m=500e3)
    r = R_EARTH_M + Float64(leo_alt_m)
    v_circ = sqrt(MU_EARTH / r)
    v_esc = sqrt(2.0 * MU_EARTH / r)
    v_hyp = sqrt((Float64(vinf_kms) * 1000.0)^2 + v_esc^2)
    return (v_hyp - v_circ) / 1000.0
end

function diameter_from_H_m(H, albedo=DEFAULT_ALBEDO)
    H = Float64(H)
    albedo = Float64(albedo)
    albedo > 0.0 || throw(ArgumentError("albedo must be positive"))
    return 1329.0 / sqrt(albedo) * 10.0^(-H / 5.0) * 1000.0
end

function physical_radius_m(obj::AbstractDict, elements::AbstractDict, albedo=DEFAULT_ALBEDO)
    ca = get(obj, "next_ca", Dict{String,Any}())
    diameter_km = get_float(ca, "diameter_km")
    if finite_positive(diameter_km)
        return 0.5 * diameter_km * 1000.0, "diameter_km"
    end
    H = get_float(ca, "h")
    H === nothing && (H = get_float(elements, "H"))
    if H !== nothing && isfinite(H)
        return 0.5 * diameter_from_H_m(H, albedo), "H_albedo"
    end
    return nothing, "missing"
end

function circular_orbit_state(radius_m, mu_ast)
    r = [Float64(radius_m), 0.0, 0.0]
    v = [0.0, sqrt(Float64(mu_ast) / Float64(radius_m)), 0.0]
    return r, v
end

function propagate_two_body_relative(r0, v0, duration_s, samples; mu_ast)
    points = Vector{Vector{Float64}}()
    min_radius = Inf
    max_radius = 0.0
    for k in 0:samples
        t = (k / samples) * Float64(duration_s)
        r, _ = kepler_universal_propagate(r0, v0, t; mu=mu_ast)
        radius = norm(r)
        min_radius = min(min_radius, radius)
        max_radius = max(max_radius, radius)
        push!(points, [r[1], r[2], r[3]])
    end
    return points, min_radius, max_radius
end

function heliocentric_orbit_polyline(relative_points, elements::AbstractDict, arrival_jd, duration_s)
    output = Vector{Vector{Float64}}()
    count = max(length(relative_points) - 1, 1)
    for (idx, relative) in enumerate(relative_points)
        # Sample the asteroid ephemeris during the local-orbit validation span.
        # The local orbit itself is modeled in a simple asteroid-centered frame;
        # this overlay is only for viewer context at solar-system scale.
        jd = Float64(arrival_jd) + ((idx - 1) / count) * (Float64(duration_s) / DAY_S)
        r_ast, _ = asteroid_rv_from_elements(elements, jd)
        r = r_ast .+ relative
        push!(output, [r[1] / AU_M, r[2] / AU_M, r[3] / AU_M])
    end
    return output
end

function capture_model(obj::AbstractDict, elements::AbstractDict, arrival_jd, arrival_vinf_kms, dv_sum_mps; density=DEFAULT_DENSITY_KG_M3, albedo=DEFAULT_ALBEDO, final_orbits=10.0)
    flags = String[]
    density = Float64(density)
    albedo = Float64(albedo)
    final_orbits = Float64(final_orbits)
    density > 0.0 || throw(ArgumentError("density must be positive"))
    final_orbits > 0.0 || throw(ArgumentError("final_orbits must be positive"))
    radius_value, source = physical_radius_m(obj, elements, albedo)
    if radius_value === nothing || !finite_positive(radius_value)
        return Dict{String,Any}(
            "model" => "density-sphere-v1",
            "density_kg_m3" => density,
            "albedo" => albedo,
            "diameter_source" => source,
            "stable_final_orbit" => false,
            "stability_flags" => ["missing physical size"],
            "requested_final_orbits" => final_orbits,
            "propagated_final_orbits" => 0.0,
            "orbit_polyline_relative_m" => Any[],
            "orbit_polyline_heliocentric_au" => Any[],
        )
    end

    radius_m = Float64(radius_value)
    volume_m3 = 4.0 / 3.0 * pi * radius_m^3
    mass_kg = density * volume_m3
    mu_ast = G_M3_KG_S2 * mass_kg
    a_ast_m = Float64(elements["a_AU"]) * AU_M
    hill_radius_m = a_ast_m * cbrt(mass_kg / (3.0 * (MU_SUN / G_M3_KG_S2)))
    target_radius_m = 3.0 * radius_m
    target_altitude_m = target_radius_m - radius_m

    if !isfinite(hill_radius_m) || hill_radius_m <= 0.0
        push!(flags, "invalid Hill radius")
    elseif target_radius_m > hill_radius_m / 3.0
        push!(flags, "target orbit outside one-third Hill radius")
    end

    circular_speed_mps = sqrt(mu_ast / target_radius_m)
    escape_speed_mps = sqrt(2.0 * mu_ast / target_radius_m)
    surface_gravity_m_s2 = mu_ast / radius_m^2
    orbit_period_s = 2.0 * pi * sqrt(target_radius_m^3 / mu_ast)
    arrival_vinf_mps = Float64(arrival_vinf_kms) * 1000.0
    periapsis_speed_mps = sqrt(arrival_vinf_mps^2 + escape_speed_mps^2)
    capture_dv_mps = abs(periapsis_speed_mps - circular_speed_mps)
    propagation_s = min(final_orbits * orbit_period_s, 30.0 * DAY_S)
    propagated_final_orbits = propagation_s / orbit_period_s
    samples = 240

    r0, v0 = circular_orbit_state(target_radius_m, mu_ast)
    relative_polyline, min_radius_m, max_radius_m = propagate_two_body_relative(
        r0, v0, propagation_s, samples; mu_ast=mu_ast
    )
    if min_radius_m <= radius_m
        push!(flags, "surface intersection")
    end
    if max_radius_m >= hill_radius_m / 3.0
        push!(flags, "propagation exceeds one-third Hill radius")
    end
    drift_m = max(abs(min_radius_m - target_radius_m), abs(max_radius_m - target_radius_m))
    if drift_m > max(1.0, 1e-6 * target_radius_m)
        push!(flags, "two-body propagation drift")
    end
    stable = isempty(flags)

    return Dict{String,Any}(
        "model" => "density-sphere-v1",
        "density_kg_m3" => density,
        "albedo" => albedo,
        "diameter_source" => source,
        "radius_m" => radius_m,
        "mass_kg" => mass_kg,
        "mu_m3_s2" => mu_ast,
        "hill_radius_m" => hill_radius_m,
        "target_orbit_radius_m" => target_radius_m,
        "target_orbit_altitude_m" => target_altitude_m,
        "surface_gravity_m_s2" => surface_gravity_m_s2,
        "orbit_period_s" => orbit_period_s,
        "circular_speed_mps" => circular_speed_mps,
        "escape_speed_mps" => escape_speed_mps,
        "arrival_vinf_kms" => Float64(arrival_vinf_kms),
        "capture_dv_kms" => capture_dv_mps / 1000.0,
        "dv_total_with_capture_kms" => Float64(dv_sum_mps) / 1000.0 + capture_dv_mps / 1000.0,
        "stable_final_orbit" => stable,
        "stability_flags" => flags,
        "requested_final_orbits" => final_orbits,
        "propagated_final_orbits" => propagated_final_orbits,
        "propagation_days" => propagation_s / DAY_S,
        "propagation_min_radius_m" => min_radius_m,
        "propagation_max_radius_m" => max_radius_m,
        "propagation_radial_drift_m" => drift_m,
        "orbit_polyline_relative_m" => relative_polyline,
        "orbit_polyline_heliocentric_au" => heliocentric_orbit_polyline(relative_polyline, elements, arrival_jd, propagation_s),
    )
end

function plan_intercept_for_object(obj::AbstractDict; tof_days_grid=collect(30:10:180), arrive_offset_hours=[-12, -6, 0, 6, 12], leo_alt_m=500e3, final_orbits=10.0)
    elements = get(obj, "elements", nothing)
    ca = get(obj, "next_ca", nothing)
    (elements isa AbstractDict && ca isa AbstractDict) || return nothing
    haskey(ca, "jd_tdb") && ca["jd_tdb"] !== nothing || return nothing
    jd_ca = Float64(ca["jd_tdb"])
    best = nothing

    for d_off in arrive_offset_hours
        jd_arr = jd_ca + Float64(d_off) / 24.0
        local r_ast_arr, v_ast_arr
        try
            r_ast_arr, v_ast_arr = asteroid_rv_from_elements(elements, jd_arr)
        catch
            continue
        end
        for tof_d in tof_days_grid
            jd_dep = jd_arr - Float64(tof_d)
            r_e_dep, v_e_dep = earth_rv_analytic(jd_dep)
            tof_s = Float64(tof_d) * DAY_S
            solution = nothing
            for prograde in (true, false)
                try
                    solution = lambert_universal(r_e_dep, r_ast_arr, tof_s; mu=MU_SUN, prograde=prograde)
                    break
                catch
                    continue
                end
            end
            solution === nothing && continue
            v1, v2 = solution
            dv_dep = norm(v1 .- v_e_dep)
            dv_arr = norm(v_ast_arr .- v2)
            dv_sum = dv_dep + dv_arr
            if best === nothing || dv_sum < best["dv_sum_mps"]
                polyline = Any[]
                try
                    polyline = sample_transfer_polyline(r_e_dep, v1, tof_s, 180)
                catch
                    polyline = Any[]
                end
                plan = Dict{String,Any}(
                    "depart_jd_tdb" => jd_dep,
                    "arrive_jd_tdb" => jd_arr,
                    "tof_days" => Float64(tof_d),
                    "dv_depart_kms" => dv_dep / 1000.0,
                    "dv_arrive_kms" => dv_arr / 1000.0,
                    "dv_sum_mps" => dv_sum,
                    "c3_km2_s2" => (dv_dep / 1000.0)^2,
                    "vinf_kms" => dv_dep / 1000.0,
                    "leo_dv_kms" => leo_departure_dv(dv_dep / 1000.0, leo_alt_m),
                    "lambert_polyline_xyz_au" => polyline,
                )
                plan["capture"] = capture_model(obj, elements, jd_arr, dv_arr / 1000.0, dv_sum; final_orbits=final_orbits)
                best = plan
            end
        end
    end
    return best
end

function attach_intercepts!(payload::AbstractDict; tof_days_grid=collect(30:10:180), arrive_offset_hours=[-12, -6, 0, 6, 12], leo_alt_m=500e3, final_orbits=10.0)
    objects = get(payload, "objects", Any[])
    for obj in objects
        obj["intercept"] = plan_intercept_for_object(
            obj;
            tof_days_grid=tof_days_grid,
            arrive_offset_hours=arrive_offset_hours,
            leo_alt_m=leo_alt_m,
            final_orbits=final_orbits,
        )
    end
    payload["intercept_note"] = Dict{String,Any}(
        "frame" => "heliocentric ecliptic (Sun μ)",
        "earth_model" => "Julia analytic J2000 Kepler model",
        "cost" => "minimize |v1−vE| + |vast−v2|",
        "dv_units" => "km/s",
        "time_scale" => "TDB JD",
        "capture_model" => "density-sphere-v1",
    )
    return payload
end

function parse_int_list(text)
    isempty(strip(text)) && return Int[]
    return [parse(Int, strip(part)) for part in split(text, ",") if !isempty(strip(part))]
end

function run_cli(args=ARGS)
    input = nothing
    output = nothing
    tof_min = 30
    tof_max = 180
    tof_step = 10
    arrive_hours = "-12,-6,0,6,12"
    leo_alt_m = 500e3
    final_orbits = 10.0
    idx = 1
    while idx <= length(args)
        arg = args[idx]
        if arg in ("--input", "--inp")
            idx += 1
            input = args[idx]
        elseif arg in ("--output", "--out")
            idx += 1
            output = args[idx]
        elseif arg == "--tof-min"
            idx += 1
            tof_min = parse(Int, args[idx])
        elseif arg == "--tof-max"
            idx += 1
            tof_max = parse(Int, args[idx])
        elseif arg == "--tof-step"
            idx += 1
            tof_step = parse(Int, args[idx])
        elseif arg == "--arrive-hours"
            idx += 1
            arrive_hours = args[idx]
        elseif arg == "--leo-alt-m"
            idx += 1
            leo_alt_m = parse(Float64, args[idx])
        elseif arg == "--final-orbits"
            idx += 1
            final_orbits = parse(Float64, args[idx])
        else
            throw(ArgumentError("unknown argument: $arg"))
        end
        idx += 1
    end

    source = input === nothing ? read(stdin, String) : read(input, String)
    payload = JSON.parse(source)
    tof_grid = collect(tof_min:tof_step:tof_max)
    offsets = parse_int_list(arrive_hours)
    attach_intercepts!(payload; tof_days_grid=tof_grid, arrive_offset_hours=offsets, leo_alt_m=leo_alt_m, final_orbits=final_orbits)
    encoded = JSON.json(payload, 2)
    if output === nothing
        print(encoded)
    else
        write(output, encoded)
    end
    return 0
end

export AU_M, DAY_S, MU_SUN, DEFAULT_ALBEDO, DEFAULT_DENSITY_KG_M3
export kepler_E_from_M, kepler_universal_propagate, lambert_universal
export diameter_from_H_m, capture_model, attach_intercepts!, run_cli

end
