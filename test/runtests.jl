using AsteroidPlannerBackend
using Test
using LinearAlgebra

@testset "orbital mechanics" begin
    M = 2.1
    e = 0.63
    E = kepler_E_from_M(M, e)
    @test E - e * sin(E) ≈ M atol = 1e-11

    position = [AU_M, 0.0, 0.0]
    velocity = [0.0, sqrt(MU_SUN / AU_M), 0.0]
    period = 2.0 * pi * sqrt(AU_M^3 / MU_SUN)
    propagated, _ = kepler_universal_propagate(position, velocity, period)
    @test norm(propagated - position) < 1_000.0

    departure = [AU_M, 0.0, 0.0]
    arrival = [0.0, 1.524 * AU_M, 0.0]
    tof = 200.0 * DAY_S
    departure_velocity, _ = lambert_universal(departure, arrival, tof)
    endpoint, _ = kepler_universal_propagate(departure, departure_velocity, tof)
    @test norm(endpoint - arrival) < 25_000.0
    @test_throws ArgumentError lambert_universal(ones(3), [1.0, 2.0, 3.0], 0.0)
end

@testset "capture model" begin
    diameter_m = diameter_from_H_m(22.0, 0.14)
    @test diameter_m > 0.0
    @test diameter_m ≈ 1329.0 / sqrt(0.14) * 10.0^(-22.0 / 5.0) * 1000.0

    elements = Dict(
        "epoch_jd_tdb" => 2.461e6,
        "a_AU" => 1.2,
        "e" => 0.1,
        "i_deg" => 5.0,
        "raan_deg" => 50.0,
        "argp_deg" => 10.0,
        "ma_deg" => 20.0,
    )
    obj = Dict(
        "next_ca" => Dict("jd_tdb" => 2.4612e6, "diameter_km" => 10.0),
        "elements" => elements,
    )
    capture = capture_model(obj, elements, 2.4612e6, 0.0, 5_000.0; final_orbits=4.0)
    @test capture["model"] == "density-sphere-v1"
    @test capture["radius_m"] ≈ 5_000.0
    @test capture["mass_kg"] ≈ 2_000.0 * (4.0 / 3.0) * pi * 5_000.0^3
    @test capture["mu_m3_s2"] ≈ 6.67430e-11 * capture["mass_kg"]
    @test capture["circular_speed_mps"] > 0.0
    @test capture["escape_speed_mps"] > capture["circular_speed_mps"]
    @test capture["capture_dv_kms"] > 0.0
    @test capture["propagation_max_radius_m"] >= capture["target_orbit_radius_m"]
    @test length(capture["orbit_polyline_relative_m"]) == 241
    @test length(capture["orbit_polyline_heliocentric_au"]) == 241
    @test capture["requested_final_orbits"] == 4.0
    @test capture["propagated_final_orbits"] ≈ 4.0
    @test capture["propagation_days"] ≈ 4.0 * capture["orbit_period_s"] / DAY_S

    faster_capture = capture_model(obj, elements, 2.4612e6, 1.0, 5_000.0; final_orbits=4.0)
    @test faster_capture["capture_dv_kms"] > capture["capture_dv_kms"]

    missing = capture_model(Dict("next_ca" => Dict()), elements, 2.4612e6, 1.0, 5_000.0; final_orbits=4.0)
    @test missing["stable_final_orbit"] == false
    @test "missing physical size" in missing["stability_flags"]
end

@testset "payload integration" begin
    payload = Dict(
        "schema" => "asteroid_intercept_planner/2.0",
        "objects" => Any[
            Dict(
                "des" => "TEST-1",
                "next_ca" => Dict(
                    "jd_tdb" => 2_461_200.5,
                    "cd_tdb" => "2026-Jul-09 00:00",
                    "dist_au" => 0.0123,
                    "diameter_km" => 1.0,
                ),
                "elements" => Dict(
                    "epoch_jd_tdb" => 2_461_000.5,
                    "a_AU" => 1.2,
                    "e" => 0.1,
                    "i_deg" => 5.0,
                    "raan_deg" => 50.0,
                    "argp_deg" => 10.0,
                    "ma_deg" => 20.0,
                ),
            ),
        ],
    )
    attach_intercepts!(
        payload;
        tof_days_grid=[120],
        arrive_offset_hours=[0],
        leo_alt_m=500e3,
        final_orbits=3.0,
    )
    intercept = payload["objects"][1]["intercept"]
    @test intercept !== nothing
    @test haskey(intercept, "depart_jd_tdb")
    @test haskey(intercept, "lambert_polyline_xyz_au")
    @test haskey(intercept, "capture")
    @test haskey(intercept["capture"], "capture_dv_kms")
    @test intercept["capture"]["requested_final_orbits"] == 3.0
end
