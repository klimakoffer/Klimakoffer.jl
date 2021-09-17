using Test
using Klimakoffer
using Printf

EXAMPLES_DIR = joinpath(pathof(Klimakoffer) |> dirname |> dirname, "examples")

@time @testset "Klimakoffer" begin
  @testset "equilibrium_temperature_1950.jl" begin
    @test_nowarn include(joinpath(EXAMPLES_DIR, "equilibrium_temperature_1950.jl"))

    @test isapprox(GlobTemp, 14.484963368768746, atol=1e-12)
  end
  @testset "equilibrium_temperature_co2.jl" begin
    @test_nowarn include(joinpath(EXAMPLES_DIR, "equilibrium_temperature_co2.jl"))

    @test isapprox(GlobTemp, 14.484963368768978, atol=1e-12)
  end
  @testset "equilibrium_temperature_co2_transient.jl" begin
    @test_nowarn include(joinpath(EXAMPLES_DIR, "equilibrium_temperature_co2_transient.jl"))

    @test isapprox(mean_temperature_yearly[end], 14.563278065969344, atol=1e-12)
  end
end
