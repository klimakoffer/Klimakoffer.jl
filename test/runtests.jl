using Test
using Klimakoffer
using Printf

EXAMPLES_DIR = joinpath(pathof(Klimakoffer) |> dirname |> dirname, "examples")

@time @testset "Klimakoffer" begin
  test_file = "equilibrium_temperature_1950.jl"
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(GlobTemp, 14.484963368768746, atol=1e-12)
  end
  test_file = "equilibrium_temperature_co2.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(GlobTemp, 14.484963368768978, atol=1e-12)
  end
  test_file = "equilibrium_temperature_co2_transient.jl"
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(sol.mean_temperature_yearly[end], 14.563278065969344, atol=1e-12)
  end
end
