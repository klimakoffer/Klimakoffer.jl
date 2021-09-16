using Test
using Klimakoffer

EXAMPLES_DIR = joinpath(pathof(Klimakoffer) |> dirname |> dirname, "examples")

@time @testset "Klimakoffer" begin
  @testset "equilibrium_temperature_1950.jl" begin
    @test_nowarn include(joinpath(EXAMPLES_DIR, "equilibrium_temperature_1950.jl"))

    @test isapprox(GlobTemp, 14.484963368768746, atol=1e-12)
  end
end
