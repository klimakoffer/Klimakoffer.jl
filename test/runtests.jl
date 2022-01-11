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

    @test isapprox(GlobTemp, 14.484963368770806, atol=1e-12)
  end

  test_file = "equilibrium_temperature_albedo.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(GlobTemp, 15.007699256045628, atol=1e-12)
  end

  test_file = "equilibrium_temperature_co2.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(GlobTemp, 14.495916499277554, atol=1e-12)
  end

  test_file = "transient_temperature_co2.txt"
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(sol.mean_temperature_yearly[end], 15.12471138630498, atol=1e-12)
  end

  @testset "Printing types to the REPL" begin
    mesh = Mesh()
    @test_nowarn show(stdout, mesh)
    println(stdout)

    model = Model(mesh, 48)
    @test_nowarn show(stdout, model)
    println(stdout)

    discretization = Discretization(mesh, model, 48)
    @test_nowarn show(stdout, discretization)
    println(stdout)
  end
end
