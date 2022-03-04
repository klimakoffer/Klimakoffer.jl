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

  test_file = "equilibrium_temperature_1950_150x75.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(GlobTemp, 14.967892384136475, atol=1e-10)
  end

  test_file = "land_coverage_from_scaled_map.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(landcover, 23.341346153846153, atol=1.5) # aktuell ist es der Wert der originalen 128x65 Karte aus dem Klimakoffer
  end

  test_file = "outline_from_world.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(res, 0, atol=0) 
  end

  test_file = "image_to_world.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(result, 0, atol=0.5) 
  end

  test_file = "clear_map.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(res, 0, atol=0) 
  end

  test_file = "upscale_world.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(res, 0, atol=0) 
  end

  test_file = "upscale_albedo.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(res, 0, atol=0) 
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

    @test isapprox(GlobTemp, 14.490442457180537, atol=1e-12)
  end

  test_file = "transient_temperature_co2.txt"
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(sol.mean_temperature_yearly[end], 15.124553126197249, atol=1e-12)
  end

  test_file = "equilibrium_temperature_sea_ice_extent.jl" 
  @testset "$test_file" begin
    println("")
    println("Running ",test_file)
    println("")
    @test_nowarn include(joinpath(EXAMPLES_DIR, test_file))

    @test isapprox(GlobTemp, -5.4260487452118795, atol=1e-12)
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
