using Test
using Klimakoffer

@time @testset "Klimakoffer" begin
  @testset "global equilibrium temperature 1950" begin
    @test isapprox(main().GlobTemp, 14.484963368768746, atol=1e-12)
  end
end
