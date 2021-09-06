using Test
using Klimakoffer

@time @testset "Klimakoffer" begin
  @testset "Aswer to the Ultimate Question" begin
    @test answer() == 42
  end
end
