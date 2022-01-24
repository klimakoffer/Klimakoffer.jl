using Klimakoffer

Klimakoffer.images_to_maps(joinpath(@__DIR__,"examples", "test_instances"), joinpath(@__DIR__,"examples", "test_instances"), true, 1, 10, 11, 10, 10, 190, 210)

world = Klimakoffer.read_geography(joinpath(@__DIR__,  "examples", "test_instances", "The_World_from_image10x6_3.dat" ), 10, 6)

ref = Klimakoffer.read_geography(joinpath(@__DIR__,"examples","test_instances", "test_image_result.dat"), 10, 6)

diff = ref-world

res = sum(abs.(diff))

if !isfile(joinpath(@__DIR__, "..",  "examples", "test_instances", "The_World_from_image128x65_1.dat" ))
    res +=1
end

result = res/(10*6)

rm(joinpath(@__DIR__, "..",  "examples", "test_instances", "The_World_from_image10x6_3.dat" ))

