using Klimakoffer

Klimakoffer.images_to_maps(joinpath(@__DIR__, "..","examples", "test_instances"), joinpath(@__DIR__, "..","examples", "test_instances"), true, 1, 164, 97, 128, 10, 190, 210)

world = Klimakoffer.read_geography(joinpath(@__DIR__, "..",  "examples", "test_instances", "The_World_from_image128x65_1.dat" ), 128, 65)

ref = Klimakoffer.read_geography(joinpath(@__DIR__, "..","examples","test_instances", "testresult_image.dat"), 128, 65)

diff = ref-world

res = sum(diff)

if !isfile(joinpath(@__DIR__, "..",  "examples", "test_instances", "The_World_from_image128x65_1.dat" ))
    res +=1
end


