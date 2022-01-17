using Klimakoffer

imagefile = "./examples/test_instances/test_world_mar.jpg"

Klimakoffer.images_to_maps("./examples/test_instances/", "./examples/test_instances/", true, 1, 164, 97, 128, 10, 190, 210)

world = Klimakoffer.read_geography("./examples/test_instances/The_World_from_image128x65_1.dat", 128, 65)

ref = Klimakoffer.read_geography("./examples/test_instances/testresult_image.dat", 128, 65)

diff = ref-world

res = sum(diff)

if !isfile("./examples/test_instances/The_World_from_image128x65_1.dat")
    res +=1
end


