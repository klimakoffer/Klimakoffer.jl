using Klimakoffer

Klimakoffer.upscale_world("./examples/test_instances/world/", "The_World128x65.dat",150)

upscaled = Klimakoffer.read_geography("./examples/test_instances/world/The_World150x76.dat", 150, 76)

ref = Klimakoffer.read_geography("./examples/test_instances/world/ref_The_World150x76.dat", 150, 76)

diff = ref-upscaled

res = sum(diff)

rm("./examples/test_instances/world/The_World150x76.dat")