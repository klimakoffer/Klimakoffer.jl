using Klimakoffer

Klimakoffer.upscale_albedo("./examples/test_instances/albedo/", "albedo128x65.dat",150)

upscaled = Klimakoffer.read_albedo("./examples/test_instances/albedo/albedo150x76.dat", 150, 76)

ref = Klimakoffer.read_albedo("./examples/test_instances/albedo/ref_albedo150x76.dat", 150, 76)

diff = ref-upscaled

res = sum(diff)

rm("./examples/test_instances/albedo/albedo150x76.dat")