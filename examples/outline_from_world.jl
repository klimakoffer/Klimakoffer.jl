using Klimakoffer

ref = Klimakoffer.read_geography(string("./examples/test_instances/world/ref_The_World_Outline128x65.dat"),128,65) #original outline from Klimakoffer

save_outline_from_path("./examples/test_instances/world/", "The_World128x65.dat", 128)

outline = Klimakoffer.read_geography("./examples/test_instances/world/The_World_Outline128x65.dat", 128, 65) # calculated outline from implemented function 

diff = ref - outline

res = sum(diff)

rm("./examples/test_instances/world/The_World_Outline128x65.dat")