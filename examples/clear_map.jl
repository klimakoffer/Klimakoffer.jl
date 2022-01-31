using Klimakoffer

world = Klimakoffer.read_geography(joinpath(@__DIR__,"..","input","world","The_World128x65.dat"),128,65)

cleared_map = Klimakoffer.clear_map(world)

clear = Klimakoffer.read_geography(joinpath(@__DIR__,"..","examples","reference_data","clear_map128x65.dat"),128,65)

diff = clear - cleared_map

res = sum(abs.(diff))