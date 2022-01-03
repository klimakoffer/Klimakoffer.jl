using Klimakoffer

world = Klimakoffer.read_geography("./input/world/The_World128x65.dat", 128, 65)

rand_long = rand(128:2048)

upscaled_world = nn_interpolation(world, rand_long)

(nlong,nlat) = size(upscaled_world)

landcover = 100*count(i->(i==1), upscaled_world)/(nlong*nlat)