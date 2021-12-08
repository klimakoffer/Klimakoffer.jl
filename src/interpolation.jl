
function upscale_inputs(nlongitude=256)
     
     nlatitude = convert(Int64,nlongitude/2+1)

     println(string("Erstelle Maps für die Auflösung: ", nlongitude, "x",nlatitude))
     upscale_albedo(nlongitude)
     upscale_map(nlongitude)
     get_outline_from_map(nlongitude)
end 

"""
    nn_interpolation(...)

* array_in is the array/map, which will be interpolated in a higher resolution.
* nlongitude is the new width of the mesh, nlatitude will be calculated as nlongitude/2 + 1
"""
function nn_interpolation(array_in, nlongitude=256)
    
     (width_in, height_in) = size(array_in)
     
     width_out = nlongitude
     height_out = convert(Int64, width_out/2 + 1)
     array_out = fill(0, (width_out, height_out))
   
     # calculate the new distance between each dot/pixle (ratio)
     # "-1" because the indicies in Julia
     ratio_w = (width_in-1)/(width_out-1)
     ratio_h = (height_in-1)/(height_out-1)

     row = ones(Float64, width_out)
     for i = 1 : length(row)         
          r = row[i] + (i-1) * ratio_w
          row[i] = r
     end
     row[end] = round(width_in) 
    
     col = ones(Float64, height_out)
     for i = 1 : length(col)   
          c = col[i] + (i-1) * ratio_h
          col[i] = c
     end
     col[end] = round(height_in) 

     width_positions = round.(Int64, row)  
     col .+= 0.0000001  # Korekktur fürs Runden, da round(2.5) sonst 2 ergibt ???!!!
     height_positions = round.(Int64, col)            
     
     wide = 1
     high = 1
     
     for b in width_positions  
          high = 1
          for h in height_positions
               array_out[wide, high] = array_in[b, h]
               high += 1
               
          end
          wide += 1
     end
     return array_out
end

function upscale_map(nlongitude=256) 

     world = Klimakoffer.read_geography(string("./input/", 128, "x", 65,"/The_World.dat"), 128, 65)
     (old_long, old_lat) = size(world)
     upscaled_map = nn_interpolation(world, nlongitude)
     (nlongitude, nlatitude) = size(upscaled_map)

     if !isdir(string("./input/", nlongitude, "x", nlatitude))
          mkdir(string("./input/", nlongitude, "x", nlatitude))
     end
     if !isfile(string("./input/", nlongitude, "x", nlatitude,"/The_World.dat"))
          touch(string("./input/", nlongitude, "x", nlatitude,"/The_World.dat")) #NO: evtl nach anderer "create file"-Funktion umschauen
     end 

     open(string("./input/", nlongitude, "x", nlatitude,"/The_World.dat"),"w") do file 
          for lat = 1:nlatitude
               for long = 1:nlongitude
                    write(file, string(upscaled_map[long,lat]))
               end
               write(file, "\n")
           end
     end

     println(string("Upscale The_World.dat from ", old_long, "x", old_lat," to ", nlongitude,"x", nlatitude))
end

"""
    bilinear_interpolation(...)

* array_in is the (albedo-)array, which will be interpolated in a higher resolution.
* nlongitude is the new width of the mesh, nlatitude will be calculated as nlongitude/2 + 1
"""
function bilinear_interpolation(array_in, nlongitude =256)

     (width_in, height_in) = size(array_in)

     width_out = nlongitude
     height_out = convert(Int64, (width_out/2 + 1)) 
     array_out = Array{Float64}(ones((width_out, height_out)))

     # calculate the new distance between each dot/pixle (ratio)
     # "-1" because the indicies in Julia
     ratio_w = (width_in-1)/(width_out-1)
     ratio_h = (height_in-1)/(height_out-1)

     row = ones(Float64, width_out)
     col = ones(Float64, height_out)

     high = 1
     wide = 1

     for b = 1 : width_out 
          high = 1 
          for h = 1 : height_out

               # the coordinates from the current pixels/knot
               y = 1 + (h-1) * ratio_h  
               x = 1 + (b-1) * ratio_w  
               
               # calculate coordinates from surrounding pixels
               x_floor = floor(Int64, x)                    
               y_floor = floor(Int64, y)                        
               x_ceil = min(ceil(Int64, x), width_in)              
               y_ceil = min(ceil(Int64, y), height_in)            

               # get values of the surrounding pixels 
               # case 1: x and y are intergers => x_ceil = x_floor, y_ceil = y_floor     
               if x_ceil == x_floor && y_ceil == y_floor 
                    val = array_in[x_ceil , y_ceil]

               # case 2: x is an integer 
               elseif x_ceil == x_floor
                    val1 = array_in[x_ceil, y_floor]
                    val2 = array_in[x_ceil, y_ceil]
                    val = val1 * (y_ceil - y) + val2 * (y- y_floor)

               #case 3: y is an integer
               elseif y_ceil == y_floor
                    val1 = array_in[x_floor, y_ceil]
                    val2 = array_in[x_ceil, y_ceil]
                    val = val1 * (x_ceil - x) + val2 * (x - x_floor)
               # case 4: x and y are floating point numbers     
               else 
                    v1 = array_in[x_floor, y_floor]
                    v2 = array_in[x_ceil, y_floor]
                    v3 = array_in[x_floor, y_ceil]
                    v4 = array_in[x_ceil, y_ceil]

                    # calculate the pixels (x, y) value
                    val1 = (v1 * (x_ceil - x) + v2 * (x - x_floor))  # first linear interpolation in x direction between v1 and v2
                    val2 = (v3 * (x_ceil - x) + v4 * (x - x_floor))  # second linear interpolation in x direction between v3 and v4
                    val = val1 * (y_ceil - y) + val2 * (y- y_floor)  # final "bilinear interpolation"
               end 

               array_out[wide, high] = round(val, digits=2)
               high += 1
          end 
          wide += 1
     end 
     return array_out
end 

"""
    upscale_albedo(...)

* scales up the albedo map 
"""
function upscale_albedo(nlongitude=256)

     albedo = Klimakoffer.read_albedo(string("./input/", 128, "x", 65,"/albedo.dat"), 128, 65)
     (oldlong, oldlat) = size(albedo)
     upscaled_albedo = bilinear_interpolation(albedo, nlongitude)
     (nlongitude, nlatitude) = size(upscaled_albedo)

     if !isdir(string("./input/", nlongitude, "x", nlatitude))
          mkdir(string("./input/", nlongitude, "x", nlatitude))
     end
     if !isfile(string("./input/", nlongitude, "x", nlatitude,"/albedo.dat"))
          touch(string("./input/", nlongitude, "x", nlatitude,"/albedo.dat"))
     end 

     open(string("./input/", nlongitude, "x", nlatitude,"/albedo.dat"),"w") do file 
          for lat = 1:nlatitude
               for long = 1:nlongitude
                    write(file, string("      ",upscaled_albedo[long,lat]))
               end
               write(file,"\n")
           end
     end
     println(string("Upscale albedo.dat from ", oldlong, "x", oldlat," to ", nlongitude,"x", nlatitude))
end

"""
    get_outline_from_map(...)

* get the outline from an existing map, so the coast is presented by ones 
"""
function get_outline_from_map(nlongitude=128)

     nlatitude = convert(Int64, nlongitude/2+1)
     world = Klimakoffer.read_geography(string("./input/", nlongitude, "x", nlatitude,"/The_World.dat"), nlongitude, nlatitude)
     world = clear_map(world)
     outline = Array{Int64}(zeros((nlongitude, nlatitude)))

     for b = 1:nlongitude
          for h = 1:nlatitude
               if world[b,h] == 1 || world[b,h] == 3 # 1 = land, 3 = ice

                    if b == 1 && h == 1 
                         if !(world[b+1,h] == 1 && world[b,h+1] == 1)
                              outline[b,h] = 1
                         end

                    elseif b == 1 && h != 1 && h != nlatitude
                         if world[b,h+1] == 0 || world[b,h-1] == 0 || world[b+1,h] == 0
                              outline[b,h] = 1
                         end

                    elseif b == 1 && h == nlatitude 
                         if !(world[b+1,h] == 1 && world[b,h-1] == 1) 
                              outline[b,h] = 1
                         end 

                    elseif h == nlatitude && b != 1 && b != nlongitude
                         if world[b-1,h] == 0 || world[b, h-1] == 0 || world[b+1,h] == 0
                              outline[b,h] = 1 
                         end 

                    elseif h == nlatitude && b == nlongitude
                         if !(world[b-1,h] == 1 && world[b,h-1] == 1)
                              outline[b,h] = 1
                         end 

                    elseif b == nlongitude && h != 1 && h != nlatitude 
                         if  world[b-1,h] == 0 || world[b,h-1] == 0 || world[b,h+1] == 0
                              outline[b,h] = 1
                         end 

                    elseif b == nlongitude && h == 1
                         if !(world[b-1,h] == 1 && world[b,h+1] == 1)
                              outline[b,h] = 1
                         end 

                    elseif h == 1 && b != nlongitude && b != 1
                         if world[b-1,h] == 0 || world[b,h+1] == 0 || world[b+1,h] == 0
                              outline[b,h] = 1
                         end 

                    else 
                         if world[b-1,h] == 0 || world[b+1,h] == 0 || world[b,h-1] == 0 || world[b,h+1] == 0
                              outline[b,h] = 1
                         end 

                    end
               end
          end
     end

     if !isdir(string("./input/", nlongitude, "x", nlatitude))
          mkdir(string("./input/", nlongitude, "x", nlatitude))
     end
     if !isfile(string("./input/", nlongitude, "x", nlatitude,"/The_World_Outline.dat"))
          touch(string("./input/", nlongitude, "x", nlatitude,"/The_World_Outline.dat"))
     end 

     open(string("./input/", nlongitude, "x", nlatitude,"/The_World_Outline.dat"),"w") do file 
          for lat = 1:nlatitude
               for long = 1:nlongitude
                    write(file, string(outline[long,lat]))
               end
               write(file, "\n")
           end
     end
end

"""
    clear_map(map)

* a sidefunction, which presents land (1) and permanent snow cover (3) with ones and 
  perennial sea ice (2), lakes, inland seas (4) and ocean (5) with zeros 
"""
function clear_map(map)
     (nlongitude, nlatitude) = size(map)
     for b =1:nlongitude
          for h=1:nlatitude
               if map[b,h] == 3 
                    map[b,h] = 1
                    
               elseif map[b,h] == 2 || map[b,h] >= 4
                    map[b,h] = 0
               end
          end 
     end 
     return map
end 

