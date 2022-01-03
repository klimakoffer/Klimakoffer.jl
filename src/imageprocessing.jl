using Images, ImageFiltering
using FileIO

# /Users/otzi/Desktop/Bachelorarbeit/world5400x2700.jpg

# TODO: maybe we can classify some day by rgb values instead of grayscale values from 0 to 255
"""
    images_to_maps(...)

* has nearly the same inputs as convert_image_to_world but there is a directorypath instead of a filepath of one image and
* a targetpath, where you save the 12 maps 
"""
function images_to_maps(dirsourcepath="/Users/otzi/Desktop/Bachelorarbeit/NASA_Bilder/",targetpath = "./input/world/monthly_maps/",blurred=1, imglong=5400, imglat=2700, nlongitude=128, limit_ocean = 10, limit_land= 140, limit_seaice=205 )
    
    if isdir(dirsourcepath)
        for (root, dirs, files) in walkdir(dirsourcepath)
            
            println(root)
            println(files)
            for file in files
                if occursin("world", string(file))
                    println(file)
                    filepath = string(root, file)
                    println(filepath)
                    world = convert_image_to_world(filepath, blurred, imglong, imglat, nlongitude, limit_ocean, limit_land, limit_seaice)
                    save_world_by_month(world, string(file), targetpath)
                end
            end
        end
    end
end 
"""
    convert_image_to_world(...)

* imagefile is the filapath of that lanscape we want to convert into our map with dicrete values for each landcover
* imglong and imglat are the longitude and latitude of the image, so the pixels ratio
* dirpath is the directory path where the new map should be save under the new filename
* nlongitude is the preferred longitude for the nearest neighbor interpolation
* limit_[...] is the limiting value for the classification from grayscale per landcover, 
*       like 0 <= ocean <= limit_ocean <= land <= limit_land <= ice <= limit_seaice <= snow <= 255
"""

# TODO: gaussian blur einbinden 
function convert_image_to_world(imagefile, blurred=1, imglong=5400, imglat=2700, nlongitude=128, limit_ocean = 10, limit_land= 190, limit_seaice=210)

    world = FileIO.load(imagefile)
    world = transpose(world) # without to transpose the NASA images, we would create a symmetrically mirrored map
    (longitude, latitude) = size(world)

    # from RGB values to grayscale (better options to classify the land cover)
    grayworld = Gray.(world)
    grayworld = real.(grayworld) .*255 # Julia's gray values multiplied by 255 to get common grayscale from 0 (black) to 255 (white)

    if blurred!=1
        grayworld = imfilter(grayworld, Kernel.gaussian(blurred))
    end

    for lat = 1:latitude
        for long = 1:longitude
            if grayworld[long, lat] > 0 && grayworld[long, lat] <= limit_ocean # ocean
                grayworld[long, lat] = 5
            elseif grayworld[long, lat] > limit_ocean && grayworld[long, lat] <= limit_land # land
                grayworld[long, lat] = 1
            elseif grayworld[long, lat] > limit_land && grayworld[long, lat] <= limit_seaice # sea ice
                grayworld[long, lat] = 2
            elseif grayworld[long, lat] > limit_seaice && grayworld[long, lat] <= 255 # snow cover (snow has the highest albedo value)
                grayworld[long, lat] = 3
            end
        end
    end
   
    return nn_interpolation(grayworld, nlongitude)
end

function save_world_by_month(array_in, img_filename, dirpath)
    if !isdir(dirpath)
        mkdir(dirpath)
    end 
    
    # Because our simulation begins at 21. march, we have to label the march with 1 
    months = Dict([("jan", 11), ("feb", 12), ("mar", 1), ("apr", 2), ("may", 3), ("jun", 4), ("jul", 5), ("aug", 6), ("sep", 7), ("oct", 8), ("nov", 9), ("dec", 10)])
    month = 0
    for (key, val) in months
        if occursin(string(key), img_filename)
            (scaledlong, scaledlat) = size(array_in)
            
            new_filename = string("The_World_from_image", scaledlong, "x", scaledlat,"_", val, ".dat") 
            new_filepath = string(dirpath, new_filename)
            if !isfile(new_filepath)
                touch(new_filepath) 
            end 

            open(new_filepath,"w") do file 
                for lat = 1:scaledlat
                    for long = 1:scaledlong
                        write(file,string(array_in[long, lat]))
                    end 
                    write(file, "\n")
                end
            end
        end
    end
end

function test_gaussian_blur(imagefile="/Users/otzi/Desktop/Bachelorarbeit/world5400x2700.jpg")
    world = FileIO.load(imagefile)
    blurred_world = imfilter(world, Kernel.gaussian(10))
    save(File{format"JPEG"}("/Users/otzi/Desktop/Bachelorarbeit/blurredworld.jpeg"), blurred_world)
end