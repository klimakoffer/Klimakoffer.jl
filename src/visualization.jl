

function showGif(field,title_in="Solar Forcing of the Earth: t = ",minval=-46.073107944313556,maxval=32.52216883281448)
    World = read_geography(joinpath(@__DIR__,"..","The_World_Outline.dat"))*(-80)
    Plots.@gif for t in 1:size(field,3)
        title = title_in * lpad(string(floor(t*365/48)), 3, ' ') * " days"
    
        data = field[:,:,t]
        data = transpose(data)
        Plots.heatmap(LinRange(-180,180,128),LinRange(-90,90,65),data,
            clims=(minval,maxval),
            yflip=true,
            title=title,
            xlabel="longitude [°]",ylabel="latitude [°]")
        Plots.contour!(LinRange(-180,180,128),LinRange(-90,90,65),transpose(World),xlim=(-180,180), ylim=(-90, 90),levels=-40:40,seriesalpha=0.05)
        
        fpath = "png/field/" * lpad(string(t), 3, '0') * ".png"
        println(fpath)
        Plots.savefig(fpath)
    end
end
