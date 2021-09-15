

function showGif(field,title_in="Solar Forcing of the Earth: t = ",minval=-46.073107944313556,maxval=32.52216883281448)
    @gif for t in 1:size(field,3)
        title = title_in * lpad(string(floor(t*365/48)), 3, ' ') * " days"
    
        data = field[:,:,t]
        data = transpose(data)
        heatmap(LinRange(-180,180,128),LinRange(-90,90,65),data,
            clims=(minval,maxval),
            yflip=true,
            title=title,
            xlabel="longitude [°]",ylabel="latitude [°]")
    
        fpath = "png/field/" * lpad(string(t), 3, '0') * ".png"
        println(fpath)
        savefig(fpath)
    end
end