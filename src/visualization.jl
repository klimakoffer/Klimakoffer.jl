using Interpolations;

function showGif(field,title_in="Solar Forcing of the Earth: t = ",minval=-46.073107944313556,maxval=32.52216883281448)
    World = read_geography(joinpath(@__DIR__,"..", "input", "The_World_Outline.dat"))*(-80)
    (width, height) = size(World) # NO hinzugefügt
    Plots.@gif for t in 1:size(field,3)
        title = title_in * lpad(string(floor(t*365/48)), 3, ' ') * " days"
    
        data = field[:,:,t]
        data = transpose(data)
        Plots.heatmap(LinRange(-180,180,width),LinRange(-90,90,height),data,
            clims=(minval,maxval),
            yflip=true,
            title=title,
            xlabel="longitude [°]",ylabel="latitude [°]")
        Plots.contour!(LinRange(-180,180,width),LinRange(-90,90,height),transpose(World),xlim=(-180,180), ylim=(-90, 90),levels=-40:40,seriesalpha=0.05)        
        fpath = "png/field/" * lpad(string(t), 3, '0') * ".png"
        println(fpath)
        Plots.savefig(fpath)
    end
end

# This is a dirty function. What was it again?
# First make it work, then make it pretty ... something like that.(?)
function apply_richardson_projection(data)

    lat = [0 , 5 , 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90] * pi/180.0
    X   = [1.0000, 0.9986, 0.9954, 0.99  , 0.9822, 0.973 , 0.96  , 0.9427, 0.9216, 0.8962, 0.8679, 0.835 , 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322]
    Y   = [0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958, 0.5571, 0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000]

    itpl_X = LinearInterpolation(lat,X);
    itpl_Y = LinearInterpolation(lat,Y);

    calc_x = (radius,longitude,latitude) -> 0.8487 * radius .* itpl_X.(abs.(latitude)) .* longitude;
    calc_y = (radius,longitude,latitude) -> 1.3523 * radius .* itpl_Y.(abs.(latitude)) .* 0.2536 .* sign.(latitude)

    nx = size(data,1)
    ny = size(data,2)

    longitude = LinRange(-pi,pi,nx) .* ones(ny)'
    latitude  = ones(nx) .* LinRange(-pi/2,pi/2,ny)'

    x = calc_x(1.0,longitude,latitude);
    y = calc_y(1.0,longitude,latitude);
    z = reverse(reverse(data),dims=1)

    x = 180.0*(x ./ maximum(x));
    y =  90.0*(y ./ maximum(y));

    return x,y,z

end

