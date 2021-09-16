struct Mesh
    nx::Int         # Number of DOFs in x (Longitude)
    ny::Int         # Number of DOFs in y (Latitude)
    dof::Int
    h::Float64
    sh2::Float64
    geom::Float64
    csc2::Array{Float64,1}
    cot::Array{Float64,1}
    area::Array{Float64,1}
end

function Mesh(nx::Int)
    ny = Int(nx/2+1) # 65
    dof = nx*ny

    h = pi / (ny-1) # Uniform grid size in radians
    sh2 = 1 / h^2

    # Trigonometric definitions for inner degrees of freedom
    csc2 = zeros(Float64,ny) # csc^2(theta), where theta is the colatitude angle (from the pole)
    cot  = zeros(Float64,ny) #   cot(theta), where theta is the colatitude angle (from the pole)
    for j=2:ny-1
        theta = h*(j-1)
        sintheta = sin(theta)
        csc2[j] = 1/sintheta^2
        cot[j] = cos(theta)/sintheta
    end

    area = zeros(Float64,ny) # Fractional area for the cells depending on latitude
    # Grid point fractional area for inner DOFs
    for j=2:ny-1
        area[j] =  sin(0.5*h)*sin(h*(j-1))/nx
    end

    # Fractional area for the poles
    area[1]  = 0.5*(1 - cos(0.5*h))
    area[ny] = area[1]

    geom = sin(0.5*h)/area[1]

    return Mesh(nx,ny,dof,h,sh2,geom,csc2,cot,area)
end
