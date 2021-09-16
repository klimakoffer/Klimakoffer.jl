struct Mesh
    nx::Int                     # Number of DOFs in x (Longitude)
    ny::Int                     # Number of DOFs in y (Latitude)
    dof::Int                    # Total number of DOFs: nx * ny
    h::Float64                  # Grid size [radians]: π/(ny-1)
    geom::Float64               # Geometrical parameter at poles
    csc2::Array{Float64,1}      # Metric term for the transformation to spherical coordinates: cosc²(θ), where θ is the colatitude angle
    cot::Array{Float64,1}       # Metric term for the transformation to spherical coordinates: cot(θ), where θ is the colatitude angle
    area::Array{Float64,1}      # Area of each cell on the surface of the sphere (*only* depends on latitude)
end

function Mesh(nx=128)
    ny = Int(nx/2+1) # 65
    dof = nx*ny

    h = pi / (ny-1) # Uniform grid size in radians

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

    return Mesh(nx,ny,dof,h,geom,csc2,cot,area)
end
