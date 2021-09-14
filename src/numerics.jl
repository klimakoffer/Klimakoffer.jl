using UnPack

#Mapping functions
##################
function index2d(k,nx)
    j = floor(Int,(k-1)/nx)+1
    # TODO: add consistency check with ny?
    return k-(j-1)*nx,j 
end

function index1d(i,j,nx)
    # TODO: add consistency check with ny?
    return i+(j-1)*nx
end

struct mesh_t
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

function mesh_init(nx::Int)
    ny = Int(nx/2+1) # 65
    dof = nx*ny

    h = pi / (ny-1) # Uniform grid size in radians
    sh2 = 1 / h^2

    # Trigonometric definitions for inner degrees of freedom
    csc2 = zeros(Float64,ny) # csc^2(theta), where theta is the colatitude angle (from the pole)
    cot  = zeros(Float64,ny) #   cot(theta), where theta is the colatitude angle (from the pole)
    for j=2:ny-1
        theta = h/(j-1)
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

    return mesh_t(nx,ny,dof,h,sh2,geom,csc2,cot,area)
end

struct model_t
    D_DiffCoeff::Array{Float64,2}
    C_HeatCapacity::Array{Float64,2}
    a_albedo::Array{Float64,2}
    SolarForcing::Array{Float64,3}
    A_coeff::Float64
    B_coeff::Float64
end

function model_init(mesh,NT)
    @unpack nx,ny = mesh

    D_DiffCoeff    = ones(Float64,nx,ny)
    C_HeatCapacity = ones(Float64,nx,ny)
    a_albedo       = ones(Float64,nx,ny)
    SolarForcing   = ones(Float64,nx,ny,NT)

    # Call fancy routines from the puppies team 
    A_coeff = 210.3  #[W/m^2]   : CO2 coefficient in the notes/paper
    B_coeff = 2.15   #[W/m^2/°C]: sensitivity of the seasonal cycle and annual change in the forcing agents

    return model_t(D_DiffCoeff,C_HeatCapacity,a_albedo,SolarForcing,A_coeff,B_coeff)
end

function main()
    #Input
    ######
    nx = 16
    NT = 48 # Number of time-steps per year

    #Initial definitions
    ####################

    #Time-discretization parameters
    dt = 1/NT

    #Mesh construction
    mesh = mesh_init(nx)

    #Static parameters
    

    #Solver variables 
    Temp = zeros(Float64,mesh.dof)
    Temp.= 5 # Magic initialization

    #Read in parameters
    ###################

    model = model_init(mesh,NT)


    # Assemble the matrix
    #####################

    A = ComputeMatrix(mesh,NT,model)

    RHS     = zeros(Float64,mesh.dof)
    LastRHS = zeros(Float64,mesh.dof)

    time_step = 2
    UpdateRHS!(RHS, mesh, NT, time_step, Temp, model, LastRHS)

    return (;A,RHS)
end

function ComputeMatrix(mesh,NT,model)
    @unpack nx,ny,dof,h,sh2,geom,csc2,cot,area = mesh
    @unpack D_DiffCoeff,C_HeatCapacity,B_coeff = model
    A    = zeros(Float64,dof,dof)
    # Inner DOFs (c coefficients are divided by h^2)
    for j=2:ny-1
        for i=1:nx

            # Compute coefficients
            c0 = 2 * sh2 * D_DiffCoeff[i,j] * (1 + csc2[j]) + 2 * C_HeatCapacity[i,j] * NT + B_coeff

            # Get diffusion coefficient 
            if (i == 1) # Periodic BC
                dPhi = (D_DiffCoeff[2,j] - D_DiffCoeff[nx,j]) / (2 * h)
            elseif (i == nx) # Periodic BC
                dPhi = (D_DiffCoeff[1,j] - D_DiffCoeff[nx-1,j]) / (2 * h)
            else # Inner DOFs
                dPhi = (D_DiffCoeff[i+1,j] - D_DiffCoeff[i-1,j]) / (2 * h)
            end

            c1 = sh2 * csc2[j] * (D_DiffCoeff[i,j] - 0.5 * h * dPhi)
            c3 = sh2 * csc2[j] * (D_DiffCoeff[i,j] + 0.5 * h * dPhi)

            dTheta= (D_DiffCoeff[i,j+1] - D_DiffCoeff[i,j-1]) / (2 * h)
            c2 = sh2 * (D_DiffCoeff[i,j] - 0.5 * h * (D_DiffCoeff[i,j] * cot[j] + dTheta))
            c4 = sh2 * (D_DiffCoeff[i,j] + 0.5 * h * (D_DiffCoeff[i,j] * cot[j] + dTheta))

            # Fill matrix A
            row_idx = index1d(i,j,nx)
            col_idx_c0 = row_idx
            col_idx_c1 = row_idx - 1
            col_idx_c3 = row_idx + 1
            col_idx_c2 = row_idx - nx
            col_idx_c4 = row_idx + nx

            A[row_idx,col_idx_c0] = c0
            A[row_idx,col_idx_c1] = -c1
            A[row_idx,col_idx_c2] = -c2
            A[row_idx,col_idx_c3] = -c3
            A[row_idx,col_idx_c4] = -c4

            # Get RHS
            # RHS[row_idx] = -(c0 * Temp[col_idx_c0] - c1 * Temp[col_idx_c0] - c2 * Temp[col_idx_c2] - c3 * Temp[col_idx_c3] - c4 * Temp[col_idx_c4]) + 4 * C_HeatCapacity[i,j]
        end
    end

    # Poles
    Tarea = area[1] + area[2]

    DT2np = zeros(Float64,nx)
    DT2sp = zeros(Float64,nx)
    for i=1:nx
        DT2np[i] = (area[1] * D_DiffCoeff[1,1]  + area[2] * D_DiffCoeff[i,2]) / Tarea 
        DT2sp[i] = (area[1] * D_DiffCoeff[1,ny] + area[2] * D_DiffCoeff[i,ny-1]) / Tarea 
    end
    MDnp = sum(DT2np)
    MDsp = sum(DT2sp) 

    GCnp = geom * MDnp + B_coeff + 2 * C_HeatCapacity[1,1] * NT     # north pole 
    GCsp = geom * MDsp + B_coeff + 2 * C_HeatCapacity[1,ny] * NT    # south pole

    for i=1:nx
        #North pole
        row_idx = i
        A[row_idx,row_idx] = GCnp
        for ii=1:nx
            A[row_idx,nx+ii] = - geom * DT2np[ii]
        end
        
        #South pole
        row_idx = dof + 1 - i
        A[row_idx,row_idx] = GCsp
        for ii=1:nx
            A[row_idx,dof-2*nx+ii] = - geom * DT2sp[ii]
        end
    end
    return A
end

function UpdateRHS!(RHS, mesh, NT, time_step, Temp, model, LastRHS)
    @unpack nx,ny = mesh
    @unpack C_HeatCapacity, SolarForcing, A_coeff = model
    for j=1:ny
        for i=1:nx
            row_idx = index1d(i,j,nx)

            RHS[row_idx] = 4 * C_HeatCapacity[i,j] * Temp[row_idx] * NT  + LastRHS[row_idx] - 2 * A_coeff
            if (time_step == 1)
                RHS[row_idx] += SolarForcing[i,j,NT] + SolarForcing[i,j,time_step]
            else
                RHS[row_idx] += SolarForcing[i,j,time_step-1] + SolarForcing[i,j,time_step]
            end
        end
    end
    LastRHS .= RHS
end
