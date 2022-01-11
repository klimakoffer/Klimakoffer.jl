
struct Discretization{DecType}
  lu_dec::DecType                       # LU decomposition of system matrix
  num_steps_year::Int64                 # Number of time steps per astronomical year
  mesh::Mesh
  model::Model                          # Physical model
  annual_temperature::Array{Float64, 2} # Annual temperature at every time step in an astronomical year
  rhs::Vector{Float64}
  last_rhs::Vector{Float64}
end


function Discretization(mesh, model, num_steps_year; run_garbage_collector = true)
    lu_dec = compute_lu_matrices(mesh, model, num_steps_year)

    if run_garbage_collector
      GC.gc()
    end

    annual_temperature = fill(5.0, mesh.dof, num_steps_year) # Magic initialization
    rhs     = zeros(mesh.dof)  # TODO: The EBM Fortran code initializes the RHS to zero... Maybe we want to initialize it differently
    last_rhs = zeros(mesh.dof)

    Discretization(lu_dec, num_steps_year, mesh, model, annual_temperature, rhs, last_rhs)
end

Base.size(discretization::Discretization) = size(discretization.mesh)

function Base.show(io::IO, discretization::Discretization)
  nx, ny = size(discretization)
  print(io, "Discretization() with ", nx, "×", ny, " degrees of freedom")
end

function compute_lu_matrices(mesh, model, num_steps_year)
  mat = compute_matrix(mesh,num_steps_year,model)
  lu_dec = lu(mat)
  return lu_dec
end