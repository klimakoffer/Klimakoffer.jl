module Klimakoffer

using LinearAlgebra: \, lu
using SparseArrays: sparse
using UnPack: @unpack

export Mesh, Model, Discretization, compute_equilibrium!

# Include additional files
include("mesh.jl")
include("model.jl")
include("discretization.jl")
include("numerics.jl")

end # module
