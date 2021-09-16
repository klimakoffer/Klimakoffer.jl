module Klimakoffer

using LinearAlgebra: \, lu
using SparseArrays: sparse
using UnPack: @unpack

export main

# Include additional files
include("mesh.jl")
include("model.jl")
include("numerics.jl")

end # module
