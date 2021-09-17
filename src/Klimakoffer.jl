module Klimakoffer

using LinearAlgebra: \, lu
using SparseArrays: sparse
using UnPack: @unpack
using Requires: @require

export Mesh, Model, Discretization, compute_equilibrium!, compute_evolution!, set_co2_concentration!

function __init__()
    # Enable features that depend on the availability of the Plots package
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
      using .Plots: Plots
      include("visualization.jl")
    end

    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
      using .Makie: Makie
    end
end

# Include additional files
include("mesh.jl")
include("model.jl")
include("discretization.jl")
include("numerics.jl")

end # module
