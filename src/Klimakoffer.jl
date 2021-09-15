module Klimakoffer

include("Klimaparameter.jl")

export calc_CO2_concentration_A, calc_heat_capacities_C, calc_diffusion_coefficients, calc_diffusion_coefficients_poles, orbital_params, read_albedo, read_world

export answer, Mesh

"""
    answer()

Return the Answer to the Ultimate Question of Life, The Universe, and Everything (see
Douglas Adams' "The Hitchhiker's Guide to the Galaxy" for more details).
"""
answer() = 42

# Include additional files
include("mesh.jl")

end # module
