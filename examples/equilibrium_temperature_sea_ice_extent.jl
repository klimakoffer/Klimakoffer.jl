using Klimakoffer

NT = 48 

mesh = Mesh()

co2_concentration = 315.0

model = Model(mesh, NT; co2_concentration = co2_concentration,compute_albedo=true, compute_sea_ice_extent = true, sea_ice_extent_year = 1979)

discretization1 = Discretization(mesh, model, NT)

GlobTemp = compute_equilibrium!(discretization1, update_heat_capacity = true ,update_solar_forcing = true)

