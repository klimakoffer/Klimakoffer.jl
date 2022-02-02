using Klimakoffer

NT = 48 

mesh = Mesh(150)

co2_concentration = 315.0

model = Model(mesh, NT; co2_concentration = co2_concentration, compute_albedo=true)

disc = Discretization(mesh, model, NT)

GlobTemp = compute_equilibrium!(disc)