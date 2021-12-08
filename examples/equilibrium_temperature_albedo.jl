using Klimakoffer

NT = 48 

mesh = Mesh()

co2_concentration = 315.0

model = Model(mesh, NT; co2_concentration = co2_concentration,compute_albedo=true)

discretization1 = Discretization(mesh, model, NT)

GlobTemp = compute_equilibrium!(discretization1)
