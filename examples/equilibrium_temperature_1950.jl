using Klimakoffer

NT = 48 # Number of time-steps per year

mesh = Mesh()

co2_concentration = 315.0 # in [ppm] from year 1950
model = Model(mesh, NT; co2_concentration = co2_concentration)

discretization = Discretization(mesh, model, NT)

GlobTemp = compute_equilibrium!(discretization)
