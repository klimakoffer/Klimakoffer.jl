using Klimakoffer

NT = 48 # Number of time-steps per year

mesh = Mesh()

model = Model(mesh, NT)

discretization = Discretization(mesh, model, NT)

GlobTemp = compute_equilibrium!(discretization)
