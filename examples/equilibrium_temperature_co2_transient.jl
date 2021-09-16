using Klimakoffer

NT = 48 # Number of time-steps per year

mesh = Mesh()

model = Model(mesh, NT)

discretization = Discretization(mesh, model, NT)

year_start = 1950
year_end = 2020
year_delta = year_end-year_start+1

co2_concentration_yearly = zeros(Float64,year_delta)
mean_temperature_yearly = zeros(Float64, year_delta)

co2_concentration_start = 315.01
co2_concentration_yearly_delta = 1.39

for years = 1950:2020
  iter = years+1-year_start
  co2_concentration_yearly[iter] = co2_concentration_start + co2_concentration_yearly_delta * (iter - 1)
end
mean_temperature_yearly[1] = compute_equilibrium!(discretization; co2_concentration=co2_concentration_yearly[1], verbose = false)   
compute_evolution!(discretization, co2_concentration_yearly, mean_temperature_yearly, year_start; verbose=true) 
