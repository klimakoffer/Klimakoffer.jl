using Klimakoffer

NT = 48 # Number of time-steps per year

year_start = 1950
year_end = 1960
year_delta = year_end - year_start + 1

co2_concentration_start = 315.0
co2_concentration_yearly_delta = 1.39

co2_concentration_yearly = zeros(Float64,year_delta)
mean_temperature_yearly = zeros(Float64, year_delta)

mesh = Mesh()
model = Model(mesh, NT; co2_concentration = co2_concentration_start)
discretization = Discretization(mesh, model, NT)

for years = year_start:year_end
  iter = years + 1 - year_start
  co2_concentration_yearly[iter] = co2_concentration_start + co2_concentration_yearly_delta * (iter - 1)
end
println("First compute initial temperature distribution when equilibrium is reached")
mean_temperature_yearly[1] = compute_equilibrium!(discretization; verbose = false)   
println("Start evolution of temperature...")
compute_evolution!(discretization, co2_concentration_yearly, mean_temperature_yearly, year_start; verbose=true) 
println("Is it getting hotter?")
