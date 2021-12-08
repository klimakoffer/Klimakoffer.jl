using Klimakoffer
using Printf

NT = 48 # Number of time-steps per year

year_start = 1950
year_end = 1951
year_delta = year_end - year_start + 1

co2_concentration_start = 315.0 # in [ppm] for year 1950
co2_concentration_yearly_delta = 1.39

co2_concentration_yearly = zeros(Float64,year_delta)
mean_temperature_yearly = zeros(Float64, year_delta)
# co2 scenario: simple linear growth per year
for years in year_start:year_end
  iter = years+1-year_start
  co2_concentration_yearly[iter] = co2_concentration_start + co2_concentration_yearly_delta * (iter - 1)
end

mesh = Mesh()
model = Model(mesh, NT; co2_concentration = co2_concentration_start)
discretization = Discretization(mesh, model, NT)

for years in year_start:year_end
  iter = years+1-year_start
  set_co2_concentration!(model, co2_concentration_yearly[iter])
  mean_temperature_yearly[iter] = compute_equilibrium!(discretization; verbose = false)
  @printf "Global Temerature in year %4i with CO2 concentration %.2f [ppm] is T=%.3f °C" years co2_concentration_yearly[iter] mean_temperature_yearly[iter]
  println("")
end
GlobTemp = mean_temperature_yearly[year_end - year_start + 1]
