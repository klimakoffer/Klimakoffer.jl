using Klimakoffer
using Printf

NT = 48 # Number of time-steps per year

mesh = Mesh()

model = Model(mesh, NT)

discretization = Discretization(mesh, model, NT)

year_start = 1950
year_end = 2020
year_delta = year_end - year_start + 1

co2_concentration_yearly = zeros(Float64,year_delta)
mean_temperature_yearly = zeros(Float64, year_delta)

co2_concentration_start = 315.01
co2_concentration_yearly_delta = 1.39


for years in 1950:2020
  iter = years+1-year_start
  co2_concentration_yearly[iter] = co2_concentration_start + co2_concentration_yearly_delta * (iter - 1)
  mean_temperature_yearly[iter] = compute_equilibrium!(discretization; co2_concentration=co2_concentration_yearly[iter], verbose = false)
  #@printf "Global Temerature in year ", years, " with CO2=", co2_concentration_yearly[iter], " [ppm] is T=", mean_temperature_yearly[iter], " °C"
  @printf "Global Temerature in year %4i with CO2 concentration %.2f [ppm] is T=%.3f °C" years co2_concentration_yearly[iter] mean_temperature_yearly[iter]
  println("")
end
