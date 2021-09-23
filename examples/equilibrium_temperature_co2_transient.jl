# In this test, we simulate the global temperature from 1958 to 2020 using the 
# historical data of CO2 concentration in the atmosphere
# Source: https://climate.nasa.gov/vital-signs/carbon-dioxide/

using Klimakoffer
using DelimitedFiles

NT = 48 # Number of time-steps per year (must be a multiple of 12)

year_start = 1958
year_end = 2020
year_delta = year_end - year_start + 1

# Read co2 concentration from NASA's data base
co2_array=readdlm(download("ftp://aftp.cmdl.noaa.gov/products/trends/co2/co2_mm_mlo.txt"); skipstart=53)

co2_concentration_at_step = zeros(Float64,year_delta*NT+1)
step=0
total_month = 0
for years = year_start:year_end
  for month=1:12
    global total_month += 1
    for s=1:NT/12
      global step += 1
      co2_concentration_at_step[step] = co2_array[total_month,4]
    end
  end
end
co2_concentration_at_step[end] = co2_array[total_month + 1,4]

mesh = Mesh()
model = Model(mesh, NT; co2_concentration = co2_concentration_at_step[1])
discretization = Discretization(mesh, model, NT)

println("First compute initial temperature distribution when equilibrium is reached")
initial_temp = compute_equilibrium!(discretization; verbose = false)   
println("Start evolution of temperature...")
sol = compute_evolution!(discretization, co2_concentration_at_step, year_start, year_end; verbose=true) 
println("Is it getting hotter?")
