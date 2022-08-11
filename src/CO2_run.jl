ENV["GKSwstype"] = "100"
import Pkg
Pkg.instantiate()

using Interpolations
using UnPack
using Plots
using LinearAlgebra
using SparseArrays


include("CO2.jl")
include("visualization.jl")
include("model.jl")
include("mesh.jl")
include("numerics.jl")
include("discretization.jl")



#Constants
mesh = Mesh()
@unpack nx,ny = mesh

NT = 48
year_start = 1958
year_end = 2020
year_delta = year_end - year_start + 1


#=
First create the location-dependent CO2 matrices, so that they match the global average and have the appropriate size.
Fill the missing values from line 60 with the respective values of the corresponding month from the South Pole and exchange all wrong values above
with the average global CO2 concentration.
=#

read_south_pole_CO2_values()

build_CO2_matrix_with_replaced_values(num,yearSP,monthSP)
CO2_bilinear_interpolation(nasa_co2_mat_mod,nlongitude,counter)
#bilinear_one_matrix(nasa_co2_mat_mod[1],nlongitude)
bisection(bilinear_co2_mat,nasa_co2_mat_mod,mesh.area,CO2_averages,CO2_averages_SP,counter)
extend_co2_matrices(mesh,final_co2_mat,year_start,year_end,counter)


#=
Case 1: CO2 matrices are time and location dependent. The data found range from the year 2003 to the year 2017. 
All other years were created using matrices that were given the respective global average value of the corresponding month for each node.
Thus, the matrices outside the range are only time-dependent and no longer location-dependent. This serves to compare all cases.
=#



println("Case 1")
mean_co2_1 = calc_mean_co2(more_years_co2_mat,mesh.area)
model = Model(mesh,NT,more_years_co2_mat[1])
discretization = Discretization(mesh,model,NT)
solsol = compute_equilibrium!(discretization)
sol = compute_evolution!(discretization,more_years_co2_mat,year_start,year_end)


#=
Case 2: CO2 matrices are time-dependent and not location-dependent. The data comes from the co2_mm_mlo file.
=#
println("Case 2")
co2_values = CO2_Array[1:end,4]
model2 = Model(mesh,NT,co2_values)
discretization2 = Discretization(mesh,model2,NT)
compute_equilibrium!(discretization2)
sol2 = compute_evolution!(discretization2,co2_values,year_start,year_end)
   

#=
Case 3: CO2 matrices are neither time nor location dependent. 
=#
println("Case 3")
ind_year_start = findfirst(CO2_Array[:,1] .== 2003)
same_co2 = zeros(Float64,year_delta*12+1)
same_co2_vec = []

for i in 1:year_delta*12+1
  same_co2[i] = co2_values[ind_year_start+2]
  push!(same_co2_vec,same_co2[i])
end
#mean_co2_2 = calc_mean_co2(same_co2_vec,mesh.area)
model3 = Model(mesh,NT,same_co2_vec[1])
discretization3 = Discretization(mesh,model3,NT)
init_temp1=compute_equilibrium!(discretization3)
sol3=compute_evolution!(discretization3,same_co2_vec,2003,2016)



#=
Case 4: CO2 matrices are not time-dependent, but location-dependent. Here we start the simulation in year 2003 because of missing data.
=#

println("Case 4")
co2_matrix = []
#co2_matrix = zeros(Float64,(nx,ny))
#for i in 1:nx
#  for j in 1:ny
#      co2_matrix[i,j] = final_co2_mat[1][i,j]
#  end
#end
for i in 1:169
  push!(co2_matrix,final_co2_mat[1])
end
mean_co2_3 = calc_mean_co2(co2_matrix,mesh.area)
model4 = Model(mesh,NT,co2_matrix[1])
discretization4 = Discretization(mesh,model4,NT)
init_temp2=compute_equilibrium!(discretization4)
sol4=compute_evolution!(discretization4,co2_matrix,2003,2016)


"""
From here on, the temperature is calculated from the past to the future. OLS and exponential regression were implemented.
"""
#constants
future_year_start_loc_diff = 2017
future_year_start = 2022
future_year_end = 2100

#=
Case 5: CO2 matrices are time and location dependent. In this case, exponential regression is used for the calculation.
=#

println("Case 5")
calc_ols(mesh,final_co2_mat,2003,2016,NT;linear=false)
calc_co2_matrices(OLS_coefficients,future_year_start_loc_diff,future_year_end,NT;linear=false)
combine_matrices(final_co2_mat,val,counter)


#= 
To achieve a more exact result, now the vector "var_co2_mat" is going to be extended with matrices of the past
=#
(row,) = size(var_co2_mat)
co2_mat_extended = []
x = findfirst(CO2_Array[:,1] .==2003)
y = findfirst(CO2_Array[:,1] .==year_start)
for i in y:x+1
  z = fill(CO2_Array[i,4],nx,ny)
  push!(co2_mat_extended,z)
end
for i in 1:row
  push!(co2_mat_extended,var_co2_mat[i])
end
mean_co2_4 = calc_mean_co2(co2_mat_extended,mesh.area)
model5 = Model(mesh,NT,co2_mat_extended[1])
discretization5 = Discretization(mesh,model5,NT)
compute_equilibrium!(discretization5)
sol5 = compute_evolution!(discretization5,co2_mat_extended,year_start,future_year_end-1)


#=
Case 6: CO2 matrices are time and location dependent. In this case, linear regression is used for the calculation.
=#

println("Case 6")
calc_ols(mesh,final_co2_mat,2003,2016,NT;linear=true)
calc_co2_matrices(OLS_coefficients,future_year_start_loc_diff,future_year_end,NT;linear=true)
combine_matrices(final_co2_mat,val,counter)

co2_mat_extended_1 = []
x = findfirst(CO2_Array[:,1] .==2003)
y = findfirst(CO2_Array[:,1] .==year_start)
for i in y:x+1
  z = fill(CO2_Array[i,4],nx,ny)
  push!(co2_mat_extended_1,z)
end
for i in 1:row
  push!(co2_mat_extended_1,var_co2_mat[i])
end
mean_co2_5 = calc_mean_co2(co2_mat_extended_1,mesh.area)
model6 = Model(mesh,NT,co2_mat_extended_1[1])
discretization6 = Discretization(mesh,model6,NT)
compute_equilibrium!(discretization6)
sol6 = compute_evolution!(discretization6,co2_mat_extended_1,year_start,future_year_end-1)


#=
Case 7: CO2 matrices are time-dependent and not location-dependent. In this case, exponential regression is used for the calculation.
=#

println("Case 7")
calc_ols(mesh,co2_values,year_start,year_end+1,NT;linear=false)
calc_co2_matrices(beta,future_year_start,future_year_end,NT;linear=false)
mat1 = combine_matrices(co2_values,val,counter)


model7 = Model(mesh,NT,mat1)
discretization7 = Discretization(mesh,model7,NT)
compute_equilibrium!(discretization7)
sol7 = compute_evolution!(discretization7,mat1,year_start,future_year_end-1)


#=
Case 8: CO2 matrices are time-dependent and not location-dependent. In this case, linear regression is used for the calculation.
=#

println("Case 8")
calc_ols(mesh,co2_values,year_start,year_end+1,NT;linear=true)
calc_co2_matrices(beta,future_year_start,future_year_end,NT;linear=true)
mat2 =combine_matrices(co2_values,val,counter)

model8 = Model(mesh,NT,mat2)
discretization8 = Discretization(mesh,model8,NT)
compute_equilibrium!(discretization8)
sol8 = compute_evolution!(discretization8,mat2,year_start,future_year_end-1)
  



 

 using Downloads

temp_array=readdlm(Downloads.download("https://data.giss.nasa.gov/gistemp/graphs/graph_data/Global_Mean_Estimates_based_on_Land_and_Ocean_Data/graph.txt"); skipstart=5)
ind_1958 = findfirst(temp_array[:,1] .== 1958)
ind_2021 = findfirst(temp_array[:,1] .== 2021)

time_scale_short = range(year_start+2/12,year_end+1+2/12,(year_end+1-year_start)*12+1)
time_scale_tiny = range(2003+2/12,2017+2/12,(2017-2003)*12+1)
time_scale_long = range(year_start+2/12,future_year_end+2/12,(future_year_end-year_start)*12+1)


plot(sol2.year_array,sol2.mean_temperature_yearly,label="KK t var x const")
plot!(sol.year_array,sol.mean_temperature_yearly,label="KK t, x var")
ylabel!("Mean annual temperature [°C]")
xlabel!("year") 
savefig("plot_time_variable.png")

plot(time_scale_short,mean_co2_1,label="KK t, x var")
plot!(time_scale_short,co2_values[1:757],label="KK t var x const")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_plot_t_variable.png")


plot(sol3.year_array,sol3.mean_temperature_yearly,label="KK t, x const")
plot!(sol4.year_array,sol4.mean_temperature_yearly,label="KK t const, x var")
ylabel!("Mean annual temperature [°C]")
xlabel!("year") 
savefig("plot_time_const.png")

plot(time_scale_tiny,same_co2_vec[1:169],label="KK t, x const")
plot!(time_scale_tiny,mean_co2_3,label="KK t const x var")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_plot_t_const.png")


plot(sol5.year_array,sol5.mean_temperature_yearly,label="KK t, x var exp")
plot!(sol6.year_array,sol6.mean_temperature_yearly,label="KK t, x var lin")
ylabel!("Mean annual temperature [°C]")
xlabel!("year") 
savefig("ols_plot_t_x_var.png")

plot(time_scale_long,mean_co2_4,label="KK t, x var exp")
plot!(time_scale_long,mean_co2_5,label="KK t, x var lin")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_ols_plot_t_x_var.png")


plot(sol7.year_array,sol7.mean_temperature_yearly,label="KK t var x const exp")
plot!(sol8.year_array,sol8.mean_temperature_yearly,label="KK t var x const lin")
ylabel!("Mean annual temperature [°C]")
xlabel!("year") 
savefig("ols_plot_t_var_x_const.png")

plot(time_scale_long,mat1,label="KK t var x const exp")
plot!(time_scale_long,mat2,label="KK t var x const lin")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_ols_plot_t_var_x_const.png")


plot(sol5.year_array,sol5.mean_temperature_yearly,label="KK t, x var exp")
plot!(sol7.year_array,sol7.mean_temperature_yearly,label="KK t var x const exp")
ylabel!("Mean annual temperature [°C]")
xlabel!("year") 
savefig("ols_exp_plot_t_var_x_const.png")

plot(time_scale_long,mean_co2_4,label="KK t, x var exp")
plot!(time_scale_long,mat1,label="KK t var x const exp")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_ols_exp_plot_t_var_x_const.png")


plot(sol6.year_array,sol6.mean_temperature_yearly,label="KK t, x var lin")
plot!(sol8.year_array,sol8.mean_temperature_yearly,label="KK t var x const lin")
ylabel!("Mean annual temperature [°C]")
xlabel!("year") 
savefig("ols_lin_plot_t_var_x_const.png")

plot(time_scale_long,mean_co2_5,label="KK t, x var lin")
plot!(time_scale_long,mat2,label="KK t var x const lin")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_ols_lin_plot_t_var_x_const.png")




#=
#plotting co2
time_scale_short = range(year_start+2/12,year_end+1+2/12,(year_end+1-year_start)*12+1)

plot(time_scale_short,mean_co2_1,label="KK t, x var")
plot!(time_scale_short,co2_values[1:757],label="KK t var x const")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_plot_t_variable.png")

time_scale_tiny = range(2003+2/12,2017+2/12,(2017-2003)*12+1)
plot(time_scale_tiny,same_co2_vec[1:169],label="KK t, x const")
plot!(time_scale_tiny,mean_co2_3,label="KK t const x var")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_plot_t_const.png")

time_scale_long = range(year_start+2/12,future_year_end+2/12,(future_year_end-year_start)*12+1) # adding 2/12 because of start in 
#time_scale_KK = time_scale[3:end]
plot(time_scale_long,mean_co2_4,label="KK t, x var lin")
plot!(time_scale_long,mat1,label="KK t var x const lin")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_t_var_x_const.png")

plot(time_scale_long,mean_co2_5,label="KK t, x var lin")
plot!(time_scale_long,mat2,label="KK t var x const lin")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_ols_lin_plot_t_var_x_const.png")

plot(time_scale_long,mean_co2_4,label="KK t, x var exp")
plot!(time_scale_long,mean_co2_5,label="KK t, x var exp")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_ols_plot_t_x_var.png")

plot(time_scale_long,mat1,label="KK t, x var lin")
plot!(time_scale_long,mat2,label="KK t var x const lin")
ylabel!("Mean annual co2 concentration [ppm]")
xlabel!("year") 
savefig("co2_ols_lin_plot_t_var_x_const.png")
=#


#plot!(sol1.year_array,temp_array[ind1:ind2-1,2].-temp_array[ind1,2].+sol1.mean_temperature_yearly[1],label="NASA") 
#plot!(sol1.year_array,temp_array[ind1:ind2-1,3].-temp_array[ind1,2].+sol1.mean_temperature_yearly[1],label="NASA smoothed")
 #plot!(sol1.year_array,temp_array[ind1:ind2-1,2].-temp_array[ind1,2].+14.0,label="NASA") 
 #plot!(sol1.year_array,temp_array[ind1:ind2-1,3].-temp_array[ind1,2].+14.0,label="NASA smoothed")





# Plotting co2 concentration 
# month_steps1 = range(1959,2020,744)
# month_steps1_1 = range(1959,2021,756)
#x = 1959:2020:744
# println(x)
#ind_1959 = findfirst(CO2_Array[:,4].==1959)
#plot!(month_steps1,co2_values[11:754],label="CO2_average-Verlauf")
#ylabel!("CO2 (ppm)")
#xlabel!("year")
#savefig("CO2-Verlauf")

