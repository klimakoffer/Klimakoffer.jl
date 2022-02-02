#Mapping functions
##################
# function index2d(k,nx)
#     j = floor(Int,(k-1)/nx)+1
#     # TODO: add consistency check with ny?
#     return k-(j-1)*nx,j 
# end

function index1d(i,j,nx)
    # TODO: add consistency check with ny?
    return i+(j-1)*nx
end


"""
    compute_equilibrium!(...)
* rel_error is the tolerance for global temperature equilibrium (default is 2e-5).
* max_years is the maximum number of annual cycles to be computed when searching for equilibrium
* update_heat_capacity : manage the monthly update function for the heat capacity
* update_solar_forcing : manage the monthly update function for the solar forcing
"""
function compute_equilibrium!(discretization; max_years=100, rel_error=2e-5, verbose=true, update_heat_capacity = false,update_solar_forcing = false)
    @unpack mesh, model, num_steps_year, annual_temperature, rhs, last_rhs  = discretization
    @unpack nx, dof = mesh
    
    average_temperature = average_temperature_old = area_weighted_average(view(annual_temperature, :, num_steps_year), mesh)

    if verbose
      println("year","  ","Average Temperature")
      println(0,"  ",average_temperature_old)
    end

    for year in 1:max_years
        average_temperature = 0.0
        for time_step in 1:num_steps_year

            update_monthly_params!(model, discretization.lu_decomposition, mesh, num_steps_year, time_step, update_heat_capacity, update_solar_forcing)

            old_time_step = (time_step == 1) ? num_steps_year : time_step - 1
            update_rhs!(rhs, mesh, num_steps_year, time_step, view(annual_temperature, :, old_time_step), model, last_rhs)
                        
            # Use in-place operation `ldiv!` instead of `\` to avoid allocations
            ldiv!(view(annual_temperature, :, time_step), discretization.lu_decomposition, rhs)

            annual_temperature[1:nx, time_step] .= annual_temperature[1, time_step]
            annual_temperature[dof-nx+1:dof, time_step] .= annual_temperature[dof, time_step]

            average_temperature += area_weighted_average(view(annual_temperature, :, time_step), mesh)
        end
        average_temperature = average_temperature/num_steps_year
        if verbose
          println(year,"  ",average_temperature)
        end
        
        if (abs(average_temperature-average_temperature_old)<rel_error)
            if verbose
              println("EQUILIBRIUM REACHED!")
            end
            break
        end
        
        average_temperature_old = average_temperature
    end

    return average_temperature
end

"""
    compute_evolution!(...)
Compute the evolution of the mean temperature with varying CO2 levels.
"""
function compute_evolution!(discretization, co2_concentration_at_step, year_start, year_end; verbose=true, update_heat_capacity = false,update_solar_forcing = false)
    @unpack mesh, model, num_steps_year, annual_temperature, rhs, last_rhs  = discretization
    @unpack nx, dof = mesh
    
    if verbose
      println("year","  ","Average Temperature")
    end
    
    max_years = year_end - year_start + 1

    # Allocate output arrays
    year_at_step = zeros(Float64,max_years*num_steps_year+1)
    mean_temperature_at_step = zeros(Float64,max_years*num_steps_year+1)
    year_array = zeros(Float64,max_years)
    mean_temperature_yearly = zeros(Float64,max_years)

    # Fill year arrays
    year_at_step[1] = year_start + 0.21644  # The simulation starts at the vernal equinox
    for year in 2:size(year_at_step,1)
        year_at_step[year] = year_at_step[year-1] + 1.0/num_steps_year
    end
    year_array[1] = year_start + 0.71644    # The simulation starts at the vernal equinox
    for year in 2:max_years
        year_array[year] = year_array[year-1] + 1.0
    end

    # Initialize the temperature array
    mean_temperature_at_step[1] = area_weighted_average(view(annual_temperature, :, num_steps_year), mesh)
    step = 1

    for year in 1:max_years
        average_temperature = 0.0
        for time_step in 1:num_steps_year
            set_co2_concentration!(model, co2_concentration_at_step[step])
            update_monthly_params!(model, discretization.lu_decomposition, mesh, num_steps_year, time_step, update_heat_capacity, update_solar_forcing)   
            old_time_step = (time_step == 1) ? num_steps_year : time_step - 1
            update_rhs!(rhs, mesh, num_steps_year, time_step, view(annual_temperature, :, old_time_step), model, last_rhs)
                        
            # Use in-place operation `ldiv!` instead of `\` to avoid allocations
            ldiv!(view(annual_temperature, :, time_step), discretization.lu_decomposition, rhs)

            annual_temperature[1:nx, time_step] .= annual_temperature[1, time_step]
            annual_temperature[dof-nx+1:dof, time_step] .= annual_temperature[dof, time_step]

            # Compute mean temperature
            step += 1
            mean_temperature_at_step[step] = area_weighted_average(view(annual_temperature, :, time_step), mesh)
            average_temperature += mean_temperature_at_step[step]
        end
        average_temperature = average_temperature/num_steps_year
        mean_temperature_yearly[year] = average_temperature
        if verbose
          println(year_start+year-1,"  ",average_temperature)
        end
    end

    return (; year_at_step, mean_temperature_at_step, year_array, mean_temperature_yearly)
end


"""
Computes mean temperature in the globe at a specific time
"""
function area_weighted_average(vector,mesh)
    @unpack nx,ny,area,dof = mesh

    average_vector = 0.0
    
    # Contribution of inner points
    for j=2:ny-1
        for i=1:nx
            row_idx = index1d(i,j,nx)
            average_vector += area[j] * vector[row_idx]
        end
    end

    # Poles
    average_vector += area[1] * (vector[1] + vector[dof]) 

    return average_vector
end

"""
compute_matrix(mesh,num_steps_year,model)
Compute the system matrix for the linear solver
"""
function compute_matrix(mesh,num_steps_year,model)
    @unpack nx,ny,dof,h,geom,csc2,cot,area = mesh
    @unpack diffusion_coeff,heat_capacity,radiative_cooling_feedback = model
    matrix = spzeros(Float64,dof,dof)
    sh2 = 1 / h^2
    # Inner DOFs (c coefficients are divided by h^2)
    for j=2:ny-1
        for i=1:nx

            # Compute coefficients
            c0 = 2 * num_steps_year + (2 * sh2 * diffusion_coeff[i,j] * (1 + csc2[j]) + radiative_cooling_feedback) / heat_capacity[i,j,2]

            # Get diffusion coefficient 
            if (i == 1) # Periodic BC
                d_phi = (diffusion_coeff[2,j] - diffusion_coeff[nx,j]) / (2 * h)
            elseif (i == nx) # Periodic BC
                d_phi = (diffusion_coeff[1,j] - diffusion_coeff[nx-1,j]) / (2 * h)
            else # Inner DOFs
                d_phi = (diffusion_coeff[i+1,j] - diffusion_coeff[i-1,j]) / (2 * h)
            end

            c1 = sh2 * csc2[j] * (diffusion_coeff[i,j] - 0.5 * h * d_phi)
            c3 = sh2 * csc2[j] * (diffusion_coeff[i,j] + 0.5 * h * d_phi)

            d_theta= (diffusion_coeff[i,j+1] - diffusion_coeff[i,j-1]) / (2 * h)
            c2 = sh2 * (diffusion_coeff[i,j] - 0.5 * h * (diffusion_coeff[i,j] * cot[j] + d_theta))
            c4 = sh2 * (diffusion_coeff[i,j] + 0.5 * h * (diffusion_coeff[i,j] * cot[j] + d_theta))

            # Fill matrix A
            row_idx = index1d(i,j,nx)
            col_idx_c0 = row_idx
            
            if (i==1) # Periodic BC
                col_idx_c1 = index1d(nx,j,nx)
                col_idx_c3 = row_idx + 1
            elseif (i==nx) # Periodic BC
                col_idx_c1 = row_idx - 1
                col_idx_c3 = index1d(1,j,nx)
            else
                col_idx_c1 = row_idx - 1
                col_idx_c3 = row_idx + 1
            end
            col_idx_c2 = row_idx - nx
            col_idx_c4 = row_idx + nx

            matrix[row_idx,col_idx_c0] = c0
            matrix[row_idx,col_idx_c1] = -c1 / heat_capacity[i,j,2]
            matrix[row_idx,col_idx_c2] = -c2 / heat_capacity[i,j,2]
            matrix[row_idx,col_idx_c3] = -c3 / heat_capacity[i,j,2]
            matrix[row_idx,col_idx_c4] = -c4 / heat_capacity[i,j,2]

            # Get rhs
            # rhs[row_idx] = -(c0 * Temp[col_idx_c0] - c1 * Temp[col_idx_c0] - c2 * Temp[col_idx_c2] - c3 * Temp[col_idx_c3] - c4 * Temp[col_idx_c4]) + 4 * heat_capacity[i,j]
        end
    end

    # Poles
    total_area = area[1] + area[2]

    diff_total_np = zeros(Float64,nx)
    diff_total_sp = zeros(Float64,nx)
    for i=1:nx
        diff_total_np[i] = (area[1] * diffusion_coeff[1,1]  + area[2] * diffusion_coeff[i,2]) / total_area 
        diff_total_sp[i] = (area[1] * diffusion_coeff[1,ny] + area[2] * diffusion_coeff[i,ny-1]) / total_area 
    end
    diff_total_sum_np = sum(diff_total_np)
    diff_total_sum_sp = sum(diff_total_sp) 

    gc_np = 2 * num_steps_year + (geom * diff_total_sum_np + radiative_cooling_feedback) / heat_capacity[1,1,2]    # north pole 
    gc_sp = 2 * num_steps_year + (geom * diff_total_sum_sp + radiative_cooling_feedback) / heat_capacity[1,ny,2]    # south pole 

    for i=1:nx
        #North pole
        row_idx = i
        matrix[row_idx,row_idx] = gc_np
        for ii=1:nx
            matrix[row_idx,nx+ii] = - geom * diff_total_np[ii] / heat_capacity[1,1,2]
        end
        
        #South pole
        row_idx = dof + 1 - i
        matrix[row_idx,row_idx] = gc_sp
        for ii=1:nx
            matrix[row_idx,dof-2*nx+ii] = - geom * diff_total_sp[ii] / heat_capacity[1,ny,2]
        end
    end
    return matrix
end

function update_rhs!(rhs, mesh, num_steps_year, time_step, temperature, model,  last_rhs)
    @unpack nx,ny = mesh
    @unpack heat_capacity, solar_forcing, radiative_cooling_co2 = model
    for j=1:ny
        for i=1:nx
            row_idx = index1d(i,j,nx)

            rhs[row_idx] = 4 * temperature[row_idx] * num_steps_year  - last_rhs[row_idx] + solar_forcing[i,j,time_step] / heat_capacity[i,j,2] - (radiative_cooling_co2[2] / heat_capacity[i,j,2] + radiative_cooling_co2[1] / heat_capacity[i,j,1])
            if (time_step == 1)
                rhs[row_idx] += solar_forcing[i,j,num_steps_year] / heat_capacity[i,j,1]
            else
                rhs[row_idx] += solar_forcing[i,j,time_step-1] / heat_capacity[i,j,1]
            end
        end
    end
    last_rhs .= rhs
end



"""
 update_monthly_params!(...)
* The routine updates monthly all important parameters of the model,
  which are necessary for the passed update settings.
  When updating the heat capacity, the LU decomposition must take place again. 
  Update start is the first week of April.
"""
function update_monthly_params!(model, lu_decomposition, mesh, num_steps_year, time_step, update_heat_capacity, update_solar_forcing)
   
    if mod(time_step,4) == 2 && time_step in (2:48)
        if update_heat_capacity == true && update_solar_forcing == true     
         update_model(model, mesh,time_step, true,true,false,true,true)
         lu_decomposition = compute_lu_matrices(mesh, model, num_steps_year)
        elseif update_heat_capacity == true && update_solar_forcing == false 
            update_model(model, mesh, time_step, true,false,false,true,false)
            lu_decomposition = compute_lu_matrices(mesh, model, num_steps_year)   
        elseif update_heat_capacity == false && update_solar_forcing == true
            update_model(model, mesh, time_step, true,true,false,false,true)  
        elseif update_heat_capacity == false && update_solar_forcing == false 
            update_model(model, mesh, time_step, false,false,false,false,false)
        end
    else
        if update_heat_capacity
            model.heat_capacity[:,:,1] = model.heat_capacity[:,:,2]
        end
    end
end