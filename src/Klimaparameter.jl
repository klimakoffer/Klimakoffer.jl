  function calc_CO2_concentration_A(CO2ppm)

    # Define base values CO2_Base and A_Base
    CO2_Base = 315.0
    A_Base = 210.3

    # Doesn't change as long as CP2ppm doesn't change
    A=A_Base-5.35*log(CO2ppm/CO2_Base)

    return A
  end


  function calc_heat_capacities_C(geography, B)

    # Depths 
    depth_atmos = 5000.         # meters
    depth_mixed_layer = 70.0    # meters
    depth_soil = 2.0            # meters
    depth_seaice = 2.5          # meters
    depth_snow = 2.0            # meters
    layer_depth = 0.5           # kilometers

    # Physical properties of atmosphere
    rho_atmos = 1.293           # kg m^-3  dry air (STP)
    csp_atmos = 1005.0          # J kg^-1 K^-1 (STP)
    scale_height = 7.6          # kilometers
        
    # Physical properties of water
    rho_water = 1000.0          # kg m^-3
    csp_water = 4186.0          # J kg^-1 K^-1

    # Physical properties of soil 
    rho_soil = 1100.0           # kg m^-3   
    csp_soil = 850.0            # J kg^-1 K^-1
      
    # Physical properties of sea ice
    rho_sea_ice = 917.0         # kg m^-3
    csp_sea_ice = 2106.0        # J kg^-1 K^-1  

    # Physical properties of snow covered surface
    rho_snow = 400.0            # kg m^-3
    csp_snow = 1900.0           # J kg^-1 K^-1
       
    # Other constants  
    sec_per_yr = 3.15576e7      # seconds per year
    days_per_yr = 365.2422      # days per year


    # atmosphere with exponentially decaying density
    sum = 0.0
    for n in 1:10
      z = (0.25 + layer_depth*real(n-1))/scale_height
      sum = sum + exp(-z)
    end

    C_atmos  	= csp_atmos*layer_depth*1000.0*rho_atmos*sum/sec_per_yr
    C_soil   	= depth_soil*rho_soil*csp_soil/sec_per_yr 
    C_seaice 	= depth_seaice*rho_sea_ice*csp_sea_ice/sec_per_yr
    C_snow   	= depth_snow * rho_snow * csp_snow/sec_per_yr
    C_mixed_layer = depth_mixed_layer*rho_water*csp_water/sec_per_yr

    # Calculate radiative relaxation times for columns
    tau_land = (C_soil + C_atmos)/B * days_per_yr    
    tau_snow = (C_snow + C_atmos)/B * days_per_yr 
    tau_sea_ice = (C_seaice + C_atmos)/B * days_per_yr  
    tau_mixed_layer = (C_mixed_layer + C_atmos)/B   

    # define heatcap
    heatcap = zeros(size(geography,2),size(geography,1))

    # Assign the correct value of the heat capacity of the columns
    for j in 1:size(geography,1)
      for i in 1:size(geography,2)
        geo  = geography[i,j]
        if geo == 1                            # land
          heatcap[i,j] = C_soil + C_atmos  
        elseif geo == 2                        # perennial sea ice
          heatcap[i,j] = C_seaice + C_atmos
        elseif geo == 3                        # permanent snow cover 
          heatcap[i,j] = C_snow + C_atmos         
        elseif geo == 4                        # lakes, inland seas
          heatcap[i,j] = C_mixed_layer/3.0 + C_atmos 
        elseif geo == 5                        # Pacific ocean 
          heatcap[i,j] = C_mixed_layer + C_atmos
        elseif geo == 6                        # Atlantic ocean 
          heatcap[i,j] = C_mixed_layer + C_atmos
        elseif geo == 7                        # Indian ocean 
          heatcap[i,j] = C_mixed_layer + C_atmos
        elseif geo == 8                        # Mediterranean 
          heatcap[i,j] = C_mixed_layer + C_atmos
        end                           
      end
    end  

    return heatcap, tau_land, tau_snow, tau_sea_ice, tau_mixed_layer
  end


  # Calculate the diffusion coefficients at finest grid level.
  function calc_diffusion_coefficients(geography,nlongitude=128,nlatitude=65)

    Keq     = 0.65 # coefficinet for diffusion at equator
    Kocean  = 0.40 # coefficient for ocean diffusion
    Kland   = 0.65 # coefficinet for land diffusion
    KlandNP = 0.28 # coefficinet for land diffusion (north pole)
    KlandSP = 0.20 # coefficinet for land diffusion (south pole)

    diffusion = zeros(Float64,nlongitude,nlatitude)

    j_equator = div(nlatitude,2) + 1

    for j = 1:nlatitude
        theta = pi*real(j-1)/real(nlatitude-1)
        colat = sin(theta)^5

        for i = 1:nlongitude
            let geo = geography[i,j]
                if geo >= 5 && Geo <= 7 # oceans
                   diffusion[i,j] = (Keq-Kocean)*colat + Kocean
                else # land, sea ice, etc
                    if j <= j_equator # northern hemisphere
                        diffusion[i,j] = (Kland-KlandNP)*colat + KlandNP
                    else # southern hemisphere
                        diffusion[i,j] = (Kland-KlandSP)*colat + KlandSP
                    end
                 end
            end
        end
    end

    return diffusion
  end

  # Calculate the area weighted mean of the diffusion at the mid-point between
  # the pole and the first ring of grid points.

  # If you need the mean diffusion at the poles just do the following:
  #   mean_diffusion_north = sum(diffusion_north)
  #   mean_diffusion_south = sum(diffusion_south)

  function calc_diffusion_coefficients_poles(diffusion,nlongitude=128,nlatitude=65)

    # Fractional areas for the poles
    angle  = pi/real(nlatitude-1)
    area_1 = 0.5*(1.0 - cos(0.5*angle))
    area_2 = sin(0.5*angle)*sin(angle)/float(nlongitude)

    total_area = area_1 + area_2

    diffusion_north = zeros(Float64,nlongitude)
    diffusion_south = zeros(Float64,nlongitude)

    for i = 1:nlongitude
        diffusion_north[i] = (area_1*diffusion[1,        1] + area_2*diffusion[i,          2])/total_area
        diffusion_south[i] = (area_1*diffusion[1,nlatitude] + area_2*diffusion[i,nlatitude-1])/total_area
    end

    return diffusion_north,diffusion_south
  end



  function read_albedo(filepath="albedo.dat",nlongitude=128,nlatitude=65)
    result = zeros(Float64,nlatitude,nlongitude)
    open(filepath) do fh
        for lat = 1:nlatitude
            if eof(fh) break end
            result[lat,:] = parse.(Float64,split(strip(readline(fh) ),r"\s+"))
        end
    end
    return result
  end

 
  function read_world(filepath="The_World.dat",nlongitude=128,nlatitude=65)
    result = zeros(Float64,nlatitude,nlongitude)
    open(filepath) do fh
        for lat = 1:nlatitude
            if eof(fh) break end
            result[lat,:] = parse.(Int8,split(strip(readline(fh) ),r""))
        end
    end
    return result
  end


