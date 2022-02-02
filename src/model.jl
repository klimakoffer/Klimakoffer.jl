using DelimitedFiles

mutable struct Model
    diffusion_coeff::Array{Float64,2}       # Diffusion coefficient: depends on the sine of latitude
    heat_capacity::Array{Float64,3}         # Heat capacity: depends on the geography (land, ocean, ice, etc.). Position [:,:,1]: Previous time step, position [:,:,2]: current time step.
    albedo::Array{Float64,2}                # Albedo coefficient: depends on the geography (land, ocean, ice, etc.)
    solar_forcing::Array{Float64,3}         # Time-dependent incoming solar radiation: depends on the orbital parameters and the albedo
    radiative_cooling_co2::Array{Float64,1} # Outgoing long-wave radiation: depends on the CO2 concentration. Position 1: Previous time step, position 2: current time step.
    radiative_cooling_feedback::Float64     # Outgoing long-wave radiation (feedback effects): models the water vapor cyces, lapse rate and cloud cover           
    co_albedo::Array{Float64,2}             # Absorbed radiatian coefficient: depends on the albedo
    compute_albedo::Bool                    # Attribute for deciding whether the albedo should be calculated instead of being read from a file
    geography::Array{Int8,2}                # Monthly land-sea-snow/sea ice distribution: depends on initial geography input if compute_sea_ice_extent is true
    solar_irradiance::Array{Float64,2}      # Solar irradiance: depends on the orbital parameters
    compute_sea_ice_extent::Bool            # Attribute for deciding whether sea ice extent should be calculated for the geography or read from given maps
    year::Int64    
  end

  function Model(mesh, num_steps_year; co2_concentration = 315.0, compute_albedo = false, compute_sea_ice_extent = false, sea_ice_extent_year = 1979) # co2_concentration in [ppm], default value from year 1950
    @unpack nx,ny = mesh
    
    # Constants
    radiative_cooling_feedback = 2.15   #[W/m^2/°C]: sensitivity of the seasonal cycle and annual change in the forcing agents
    
    radiative_cooling_co2 = zeros(Float64,2)
    radiative_cooling_co2[1] = calc_radiative_cooling_co2(co2_concentration) # Initialize both entries of radiative_cooling_co2 to the same value
    radiative_cooling_co2[2] = radiative_cooling_co2[1]

    # Read parameters
    
    geography = read_geography(joinpath(@__DIR__, "..", "input", "The_World.dat"),nx,ny)
    
    compute_sea_ice_extent = compute_sea_ice_extent
    year = sea_ice_extent_year

    compute_albedo = compute_albedo
    if compute_albedo == false
      albedo = read_albedo(joinpath(@__DIR__, "..", "input", "albedo.dat"),nx,ny)
    else
  
      albedo = calc_albedo(geography,nx,ny)
    end 


  
    diffusion_coeff = calc_diffusion_coefficients(geography,nx,ny)

    # Compute heat capacity (we initialize the heat capacity in the previous time step as the current heat capacity)
    heat_capacity = zeros(nx,ny,2)
    calc_heat_capacity!(view(heat_capacity,:,:,1),geography,radiative_cooling_feedback) # TODO: Remove unused variables
    heat_capacity[:,:,2] = heat_capacity[:,:,1]

    co_albedo = 1.0 .- albedo
    
    ecc = eccentricity(1950)
    ob = obliquity(1950)
    per = perihelion(1950)
     
    
    solar_irradiance = calc_solar_forcing_params(nlatitude=ny, ntimesteps=num_steps_year, ecc=ecc, ob=ob, per=per)
    solar_forcing = calc_solar_forcing(co_albedo, solar_irradiance, nx, ny, num_steps_year)
  
    
    return Model(diffusion_coeff, heat_capacity, albedo, solar_forcing, radiative_cooling_co2, radiative_cooling_feedback,co_albedo,compute_albedo,geography,solar_irradiance,compute_sea_ice_extent,year)
  end
  
  function set_co2_concentration!(model, co2_concentration)
    model.radiative_cooling_co2[1] = model.radiative_cooling_co2[2]
    model.radiative_cooling_co2[2] = calc_radiative_cooling_co2(co2_concentration)
  end
  
  
 """
  update_model()
  Gives the possibility to interactively decide which parameters should be updated in a simulation
 """
 function update_model(model, mesh, time_step, update_geography=true, update_albedo=true, update_diffusion=true,update_heat_capacity=true, update_solar_forcing=true)
  @unpack nx,ny = mesh
  
  if update_geography == true
    set_geography!(model,mesh, time_step)
  end
  
  if update_albedo == true
    set_albedo!(model,mesh)
    set_co_albedo!(model)
  end
  
  if update_diffusion == true
    set_diffusion_coeff!(model, mesh)
  end
  if update_heat_capacity == true
    set_heat_capacity!(model)
  end
  if update_solar_forcing == true
    set_solar_forcing!(model, mesh)
  end
  end
  
  
  # TODO: this routine works only for NT = 48

  function set_geography!(model, mesh, time_step)
  @unpack nx,ny = mesh

  # monthly geography can be updated by reading the input maps or by computing the sea ice extent
  # first time step is the last week of March as the simulation starts at the vernal equinox
   if model.compute_sea_ice_extent == false
      if time_step in (1,46,47,48)
        model.geography = read_geography(joinpath(@__DIR__, "..", "input","world", "year2004", string("The_World_from_image", nx, "x", ny,"_1")),nx,ny) 
    
      elseif mod(time_step,4) == 2 && time_step in (2:45)
       month::Int64 = 1
        for nt in 1:time_step
          if mod(nt,4) == 1
           month +=1
          end
        end
        model.geography = read_geography(joinpath(@__DIR__, "..", "input","world", "year2004", string("The_World_from_image", nx, "x", ny,"_",month)),nx,ny)
      end
   else  
     annual_sea_ice_extent = read_sea_ice_extent_N(model.year)
     calc_geography_per_month_extent!(model.geography,annual_sea_ice_extent, time_step, nx, ny)
   end   
  end

  
  function set_albedo!(model,mesh)
  @unpack nx,ny = mesh
  if model.compute_albedo == false
    model.albedo = read_albedo(joinpath(@__DIR__, "..", "input", "albedo.dat"),nx,ny)
  else
    model.albedo = calc_albedo(model.geography,nx,ny)
  end 
  end
  
  function set_co_albedo!(model)
  model.co_albedo = 1.0 .- model.albedo
  end
  
  function set_diffusion_coeff!(model, mesh)
  @unpack nx,ny = mesh
  model.diffusion_coeff = calc_diffusion_coefficients(model.geography,nx,ny)
  end
  
  function set_heat_capacity!(model)
  model.heat_capacity[:,:,1] = model.heat_capacity[:,:,2]
  calc_heat_capacity!(view(model.heat_capacity,:,:,2), model.geography, model.radiative_cooling_feedback)
  end
  
  function set_solar_forcing!(model, mesh)
  @unpack nx,ny = mesh
  solar_irradiance = model.solar_irradiance
  model.solar_forcing = calc_solar_forcing(model.co_albedo, solar_irradiance, nx, ny)
  end
  
  
  Base.size(model::Model) = size(model.heat_capacity)
  
  function Base.show(io::IO, model::Model)
  nx, ny = size(model)
  print(io, "Model() with ", nx, "×", ny, " degrees of freedom")
  end
  
  
  """
  calc_radiative_cooling_co2()
  Computes the CO2 parameter depending on the CO2 concentration
  
  * Default CO2 concentration is 315 ppm (equivalent to year 1950)
  """
  function calc_radiative_cooling_co2(co2_concentration=315.0)
  
  # Define base values for co2_concentration and radiative_cooling_co2
  co2_concentration_base = 315.0
  radiative_cooling_co2_base = 210.3
  
  # Doesn't change as long as CP2ppm doesn't change
  radiative_cooling_co2=radiative_cooling_co2_base-5.35*log(co2_concentration/co2_concentration_base)
  
  return radiative_cooling_co2
  end
  
  """
  _calc_insolation()
  auxiliar function of calc_solar_forcing
  """
  function _calc_insolation(dt, ob, ecc, per, nt, nlat, siny, cosy, tany, s0)
  
  eccfac = 1.0 - ecc^2
  rzero  = (2.0*pi)/eccfac^1.5
  
  #  Solve the orbital equation for lambda as a function of time with
  #  a fourth-order Runge-Kutta method
  lambda = zeros(Float64,nt+1)
  solar = zeros(Float64,nlat,nt)
  
  for n in 2:nt+1
    nu = lambda[n-1] - per
    t1 = dt*(rzero*(1.0 - ecc * cos(nu))^2)
    t2 = dt*(rzero*(1.0 - ecc * cos(nu+0.5*t1))^2)
    t3 = dt*(rzero*(1.0 - ecc * cos(nu+0.5*t2))^2)
    t4 = dt*(rzero*(1.0 - ecc * cos(nu + t3))^2)
    lambda[n] = lambda[n-1] + (t1 + 2.0*t2 + 2.0*t3 + t4)/6.0
  end
  
  #  Compute the average daily solar irradiance as a function of
  #  latitude and longitude (time)
  
  for n in 1:nt
    nu = lambda[n] - per
    rhofac = ((1.0- ecc*cos(nu))/eccfac)^2
    sindec = sin(ob)*sin(lambda[n])
    cosdec = sqrt(1.0-sindec^2)
    tandec = sindec/cosdec
    for j in 1:nlat
      z = -tany[j]*tandec
      if z >= 1.0    # polar latitudes when there is no sunrise (winter)
        solar[j,n] = 0.0
      else
        if z <= -1.0                # when there is no sunset (summer)
          solar[j,n] = rhofac * s0 * siny[j] * sindec
        else
          h_zero = acos(z)
          solar[j,n] = rhofac/pi*s0*(h_zero*siny[j]*sindec+cosy[j]*cosdec*sin(h_zero))
        end
      end
    end
  end
  
  return lambda, solar
  end
  
  """
  calc_solar_forcing_params()
  
  * Default s0 is 1371.685 [W/m²] (Current solar constant)
  * Default orbital parameters of correspond to year 1950 AD:
    ecc = 0.016740             
    ob  = 0.409253
    per = 1.783037  
  """
  
  function calc_solar_forcing_params(yr=0; nlatitude=65, ntimesteps=48, solar_cycle=false, s0=1371.685, orbital=false, ecc=0.016740, ob=0.409253, per=1.783037)
  # Calculate the sin, cos, and tan of the latitudes of Earth from the
  # colatitudes, calculate the insolation
  
  dy = pi/(nlatitude-1.0)
  dt = 1.0 / ntimesteps
  
  siny = zeros(Float64,nlatitude)
  cosy = zeros(Float64,nlatitude)
  tany = zeros(Float64,nlatitude)
  
  if yr == 0
    for j in 1:nlatitude
      lat = pi/2.0 - dy*(j-1)    # latitude in radians
      siny[j] = sin(lat)
      if j == 1
        cosy[j] = 0.0
        tany[j] = 1000.0
      elseif j == nlatitude
        cosy[j] = 0.0
        tany[j] = -1000.0
      else
        cosy[j] = cos(lat)
        tany[j] = tan(lat)
      end
    end
    lambda, solar = _calc_insolation(dt, ob, ecc, per, ntimesteps, nlatitude, siny, cosy, tany, s0)
  elseif yr > 0 && (solar_cycle || orbital)
    lambda, solar = _calc_insolation(dt, ob, ecc, per, ntimesteps, nlatitude, siny, cosy, tany, s0)
  end
  
  return solar
  
  end
  
  
  function calc_solar_forcing(co_albedo, solar, nlongitude=128, nlatitude=65, ntimesteps=48)
  
  solar_forcing = zeros(Float64,nlongitude,nlatitude,ntimesteps)
  # calcualte the seasonal forcing
  for ts in 1:ntimesteps
    for j in 1:nlatitude
      for i in 1:nlongitude
        solar_forcing[i,j,ts] = solar[j,ts]*co_albedo[i,j]
      end
    end
  end
  return solar_forcing
  
  end
  
  function calc_heat_capacity!(heatcap,geography,radiative_cooling_feedback=2.15)
  
  # Depths 
  depth_atmos = 5000.         # meters #TODO: not used???
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
  
  c_atmos  	= csp_atmos*layer_depth*1000.0*rho_atmos*sum/sec_per_yr
  c_soil   	= depth_soil*rho_soil*csp_soil/sec_per_yr 
  c_seaice 	= depth_seaice*rho_sea_ice*csp_sea_ice/sec_per_yr
  c_snow   	= depth_snow * rho_snow * csp_snow/sec_per_yr
  c_mixed_layer = depth_mixed_layer*rho_water*csp_water/sec_per_yr
  
  # Assign the correct value of the heat capacity of the columns
  for j in 1:size(geography,2)
    for i in 1:size(geography,1)
      geo  = geography[i,j]
      if geo == 1                            # land
        heatcap[i,j] = c_soil + c_atmos  
      elseif geo == 2                        # perennial sea ice
        heatcap[i,j] = c_seaice + c_atmos
      elseif geo == 3                        # permanent snow cover 
        heatcap[i,j] = c_snow + c_atmos         
      elseif geo == 4                        # lakes, inland seas
        heatcap[i,j] = c_mixed_layer/3.0 + c_atmos 
      elseif geo in (5, 6, 7, 8)
        # 5: Pacific ocean
        # 6: Atlantic ocean
        # 7: Indian ocean
        # 8: Mediterranean
        heatcap[i,j] = c_mixed_layer + c_atmos
      end                           
    end
  end  
  
  end
  
  """
  calc_diffusion_coefficients()
  Calculate the diffusion coefficients at finest grid level.
  """
  function calc_diffusion_coefficients(geography,nlongitude=128,nlatitude=65)
  
  coeff_eq     = 0.65 # coefficinet for diffusion at equator
  coeff_ocean  = 0.40 # coefficient for ocean diffusion
  coeff_land   = 0.65 # coefficinet for land diffusion
  coeff_landNP = 0.28 # coefficinet for land diffusion (north pole)
  coeff_landSP = 0.20 # coefficinet for land diffusion (south pole)
  
  diffusion = zeros(Float64,nlongitude,nlatitude)
  
  j_equator = div(nlatitude,2) + 1
  
  for j = 1:nlatitude
      theta = pi*real(j-1)/real(nlatitude-1)
      colat = sin(theta)^5
  
      for i = 1:nlongitude
          let geo = geography[i,j]
              if geo >= 5 && geo <= 7 # oceans
                  diffusion[i,j] = (coeff_eq-coeff_ocean)*colat + coeff_ocean
              else # land, sea ice, etc
                  if j <= j_equator # northern hemisphere
                      diffusion[i,j] = (coeff_land-coeff_landNP)*colat + coeff_landNP
                  else # southern hemisphere
                      diffusion[i,j] = (coeff_land-coeff_landSP)*colat + coeff_landSP
                  end
                end
          end
      end
  end
  
  return diffusion
  end
  
  function read_albedo(filepath="albedo.dat",nlongitude=128,nlatitude=65)
  result = zeros(Float64,nlongitude,nlatitude)
  open(filepath) do fh
      for lat = 1:nlatitude
          if eof(fh) break end
         result[:,lat] = parse.(Float64,split(strip(readline(fh) ),r"\s+"))
     end
  end
  return result
  end
  
  #Calculate the albedo depending on the geography

  function calc_albedo(geography,nlongitude=128,nlatitude=65)
  result = zeros(Float64,nlongitude,nlatitude)
  dtheta= pi/(nlatitude-1.0)
  
  for j in 1:nlatitude
    theta = 0.5*pi-dtheta*(j-1)
    sintheta = sin(theta)
    legendrepol = (3*sintheta^2-1)*0.5
    for i in 1:nlongitude
      let geo = geography[i,j]
        if geo == 1 # land
          result[i,j] = 0.3+0.12*legendrepol
        elseif geo == 2 # sea ice
          result[i,j] = 0.6
        elseif geo == 3 # snow cover
          result[i,j] = 0.75
        elseif geo in (5,6,7,8) # oceans
          result[i,j] = 0.29+0.12*legendrepol
        end
      end
    end
  end
  return result
  end
  
  function read_geography(filepath="The_World.dat",nlongitude=128,nlatitude=65)
  result = zeros(Int8,nlongitude,nlatitude)
  open(filepath) do fh
      for lat = 1:nlatitude
          if eof(fh) break end
          result[:,lat] = parse.(Int8,split(strip(readline(fh) ),r""))
      end
  end
  return result
  end
  
  
  # Table 4 from Berger 1978, Long-Term Variations of daily Insolation and Quaternary Climatic Changes
  # Higher accuracy values from Gary Russel, https://data.giss.nasa.gov/modelE/ar5plots/srorbpar.html
  const eccentricity_parameters_amplitude = [0.01860798, 0.01627522, -0.01300660, 0.00988829,
                                           -0.00336700, 0.00333077, -0.00235400, 0.00140015,
                                           0.00100700, 0.00085700, 0.00064990, 0.00059900,
                                           0.00037800, -0.00033700, 0.00027600, 0.00018200,
                                           -0.00017400, -0.00012400, 0.00001250]
  const eccentricity_parameters_mean_rate = [4.2072050, 7.3460910, 17.8572630, 17.2205460, 16.8467330,
                                           5.1990790, 18.2310760, 26.2167580, 6.3591690, 16.2100160,
                                           3.0651810, 16.5838290, 18.4939800, 6.1909530, 18.8677930,
                                           17.4255670, 6.1860010, 18.4174410, 0.6678630]
  const eccentricity_parameters_phase = [28.620089, 193.788772, 308.307024, 320.199637, 279.376984,
                                       87.195000, 349.129677, 128.443387, 154.143880, 291.269597,
                                       114.860583, 332.092251, 296.414411, 145.769910, 337.237063,
                                       152.092288, 126.839891, 210.667199, 72.108838]
  
  # Table 1 from Berger 1978, Long-Term Variations of daily Insolation and Quaternary Climatic Changes
  # Higher accuracy values from Gary Russel, https://data.giss.nasa.gov/modelE/ar5plots/srorbpar.html
  const obliquity_parameters_amplitude = [-2462.2214466, -857.3232075, -629.3231835, -414.2804924,
                                        -311.7632587, 308.9408604, -162.5533601, -116.1077911,
                                        101.1189923, -67.6856209, 24.9079067, 22.5811241,
                                        -21.1648355, -15.6549876, 15.3936813, 14.6660938,
                                        -11.7273029, 10.2742696, 6.4914588, 5.8539148, -5.4872205,
                                        -5.4290191, 5.160957, 5.0786314, -4.0735782, 3.7227167,
                                        3.3971932, -2.8347004, -2.6550721, -2.5717867, -2.4712188,
                                        2.462541, 2.2464112, -2.0755511, -1.9713669, -1.8813061,
                                        -1.8468785, 1.8186742, 1.7601888, -1.5428851, 1.4738838,
                                        -1.4593669, 1.4192259, -1.181898, 1.1756474, -1.1316126,
                                        1.0896928]
  const obliquity_parameters_mean_rate = [31.609974, 32.620504, 24.172203, 31.983787, 44.828336,
                                        30.973257, 43.668246, 32.246691, 30.599444, 42.681324,
                                        43.836462, 47.439436, 63.219948, 64.230478, 1.01053,
                                        7.437771, 55.782177, 0.373813, 13.218362, 62.583231,
                                        63.593761, 76.43831, 45.815258, 8.448301, 56.792707,
                                        49.747842, 12.058272, 75.27822, 65.241008, 64.604291,
                                        1.647247, 7.811584, 12.207832, 63.856665, 56.15599,
                                        77.44884, 6.801054, 62.209418, 20.656133, 48.344406,
                                        55.14546, 69.000539, 11.07135, 74.291298, 11.047742,
                                        0.636717, 12.844549]
  const obliquity_parameters_phase = [251.9025, 280.8325, 128.3057, 292.7252, 15.3747, 263.7951,
                                    308.4258, 240.0099, 222.9725, 268.7809, 316.7998, 319.6024,
                                    143.805, 172.7351, 28.93, 123.5968, 20.2082, 40.8226, 123.4722,
                                    155.6977, 184.6277, 267.2772, 55.0196, 152.5268, 49.1382,
                                    204.6609, 56.5233, 200.3284, 201.6651, 213.5577, 17.0374,
                                    164.4194, 94.5422, 131.9124, 61.0309, 296.2073, 135.4894,
                                    114.875, 247.0691, 256.6114, 32.1008, 143.6804, 16.8784,
                                    160.6835, 27.5932, 348.1074, 82.6496]
  
  # Table 5 from Berger 1978, Long-Term Variations of daily Insolation and Quaternary Climatic Changes
  # Higher accuracy values from Gary Russel, https://data.giss.nasa.gov/modelE/ar5plots/srorbpar.html
  const perihelion_parameters_amplitude = [7391.022589, 2555.1526947, 2022.7629188, -1973.6517951,
                                         1240.2321818, 953.8679112, -931.7537108, 872.3795383,
                                         606.3544732, -496.0274038, 456.9608039, 346.946232,
                                         -305.8412902, 249.6173246, -199.10272, 191.0560889,
                                         -175.2936572, 165.9068833, 161.1285917, 139.7878093,
                                         -133.5228399, 117.0673811, 104.6907281, 95.3227476,
                                         86.7824524, 86.0857729, 70.5893698, -69.9719343,
                                         -62.5817473, 61.5450059, -57.9364011, 57.1899832,
                                         -57.0236109, -54.2119253, 53.2834147, 52.1223575,
                                         -49.0059908, -48.3118757, -45.4191685, -42.235792,
                                         -34.7971099, 34.4623613, -33.8356643, 33.6689362,
                                         -31.2521586, -30.8798701, 28.4640769, -27.1960802,
                                         27.0860736, -26.3437456, 24.725374, 24.6732126, 24.4272733,
                                         24.0127327, 21.7150294, -21.5375347, 18.1148363,
                                         -16.9603104, -16.1765215, 15.5567653, 15.4846529,
                                         15.2150632, 14.5047426, -14.3873316, 13.1351419,
                                         12.8776311, 11.9867234, 11.9385578, 11.7030822, 11.6018181,
                                         -11.2617293, -10.4664199, 10.433397, -10.2377466,
                                         10.1934446, -10.1280191, 10.0289441, -10.0034259]
  const perihelion_parameters_mean_rate = [31.609974, 32.620504, 24.172203, 0.636717, 31.983787,
                                         3.138886, 30.973257, 44.828336, 0.991874, 0.373813,
                                         43.668246, 32.246691, 30.599444, 2.147012, 10.511172,
                                         42.681324, 13.650058, 0.986922, 9.874455, 13.013341,
                                         0.262904, 0.004952, 1.142024, 63.219948, 0.205021,
                                         2.151964, 64.230478, 43.836462, 47.439436, 1.384343,
                                         7.437771, 18.829299, 9.500642, 0.431696, 1.16009,
                                         55.782177, 12.639528, 1.155138, 0.168216, 1.647247,
                                         10.884985, 5.610937, 12.658184, 1.01053, 1.983748,
                                         14.023871, 0.560178, 1.273434, 12.021467, 62.583231,
                                         63.593761, 76.43831, 4.28091, 13.218362, 17.818769,
                                         8.359495, 56.792707, 8.448301, 1.978796, 8.863925,
                                         0.186365, 8.996212, 6.771027, 45.815258, 12.002811,
                                         75.27822, 65.241008, 18.870667, 22.009553, 64.604291,
                                         11.498094, 0.578834, 9.237738, 49.747842, 2.147012,
                                         1.196895, 2.133898, 0.173168]
  const perihelion_parameters_phase = [251.9025, 280.8325, 128.3057, 348.1074, 292.7252, 165.1686,
                                     263.7951, 15.3747, 58.5749, 40.8226, 308.4258, 240.0099,
                                     222.9725, 106.5937, 114.5182, 268.7809, 279.6869, 39.6448,
                                     126.4108, 291.5795, 307.2848, 18.93, 273.7596, 143.805,
                                     191.8927, 125.5237, 172.7351, 316.7998, 319.6024, 69.7526,
                                     123.5968, 217.6432, 85.5882, 156.2147, 66.9489, 20.2082,
                                     250.7568, 48.0188, 8.3739, 17.0374, 155.3409, 94.1709, 221.112,
                                     28.93, 117.1498, 320.5095, 262.3602, 336.2148, 233.0046,
                                     155.6977, 184.6277, 267.2772, 78.9281, 123.4722, 188.7132,
                                     180.1364, 49.1382, 152.5268, 98.2198, 97.4808, 221.5376,
                                     168.2438, 161.1199, 55.0196, 262.6495, 200.3284, 201.6651,
                                     294.6547, 99.8233, 213.5577, 154.1631, 232.7153, 138.3034,
                                     204.6609, 106.5938, 250.4676, 332.3345, 27.3039]
  
  
  
  # Eqn. 4 from Berger 1978, Long-Term Variations of daily Insolation and Quaternary Climatic Changes
  function e_sincos_pi(year)
  M    = eccentricity_parameters_amplitude
  g    = eccentricity_parameters_mean_rate
  beta = eccentricity_parameters_phase
  
  arg = deg2rad.(g/3600) * (year-1950) .+ deg2rad.(beta)
  
  return sum(M .* sin.(arg)), sum(M .* cos.(arg))
  end
  
  eccentricity(year, e_sin_pi, e_cos_pi) = sqrt(e_sin_pi^2 + e_cos_pi^2)
  
  eccentricity(year) = eccentricity(year, e_sincos_pi(year)...)
  
  
  # Eqn. 1 from Berger 1978, Long-Term Variations of daily Insolation and Quaternary Climatic Changes
  function obliquity(year)
  A     = obliquity_parameters_amplitude
  f     = obliquity_parameters_mean_rate
  delta = obliquity_parameters_phase
  
  epsilon_star = 23.320556
  
  return deg2rad(epsilon_star + sum(A/3600 .* cos.(deg2rad.(f/3600) * (year-1950) .+ deg2rad.(delta))))
  end
  
  
  # Eqn. 6 from Berger 1978, Long-Term Variations of daily Insolation and Quaternary Climatic Changes
  function perihelion(year)
  F     = perihelion_parameters_amplitude
  f     = perihelion_parameters_mean_rate
  delta = perihelion_parameters_phase
  
  pi_ = atan(e_sincos_pi(year)...)
  psi_tilde = 50.439273
  zeta = 3.392506
  
  psi = psi_tilde/3600 * (year-1950) + zeta + sum(F/3600 .* sin.(deg2rad.(f/3600) * (year-1950) .+ deg2rad.(delta)))
  
  return (pi_ + deg2rad(psi)) % (2*pi)
  end


"""
read_sea_ice_extent_N()
* Default year: 1979
* Possible years to pick annual sea ice extent for the Northern Hemisphere : 1979 to 2020
* returns a vector with sea ice extent per month from March to February

"""
function read_sea_ice_extent_N(year=1979)

  matrix = readdlm(joinpath(@__DIR__,"..","input/sea_ice_data/","N_ALL_extent_v3.0.dat"),comments=true)
  yearindex = findfirst(isequal(year),matrix[:,1])
  annual_sea_ice_extent = matrix[yearindex,2:13]

  return annual_sea_ice_extent
end



"""
calc_geography_per_month_extent!()
* This routine a new geography matrix by determining sea ice extent with a monthly rate of change
  and mapping extent or reduction to adjacent sea ice cells and ocean cells on the Northern Hemisphere.

* This routine assumes that the sea ice extent of the geography file corresponds to 
  the measured sea ice at the vernal equinox.
 
* There is a bias in the order sea ice cells are added/removed.
  sea ice extent: starts from top left to bottom right
  sea ice reduction: starts from bottom left to top right

"""
function calc_geography_per_month_extent!(geography,annual_sea_ice_extent, time_step, nlongitude=128, nlatitude =65)

  month::Int64 = 0


 # Gets the right month when the simulation starts at the vernal equinox
 # TODO: this part works only for NT = 48

 if time_step in (1,46,47,48)
  month = 1
elseif mod(time_step,4) == 2
    month = 1
    for nt in 1:time_step
      if mod(nt,4) == 1
        month +=1
      end
    end   
 end

  nlat::Int64 = (nlatitude-1)*0.5
  geography_start = read_geography(joinpath(@__DIR__, "..", "input", "The_World.dat"),128,65)
  geo = geography_start[:,1:nlat]

  if month == 0 || month == 1
    geography = geography_start
    return geography
  else
 

  sea_ice_index_unit = size(findall(isequal(2),geo),1)/annual_sea_ice_extent[1]
  
  

   # collect submatrix and indices for the Northern Hemisphere

   G = geography[:,1:nlat]
   N = reverse(G, dims = 2)
   IndG = CartesianIndices(G)
   IndN = CartesianIndices(N)
   Gfirst, Glast = first(IndG), last(IndG)
   Ifirst, Ilast = first(IndN), last(IndN)
   G1 = oneunit(Gfirst)
   I1 = oneunit(Ifirst)


   # get all sea ice cells for the Northern Hemispehre and determine how many cells will be affected by the extent/reduction

   sea_ice_index = findall(isequal(2),G)
   old_counter = size(sea_ice_index,1)
   new_counter = sea_ice_index_unit*annual_sea_ice_extent[month]
   sea_ice_cells_delta = abs(floor(Int,new_counter - old_counter))


   # case for an extent of sea ice
   # adjacent ocean cells (geo=5,6,7,8) will turn to sea ice cells
   if sea_ice_cells_delta == 0
     return geography
   end
   if new_counter > old_counter
      for I in IndG
          if G[I[1], I[2]] == 2
              if sea_ice_cells_delta <= 0
                  break
              else
                  for J in max(Gfirst, I-G1):min(Glast, I+G1)
                      if G[J] in (5,6,7,8)
                          G[J[1], J[2]] = 2
                          sea_ice_cells_delta -= 1
                         if sea_ice_cells_delta <= 0
                             break
                         end
                      end
                  end
              end
          end
      end
      geography[:,1:nlat] = G

   # case for a reduction of sea ice
   # adjacent sea ice cells will turn to ocean cells
   elseif new_counter < old_counter
      for I in IndN
          if N[I[1], I[2]] == 2 && sea_ice_cells_delta > 0
              N[I[1], I[2]] = 5
              sea_ice_cells_delta -= 1
              if sea_ice_cells_delta <= 0
                  break
              else
                  for J in max(Ifirst, I-I1):min(Ilast, I+I1)
                      if N[J] == 2
                          N[J[1], J[2]] = 5
                          sea_ice_cells_delta -= 1
                          if sea_ice_cells_delta <= 0
                              break
                          end
                      end
                  end
  
              end  
              
          end
      end 
      geography[:,1:nlat] = reverse(N, dims = 2)
   end
 
   return geography
  end

end    
