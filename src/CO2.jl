using DelimitedFiles
using DataFrames
using Statistics
using CSV

const CSV_DATA = joinpath(@__DIR__, "..", "input", "preprocessed") #change preprocessed to 8_Tage, if you want to use a 
const filepaths = readdir(CSV_DATA, join=true)
CO2_Array = readdlm(joinpath(@__DIR__, "..", "input", "co2_mm_mlo.txt"), comments=true)
const antarctic_CO2 = readdlm(joinpath(@__DIR__, "..", "input", "mlo_spo_monthly_mean.csv"), comments=true)

function read_south_pole_CO2_values()
    (row, column) = size(antarctic_CO2)
    global yearSP = []
    global monthSP = []
    for n in 1:row
        antarctic_line = replace(antarctic_CO2[n, 1], "," => " ")
        z = split(antarctic_line, " ")
        t = parse(Int64, z[1, 1])
        r = parse(Int64, z[2, 1])
        push!(yearSP, t)
        push!(monthSP, r)
    end
    return yearSP, monthSP
end



" 
    The function count_files counts all files in an explicit folder.
        
"


function count_files()
    global counter = 0
    for file in filepaths
        if isfile(file)
            counter += 1
        end
   end
    return counter
end

function max_element_array()
max_val = zeros(counter,2)
count = 0
global max_row_column = 0
for file in filepaths
    if isfile(file)
        count +=1
    end
    read = CSV.read(file,DataFrame)
    (rows,columns) = size(read)
    max_val[count,1] = rows
    max_val[count,2] = columns
end
row_column, pos = findmax(max_val, dims = 1)
max_row_column = Int.(row_column)

return max_row_column
end


"
    The dataset given by NASA uses 2 degree steps for the latitude and 2.5 degree steps for the longitude, so the 
    maximum amout of rows are 90 - NASA uses 91 rows - and the maximum amout of columns are 144.
    The function below reads in the datasets given by NASA and then creates a matrix with size 91x144 (for each month).
    After that we cut out rows and columns so that we can use it. In addition, there is a function called CO2_bilinear_interpolation
    to cut out the right rows and columns, which are not neccessary. Also the function build_CO2_matrix_with_replaced_values
    changes the values, which are obviously incorrect, with average values from co2_mm_mlo.
      
"

num = 0
function build_CO2_matrix_with_replaced_values(num, vector1_in, vector2_in, vector3_in)
    global final_CO2 = []
    global nasa_co2_mat = []
    global index = []
    global CO2_averages = []
    global CO2_averages_SP = []
    global nasa_co2_mat_mod = []
    for path in filepaths
        global num += 1
        open(path) do io
            co2_mat = zeros(vector3_in[1,1], vector3_in[1,2]) #choose a matrix, which is big enough - 
            split_path = splitpath(path)
            filename = split(split_path[length(split_path)], "_")
            year = parse(Int, filename[1])
            push!(index, year)
            month = parse(Int, filename[2])
            indexYear = findall(x -> x == year, CO2_Array[:, 1]) #try to keep x variable, so that it can be changed
            indexMonth = findall(x -> x == month, CO2_Array[:, 2])
            indexYearSP = findall(x -> x == year, vector1_in)
            indexMonthSP = findall(x -> x == month, vector2_in)
            for t in indexYear, s in indexMonth
                if t == s
                    global month_average = CO2_Array[t, 4] * 1 / (10^6)  #change CO2 concentration in co2_mm_mlo to receive different results
                    push!(CO2_averages, month_average)
                    break
                end
            end

            for p in indexYearSP, q in indexMonthSP
                if p == q
                    global month_average_SP = antarctic_CO2[p, 5] * 1 / (10^6)
                    push!(CO2_averages_SP, month_average_SP)
                    break
                end
            end

            data = readdlm(io) #1.Schritt -> 2002.09.01.csv

            for i in 1:vector3_in[1,1]
                vector_line = split(data[i], ",")
                df = DataFrame(y=vector_line)
                df.x = parse.(Float64, df.y)
                matrix_line = permutedims(df.x)
                for j in 1:length(matrix_line)
                    co2_mat[i, j] = matrix_line[1, j]
                end
            end

            push!(nasa_co2_mat, 100000 * co2_mat)

            for n in 1:vector3_in[1,1]
                for k in 1:vector3_in[1,2]
                    if n < 60
                        if (co2_mat[n, k] == -9999.0 || co2_mat[n, k] == 0.0)
                            co2_mat[n, k] = month_average
                        end

                    else
                        if (co2_mat[n, k] == -9999.0 || co2_mat[n, k] == 0.0)
                            co2_mat[n, k] = month_average_SP
                        end
                    end
                end
            end
            co2_mat_t = transpose(co2_mat)
            push!(nasa_co2_mat_mod, co2_mat_t)
        end
    end
    return nasa_co2_mat, nasa_co2_mat_mod, CO2_averages, index, CO2_averages_SP
    #final_CO2 is a matrix with manually deleted rows and columns and wrong values changed with month_average
    #arr_no_change is a matrix with wrong values and without deleted rows and columns
    #arr_change is a matrix with changed values and without deleted rows and columns 
end

count_files()

function calc_only_annual_year(array_in, vector_in)
    global mean_year_co2_mat = []
    global sum_all = []
    for a in 0:13
        y = findfirst(x -> x == (2003 + a), vector_in)
        vec = Vector(0:11)
        mean_sum = (1 / 12) * (array_in[y+vec[1]] .+ array_in[y+vec[2]] .+ array_in[y+vec[3]] .+ array_in[y+vec[4]] .+ array_in[y+vec[5]] .+
                               array_in[y+vec[6]] .+ array_in[y+vec[7]] .+ array_in[y+vec[8]] .+ array_in[y+vec[9]] .+ array_in[y+vec[10]] .+ array_in[y+vec[11]] .+
                               array_in[y+vec[12]])
        push!(sum_all, mean_sum)
    end
    vec1 = Vector(1:14)
    mean_sum_all_years = (1 / 14) * (sum_all[vec1[1]] .+ sum_all[vec1[2]] .+ sum_all[vec1[3]] .+ sum_all[vec1[4]] .+ sum_all[vec1[5]] .+ sum_all[vec1[6]] .+ 
                        sum_all[vec1[7]] .+ sum_all[vec1[8]] .+ sum_all[vec1[9]] .+ sum_all[vec1[10]] .+ sum_all[vec1[11]] .+ sum_all[vec1[12]] .+ sum_all[vec1[13]] .+
                        sum_all[vec1[14]])
    push!(mean_year_co2_mat, mean_sum_all_years)
    return mean_year_co2_mat, sum_all
end

nlongitude = 128

function co2_bilinear_interpolation(array_in, nlongitude, counter)
    global bilinear_co2_mat = []
    for n in 1:counter
        (width_in, height_in) = size(array_in[1])

        width_out = nlongitude #convert(Int64,(nlongitude/2+1))
        height_out = convert(Int64, (nlongitude / 2 + 1))
        array_out = Array{Float64}(ones((width_out, height_out)))

        # calculate the new distance between each dot/pixle (ratio)
        # "-1" because the indicies in Julia
        ratio_w = (width_in - 1) / (width_out - 1)
        ratio_h = (height_in - 1) / (height_out - 1)

        high = 1
        wide = 1

        for w = 1:width_out
            high = 1
            for h = 1:height_out

                # the coordinates from the current pixel/knot
                y = 1 + (h - 1) * ratio_h
                x = 1 + (w - 1) * ratio_w

                # calculate coordinates from surrounding pixels/knots
                x_floor = floor(Int64, x)
                y_floor = floor(Int64, y)
                x_ceil = min(ceil(Int64, x), width_in)
                y_ceil = min(ceil(Int64, y), height_in)

                # get values of the surrounding pixels 
                # case 1: x and y are intergers => x_ceil = x_floor, y_ceil = y_floor     
                if x_ceil == x_floor && y_ceil == y_floor
                    val = array_in[n][x_ceil, y_ceil]

                    # case 2: x is an integer 
                elseif x_ceil == x_floor
                    val1 = array_in[n][x_ceil, y_floor]
                    val2 = array_in[n][x_ceil, y_ceil]
                    val = val1 * (y_ceil - y) + val2 * (y - y_floor)

                    #case 3: y is an integer
                elseif y_ceil == y_floor
                    val1 = array_in[n][x_floor, y_ceil]
                    val2 = array_in[n][x_ceil, y_ceil]
                    val = val1 * (x_ceil - x) + val2 * (x - x_floor)
                    # case 4: x and y are floating point numbers     
                else
                    v1 = array_in[n][x_floor, y_floor]
                    v2 = array_in[n][x_ceil, y_floor]
                    v3 = array_in[n][x_floor, y_ceil]
                    v4 = array_in[n][x_ceil, y_ceil]

                    # calculate the pixels (x, y) value
                    val1 = (v1 * (x_ceil - x) + v2 * (x - x_floor))  # first linear interpolation in x direction between v1 and v2
                    val2 = (v3 * (x_ceil - x) + v4 * (x - x_floor))  # second linear interpolation in x direction between v3 and v4
                    val = val1 * (y_ceil - y) + val2 * (y - y_floor)  # final "bilinear interpolation"
                end

                array_out[wide, high] = round(val, digits=8)
                high += 1
            end
            wide += 1
        end
        push!(bilinear_co2_mat, array_out)
    end
end

function bilinear_one_matrix(array_in, nlongitude)

    (width_in, height_in) = size(array_in)

    width_out = nlongitude #convert(Int64,(nlongitude/2+1))
    height_out = convert(Int64, (nlongitude / 2 + 1))
    global array_out = Array{Float64}(ones((width_out, height_out)))

    # calculate the new distance between each dot/pixle (ratio)
    # "-1" because the indicies in Julia
    ratio_w = (width_in - 1) / (width_out - 1)
    ratio_h = (height_in - 1) / (height_out - 1)

    high = 1
    wide = 1

    for w = 1:width_out
        high = 1
        for h = 1:height_out

            # the coordinates from the current pixel/knot
            y = 1 + (h - 1) * ratio_h
            x = 1 + (w - 1) * ratio_w

            # calculate coordinates from surrounding pixels/knots
            x_floor = floor(Int64, x)
            y_floor = floor(Int64, y)
            x_ceil = min(ceil(Int64, x), width_in)
            y_ceil = min(ceil(Int64, y), height_in)

            # get values of the surrounding pixels 
            # case 1: x and y are intergers => x_ceil = x_floor, y_ceil = y_floor     
            if x_ceil == x_floor && y_ceil == y_floor
                val = array_in[x_ceil, y_ceil]

                # case 2: x is an integer 
            elseif x_ceil == x_floor
                val1 = array_in[x_ceil, y_floor]
                val2 = array_in[x_ceil, y_ceil]
                val = val1 * (y_ceil - y) + val2 * (y - y_floor)

                #case 3: y is an integer
            elseif y_ceil == y_floor
                val1 = array_in[x_floor, y_ceil]
                val2 = array_in[x_ceil, y_ceil]
                val = val1 * (x_ceil - x) + val2 * (x - x_floor)
                # case 4: x and y are floating point numbers     
            else
                v1 = array_in[x_floor, y_floor]
                v2 = array_in[x_ceil, y_floor]
                v3 = array_in[x_floor, y_ceil]
                v4 = array_in[x_ceil, y_ceil]

                # calculate the pixels (x, y) value
                val1 = (v1 * (x_ceil - x) + v2 * (x - x_floor))  # first linear interpolation in x direction between v1 and v2
                val2 = (v3 * (x_ceil - x) + v4 * (x - x_floor))  # second linear interpolation in x direction between v3 and v4
                val = val1 * (y_ceil - y) + val2 * (y - y_floor)  # final "bilinear interpolation"
            end

            array_out[wide, high] = round(val, digits=8)
            high += 1
        end
        wide += 1
    end
    return array_out
end

#etol=(10^(-5))

function bisection(array1, array2, vector1_in, vector2_in, vector3_in, counter; etol=(10^(-6))) #array1 bilinear, array2 big matrix with changed values

    global final_co2_mat = []

    for t in 1:counter
        lower = 0.0002
        upper = 0.000650
        global l = 1
        #Define function(lower) for the if block
        array4 = replace(array2[t], vector2_in[t] => lower, vector3_in[t] => lower)
        sm_ar_4 = bilinear_one_matrix(array4, nlongitude)
        sm_ar_4_sum = sum(sm_ar_4, dims=1)
        sm_ar_4_sum[1] = sm_ar_4_sum[1] ./ 128
        sm_ar_4_sum[65] = sm_ar_4_sum[65] ./ 128
        func1 = dot(sm_ar_4_sum, vector1_in)

        #Define function(upper) to check if function(lower)*function(upper) <0
        array6 = replace(array2[t], vector2_in[t] => upper, vector3_in[t] => upper)
        sm_ar_6 = bilinear_one_matrix(array6, nlongitude)
        sm_ar_6_sum = sum(sm_ar_6, dims=1)
        sm_ar_6_sum[1] = sm_ar_6_sum[1] ./ 128
        sm_ar_6_sum[65] = sm_ar_6_sum[65] ./ 128
        func3 = dot(sm_ar_6_sum, vector1_in)


        #Calculate the first middle value
        section1 = (upper + lower) / 2 

        #Build the co2 matrix with values from section1 from above
        array5 = replace(array2[t], vector2_in[t] => section1, vector3_in[t] => section1)
        sm_ar_5 = bilinear_one_matrix(array5, nlongitude)
        sm_ar_5_sum = sum(sm_ar_5, dims=1)
        sm_ar_5_sum[1] = sm_ar_5_sum[1] ./ 128
        sm_ar_5_sum[65] = sm_ar_5_sum[65] ./ 128
        func2 = dot(sm_ar_5_sum, vector1_in)

        #Check whether already existing co2 matrix is in etol. Therefore calculate the global mean
        z = sum(array1[t], dims=1)
        z[1] = z[1] ./ 128
        z[65] = z[65] ./ 128
        calc_average_CO2 = dot(z, vector1_in)
        y = calc_average_CO2
        println(y)

        #Filter the mean co2 value for the specific year and month
        middle_value = vector2_in[t]
        println(middle_value)

        if (func1 - middle_value)*(func3 - middle_value) > 0
            println("Choose different borders for co2 concentration!")
        end
        func2 = func3
        if abs(y - middle_value) < etol
            println("We use the already calculated co2 matrix!")
        else
            println("We are going to calculate the needed co2 matrix!")
            func_section1 = 1
            while abs(func_section1 - middle_value) > etol
                section1 = (upper+lower)/2
                array5 = replace(array2[t], vector2_in[t] => section1, vector3_in[t] => section1)
                sm_ar_5= bilinear_one_matrix(array5, nlongitude)
                sm_ar_5_sum = sum(sm_ar_5, dims=1)
                sm_ar_5_sum[1] = sm_ar_5_sum[1] ./ 128
                sm_ar_5_sum[65] = sm_ar_5_sum[65] ./ 128
                func_section1 = dot(sm_ar_5_sum, vector1_in)

                if (func2 - vector2_in[t]) * (func_section1 - vector2_in[t]) < 0
                    lower = section1 
                else
                    upper = section1
                end
                l += 1
            end
        end
        if l == 1
            sm_ar_5 = array1[t]
        end
        array_out_1 = 10^6 .* sm_ar_5
        push!(final_co2_mat, array_out_1)
        println(l)
    end
    return final_co2_mat
end


function extend_co2_matrices(mesh, array_in, year_start, year_end, counter)
    @unpack nx, ny = mesh
    global more_years_co2_mat = []
    g = findfirst(CO2_Array[:, 1] .== 2003)
    h = findfirst(CO2_Array[:, 1] .== year_start)

    for i in h:g+1
        z = fill(CO2_Array[i, 4], nx, ny)
        push!(more_years_co2_mat, z)
    end

    for i in 1:counter
        push!(more_years_co2_mat, array_in[i])
    end

    j = findfirst(CO2_Array[:, 1] .== 2017)
    k = findfirst(CO2_Array[:, 1] .== year_end + 1)

    for i in (j+2):k+2
        y = fill(CO2_Array[i, 4], nx, ny)
        push!(more_years_co2_mat, y)
    end

    return more_years_co2_mat
end

# year_end max till 2021 because of missing values for dataset co2_mm_mlo
# year_end max till 2016 for matrices

function calc_ols(mesh, array_in, year_st, year_en, num_steps_years; linear=true)
    @unpack nx, ny = mesh
    (row,) = size(array_in)
    default = 2003
    global beta = zeros(Float64, 2)
    years = year_en - year_st #+1
    month_steps = range(year_st, year_en, years * 12 + 1)  # 12 because we have 12 months in a year
    local co2_vector = []
    global OLS_coefficients = []
    if linear
        n = 1
        if row < 200
            x_matrix = ones(Float64, years * 12 + 1, 2)
            for i in 1:years*12+1
                x_matrix[i, 2] = month_steps[i]
            end
            invers = transpose(x_matrix) * x_matrix
            while n <= nx * ny
                for k in ((year_st-default)*12)+1:12*years+1
                    vector = reshape(array_in[k], nx * ny)
                    push!(co2_vector, vector[n])
                end
                co2_vector_float = convert(Array{Float64,1}, co2_vector)
                beta = invers \ (transpose(x_matrix) * co2_vector_float)
                push!(OLS_coefficients, beta)
                deleteat!(co2_vector, 1:length(co2_vector))
                deleteat!(co2_vector_float, 1:length(co2_vector_float))
                n += 1
            end
        else
            if year_en == 2021
                x_matrix = ones(Float64, years * 12 + 1, 2)
                for i in 1:years*12+1
                    x_matrix[i, 2] = month_steps[i]
                end
                invers = transpose(x_matrix) * x_matrix
                beta = invers \ (transpose(x_matrix) * array_in[1:(years*12+1)])
            else
                x_matrix = ones(Float64, years * 12, 2)
                for i in 1:years*12
                    x_matrix[i, 2] = month_steps[i]
                end
                invers = transpose(x_matrix) * x_matrix
                beta = invers \ (transpose(x_matrix) * array_in[1:years*12])
            end
        end
        #return beta, OLS_coefficients
    else
        n = 1
        if row < 200
            x_matrix = ones(Float64, years * 12 + 1, 2)
            for i in 1:years*12+1
                x_matrix[i, 2] = month_steps[i]
            end
            invers = transpose(x_matrix) * x_matrix
            while n <= nx * ny
                for k in ((year_st-default)*12)+1:12*years+1
                    vector = reshape(array_in[k], nx * ny)
                    push!(co2_vector, vector[n])
                end
                co2_vector_float = convert(Array{Float64,1}, co2_vector)
                beta = invers \ (transpose(x_matrix) * log.(co2_vector_float))
                push!(OLS_coefficients, beta)
                deleteat!(co2_vector, 1:length(co2_vector))
                deleteat!(co2_vector_float, 1:length(co2_vector_float))
                n += 1
            end
        else
            if year_en == 2021
                x_matrix = (ones(Float64, years * 12 - 1, 2))
                for i in 1:years*12-1
                    x_matrix[i, 2] = month_steps[i]
                end
                invers = transpose(x_matrix) * x_matrix
                beta = invers \ (transpose(x_matrix) * log.(array_in[11:(years*12-1)+10]))
            else
                x_matrix = (ones(Float64, years * 12, 2))
                for i in 1:years*12
                    x_matrix[i, 2] = month_steps[i]
                end
                invers = transpose(x_matrix) * x_matrix
                beta = invers \ (transpose(x_matrix) * log.(array_in[11:years*12+10]))
            end
        end
        return beta, OLS_coefficients

    end
end

function calc_co2_matrices(array_in, future_year_str, future_year_en, num_steps_year; linear=true)
    years = future_year_en - future_year_str
    #num_steps_months = num_steps_year/4
    #x = range(future_year_str + 0.25, future_year_en + 0.25, 12 * years + 1)

    global val = []
    (row,) = size(array_in)
    g = zeros(Float64, (128, 65))
    if linear
        if row == 2 && future_year_str == 2022
            x = range(future_year_str, future_year_en, 12 * years + 1)
            f = zeros(Float64, 12 * years + 4)
            f[1] = array_in[1] + array_in[2] * (2021 + 11 / 12)
            push!(val, f[1])
            for i in 2:12*years+2
                f[i] = array_in[1] + array_in[2] * x[i-1]
                push!(val, f[i])
            end
            f[12*years+3] = array_in[1] + array_in[2] * (future_year_en + 1 / 12)
            f[12*years+4] = array_in[1] + array_in[2] * (future_year_en + 2 / 12)
            push!(val, f[12*years+3])
            push!(val, f[12*years+4])


            #elseif row == 2 && future_year_str !== 2022
            #    f = zeros(Float64, num_steps_year * years + 1)
            #    for i in 1:num_steps_year*years+1
            #        f[i] = array_in[1] + array_in[2] * x[i]
            #        push!(val, f[i])
            #    end


        else
            x = range(future_year_str + 2 / 12, future_year_en + 2 / 12, 12 * years + 1)
            f = zeros(Float64, row)
            for i in 1:12*years+1
                for j in 1:row
                    f[j] = array_in[j][1] + array_in[j][2] * x[i]
                end
                g[:, :] = reshape(f, 128, 65)
                push!(val, g[:, :])
            end
        end
    else
        if row == 2 && future_year_str == 2022
            x = range(future_year_str, future_year_en, 12 * years + 1)
            f = zeros(Float64, 12 * years + 4)
            f[1] = exp(array_in[1] + array_in[2] * (2021 + 11 / 12))
            push!(val, f[1])
            for i in 2:12*years+2
                f[i] = exp(array_in[1] + array_in[2] * x[i-1])
                push!(val, f[i])
            end
            f[12*years+3] = exp(array_in[1] + array_in[2] * (future_year_en + 1 / 12))
            f[12*years+4] = exp(array_in[1] + array_in[2] * (future_year_en + 2 / 12))
            push!(val, f[12*years+3])
            push!(val, f[12*years+4])

            #   elseif row == 2 && future_year_str !== 2022
            #       f = zeros(Float64, num_steps_year * years + 1)
            #       for i in 1:num_steps_year*years+1
            #           f[i] = exp(array_in[1] + array_in[2] * x[i])
            #           push!(val, f[i])
            #       end
        else
            x = range(future_year_str + 2 / 12, future_year_en + 2 / 12, 12 * years + 1)
            f = zeros(Float64, row)
            for i in 1:12*years+1
                for j in 1:row
                    f[j] = exp(array_in[j][1] + array_in[j][2] * x[i])
                end
                g[:, :] = reshape(f, 128, 65)
                push!(val, g[:, :])
            end
        end
    end
    return val
end

function combine_matrices(array1_in, array2_in, counter)
    global var_co2_mat = []
    (row,) = size(array2_in)
    mat = isa(array2_in[1], Matrix)
    if mat
        for i in 1:counter
            push!(var_co2_mat, array1_in[i])
        end
        for i in 1:row
            push!(var_co2_mat, array2_in[i])
        end
    else
        for i in 1:length(array1_in)
            push!(var_co2_mat, array1_in[i])
        end
        for i in 1:length(array2_in)
            push!(var_co2_mat, array2_in[i])
        end
        return var_co2_mat
    end
end

function calc_mean_co2(array_in, vector_in)
    (row,) = size(array_in)
    (nx, ny) = size(array_in[1])
    mean_co2 = []
    for i in 1:row
        sum_over_col = sum(array_in[i], dims=1)
        #update the first and last value of the sum because the poles are only a point
        sum_over_col[1] = sum_over_col[1] ./ nx
        sum_over_col[65] = sum_over_col[65] ./ nx
        #calc the mean co2
        mean_val_co2 = dot(sum_over_col, vector_in)
        push!(mean_co2, mean_val_co2)
    end
    return mean_co2
end
