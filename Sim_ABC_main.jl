include("./Sim_ABC_functions.jl")


#region, Step 1: Simulating of the pseudo-data (from the model)

using Distributed

n_proc = 10
addprocs(10, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end

@everywhere include("./Drylands_ABC_functions.jl")

@everywhere function Run_sim_model_EWS(N_sim)

    pseudo_param = CSV.read("./Data/Pseudo_parameters_Eby.csv", DataFrame, header=1, delim=';')[((N_sim-1)*200+1):(N_sim*200), :]

    param = Get_classical_param_Eby()
    fraction_cover = [0.8, 0.2]
    size_landscape = 50
    ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros(200, size(pseudo_param)[2] + n_metric)
    summary_stat_table[:, 1:2] .= pseudo_param

    @inbounds for i in 1:size(pseudo_param)[1]
        param[1:2] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
            keep_landscape=false)
        if d[size(d)[1], 1] < 0.9 #we don't keep very large cover not useful, not representative of data

            d, land = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
                keep_landscape=true, burning_phase=1000)

            if d[size(d)[1], 1] > 0.05 && d[size(d)[1], 1] < 1 #we only compute sumstat when we are sure of the outcome
                summary_stat_table[i, 3:size(summary_stat_table)[2]] = Get_summary_stat(land)
            else
                sumstat_eby2[index, :] .= zeros(9, 1)
            end
        end
    end

    CSV.write("../Data/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end

pmap(Run_sim_model_EWS, 1:1000)







#endregion
#region, Step 2: Getting the summary statistics of the empirical sites

list_land = readdir("./Data/landscapes/")

summary_stat = zeros(length(list_land), 9)


index = 1
for file_land in list_land
    print(index)
    summary_stat[index, :] = Get_summary_stat(Matrix{Int64}(CSV.read("./Data/Data_Biocom/" * file_land, DataFrame, header=0)))
    index += 1
end

CSV.write("./Data/Data_Biocom/Summary_stats_data.csv", Tables.table(summary_stat), writeheader=false)



#endregion
#region, Step 3: Spatial resolution in simulations



using Distributed



addprocs(15, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end


@everywhere include("./Drylands_ABC_functions.jl")


@everywhere function Run_sim_model_EWS(N_sim)

    n_param = 2
    pooling_vec = [0.25 1 / 3 0.5 1]
    n_keep = 100
    pseudo_param = Matrix{Float64}(CSV.read("./Data/Spatial_resolution/Parameter_Eby.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :])
    param = Get_classical_param_Eby()
    fraction_cover = [0.8, 0.2]
    size_landscape = 50
    ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros((length(pooling_vec) + 1) * n_keep, size(pseudo_param)[2] + n_metric) #for each site = 

    for n_ in 1:n_keep
        for k in 1:n_param
            summary_stat_table[((n_-1)*length(pooling_vec)+1):(n_*length(pooling_vec)), k] .= pseudo_param[n_, k]
        end
    end

    for i in 1:size(pseudo_param)[1]
        param[1:2] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land), keep_landscape=false, n_snapshot=15, n_time_bw_snap=50)

        if d[200, 1] < 0.05 || d[200, 1] > 0.9 #in case there are missing issues
            for x in 1:length(pooling_vec)
                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = zeros(9)
            end
        else
            d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true, n_snapshot=15, n_time_bw_snap=50)

            for (x, pooling_intensity) in enumerate(pooling_vec)

                if pooling_intensity < 1
                    land_pooled = zeros(Int(1 / pooling_intensity) * size(land)[1], Int(1 / pooling_intensity) * size(land)[2], size(land)[3])
                elseif pooling_intensity == 1
                    land_pooled = copy(land)
                else
                    land_pooled = zeros(convert(Int64, floor(size(land)[1] / pooling_intensity)), convert(Int64, floor(size(land)[1] / pooling_intensity)), size(land)[3])
                end

                for third_sim in 1:size(land)[3] #we pool all saved landscapes 
                    if pooling_intensity < 1
                        land_pooled[:, :, third_sim] = inverse_pooling(land[:, :, third_sim], Int(1 / pooling_intensity)) #increasing resolution
                    else
                        land_pooled[:, :, third_sim] = pooling(land[:, :, third_sim], Int(pooling_intensity)) #decreasing resolution
                    end

                end
                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled)

            end
        end
    end

    CSV.write("./Data/Spatial_resolution/Simu/Eby/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:40)






@everywhere function Run_sim_model_EWS(N_sim)

    n_param = 2
    pooling_vec = [0.25 1 / 3 0.5 1]
    n_keep = 100
    pseudo_param = CSV.read("./Data/Spatial_resolution/Parameter_Eby_feedback.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :]

    param = Get_classical_param_Eby()
    fraction_cover = [0.8, 0.2]
    size_landscape = 50
    ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros((length(pooling_vec) + 1) * n_keep, size(pseudo_param)[2] + n_metric) #for each site = 

    for n_ in 1:n_keep
        for k in 1:n_param
            summary_stat_table[((n_-1)*length(pooling_vec)+1):(n_*length(pooling_vec)), k] .= pseudo_param[n_, k]
        end
    end


    for i in 1:size(pseudo_param)[1]

        param[1:2] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land), keep_landscape=false,
            n_snapshot=15, n_time_bw_snap=50, intensity_feedback=6)

        if d[200, 1] < 0.05 || d[200, 1] > 0.9 #in case there are missing issues
            for x in 1:length(pooling_vec)
                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = zeros(9)
            end
        else
            d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true,
                n_snapshot=15, n_time_bw_snap=50, intensity_feedback=6)

            for (x, pooling_intensity) in enumerate(pooling_vec)

                if pooling_intensity < 1
                    land_pooled = zeros(Int(1 / pooling_intensity) * size(land)[1], Int(1 / pooling_intensity) * size(land)[2], size(land)[3])
                elseif pooling_intensity == 1
                    land_pooled = copy(land)
                else
                    land_pooled = zeros(convert(Int64, floor(size(land)[1] / pooling_intensity)), convert(Int64, floor(size(land)[1] / pooling_intensity)), size(land)[3])
                end

                for third_sim in 1:size(land)[3] #we pool all saved landscapes 

                    if pooling_intensity < 1
                        land_pooled[:, :, third_sim] = inverse_pooling(land[:, :, third_sim], Int(1 / pooling_intensity)) #increasing resolution
                    else
                        land_pooled[:, :, third_sim] = pooling(land[:, :, third_sim], Int(pooling_intensity)) #decreasing resolution
                    end

                end

                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled)

            end
        end
    end

    CSV.write("./Data/Spatial_resolution/Simu/Eby_feedback/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:40)








@everywhere function Run_sim_model_EWS(N_sim)

    n_param = 7
    pooling_vec = [0.25 1 / 3 0.5 1]
    n_keep = 100
    pseudo_param = CSV.read("./Data/Spatial_resolution/Parameter_Kefi.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :]

    param = Get_classical_param()
    fraction_cover = [0.8, 0.1, 0.1]
    size_landscape = 50
    ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros((length(pooling_vec) + 1) * n_keep, size(pseudo_param)[2] + n_metric) #for each site = 

    for n_ in 1:n_keep
        for k in 1:n_param
            summary_stat_table[((n_-1)*length(pooling_vec)+1):(n_*length(pooling_vec)), k] .= pseudo_param[n_, k]
        end
    end

    for i in 1:size(pseudo_param)[1]
        param[1:7] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_drylands(time=200, param=copy(param), landscape=copy(ini_land), keep_landscape=false, n_snapshot=15, n_time_bw_snap=50)


        if d[200, 2] < 0.05 || d[200, 2] > 0.9 #in case there are missing issues
            for x in 1:length(pooling_vec)
                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = zeros(9)
            end
        else

            d, land = IBM_drylands(time=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true, n_snapshot=15, n_time_bw_snap=50)

            for (x, pooling_intensity) in enumerate(pooling_vec)

                if pooling_intensity < 1
                    land_pooled = zeros(Int(1 / pooling_intensity) * size(land)[1], Int(1 / pooling_intensity) * size(land)[2], size(land)[3])
                elseif pooling_intensity == 1
                    land_pooled = copy(land)
                else
                    land_pooled = zeros(convert(Int64, floor(size(land)[1] / pooling_intensity)), convert(Int64, floor(size(land)[1] / pooling_intensity)), size(land)[3])
                end


                for third_sim in 1:size(land)[3] #we pool all saved landscapes 

                    if pooling_intensity < 1
                        land_pooled[:, :, third_sim] = inverse_pooling(land[:, :, third_sim], Int(1 / pooling_intensity)) #increasing resolution
                    else
                        land_pooled[:, :, third_sim] = pooling(land[:, :, third_sim], Int(pooling_intensity)) #decreasing resolution
                    end

                end

                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled)

            end
        end
    end

    CSV.write("./Data/Spatial_resolution/Simu/Kefi/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)




end


pmap(Run_sim_model_EWS, 1:40)






@everywhere function Run_sim_model_EWS(N_sim)

    n_param = 8
    pooling_vec = [0.25 1 / 3 0.5 1]
    n_keep = 100
    pseudo_param = CSV.read("./Data/Spatial_resolution/Parameter_Schneider.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :]

    param = Get_classical_param_Schneider()
    fraction_cover = [0.8, 0.1, 0.1]
    size_landscape = 50
    ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros((length(pooling_vec) + 1) * n_keep, size(pseudo_param)[2] + n_metric) #for each site = 

    for n_ in 1:n_keep
        for k in 1:n_param
            summary_stat_table[((n_-1)*length(pooling_vec)+1):(n_*length(pooling_vec)), k] .= pseudo_param[n_, k]
        end
    end

    for i in 1:size(pseudo_param)[1]
        param[1:8] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_drylands_Schneider(time=200, param=copy(param), landscape=copy(ini_land), keep_landscape=false, n_snapshot=15, n_time_bw_snap=50)

        if d[200, 2] < 0.05 || d[200, 2] > 0.9 #in case there are missing issues
            for x in 1:length(pooling_vec)
                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = zeros(9)
            end
        else
            d, land = IBM_drylands_Schneider(time=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true, n_snapshot=15, n_time_bw_snap=50)

            for (x, pooling_intensity) in enumerate(pooling_vec)

                if pooling_intensity < 1
                    land_pooled = zeros(Int(1 / pooling_intensity) * size(land)[1], Int(1 / pooling_intensity) * size(land)[2], size(land)[3])
                elseif pooling_intensity == 1
                    land_pooled = copy(land)
                else
                    land_pooled = zeros(convert(Int64, floor(size(land)[1] / pooling_intensity)), convert(Int64, floor(size(land)[1] / pooling_intensity)), size(land)[3])
                end

                for third_sim in 1:size(land)[3] #we pool all saved landscapes 

                    if pooling_intensity < 1
                        land_pooled[:, :, third_sim] = inverse_pooling(land[:, :, third_sim], Int(1 / pooling_intensity)) #increasing resolution
                    else
                        land_pooled[:, :, third_sim] = pooling(land[:, :, third_sim], Int(pooling_intensity)) #decreasing resolution
                    end

                end

                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled)

            end
        end
    end

    CSV.write("./Data/Spatial_resolution/Simu/Schneider/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:40)


#endregion
#region, Step 4: Consequence of the spatial resolution on our ability to infer critical transitions


## Kefi model




using Distributed



addprocs(4, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end


@everywhere include("./Drylands_ABC_functions.jl")


@everywhere function Run_sim_model_EWS(model_id)


    if model_id == 1


        n_param = 7
        pooling_vec = [0.25 1 / 3 0.5 1]
        param = Get_classical_param()
        fraction_cover = [0.8, 0.1, 0.1]
        size_landscape = 50
        ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)

        b_seq = vcat(collect(0.39:0.01:0.6), collect(0.65:0.05:0.95)...)

        mat_sum_stat = zeros(length(b_seq) * length(pooling_vec), 11)

        for (y, b_val) in enumerate(b_seq)

            param[5] = b_val #changing conditions

            d, land = IBM_drylands(time=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true, n_snapshot=15, n_time_bw_snap=50)

            if sum(land[:, :, size(land)[3]]) != 0

                for (x, pooling_intensity) in enumerate(pooling_vec)

                    if pooling_intensity < 1
                        land_pooled = zeros(Int(1 / pooling_intensity) * size(land)[1], Int(1 / pooling_intensity) * size(land)[2], size(land)[3])
                    elseif pooling_intensity == 1
                        land_pooled = copy(land)
                    else
                        land_pooled = zeros(convert(Int64, floor(size(land)[1] / pooling_intensity)), convert(Int64, floor(size(land)[1] / pooling_intensity)), size(land)[3])
                    end

                    for third_sim in 1:size(land)[3] #we pool all saved landscapes 
                        if pooling_intensity < 1
                            land_pooled[:, :, third_sim] = inverse_pooling(land[:, :, third_sim], Int(1 / pooling_intensity)) #increasing resolution
                        else
                            land_pooled[:, :, third_sim] = pooling(land[:, :, third_sim], Int(pooling_intensity)) #decreasing resolution
                        end

                    end
                    mat_sum_stat[(y-1)*length(pooling_vec)+x, :] = vcat(Get_summary_stat(land_pooled), b_val, pooling_intensity...)
                end
            else
                mat_sum_stat[(y-1)*length(pooling_vec)+x, :] = vcat(zeros(9), b_val, pooling_intensity...)
            end
        end
        CSV.write("./Data/Spatial_resolution/Trends/Trends_Kefi.csv", Tables.table(mat_sum_stat), writeheader=false)


    elseif model_id == 2
        ##Schneider model 

        n_param = 8
        pooling_vec = [0.25 1 / 3 0.5 1]
        param = Get_classical_param_Schneider()
        fraction_cover = [0.8, 0.1, 0.1]
        size_landscape = 50
        ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)

        b_seq = collect(0.39:0.01:0.95)

        mat_sum_stat = zeros(length(b_seq) * length(pooling_vec), 11)

        for (y, b_val) in enumerate(b_seq)

            param[5] = b_val #changing conditions

            d, land = IBM_drylands_Schneider(time=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true, n_snapshot=15, n_time_bw_snap=50)

            if sum(land[:, :, size(land)[3]]) != 0

                for (x, pooling_intensity) in enumerate(pooling_vec)

                    if pooling_intensity < 1
                        land_pooled = zeros(Int(1 / pooling_intensity) * size(land)[1], Int(1 / pooling_intensity) * size(land)[2], size(land)[3])
                    elseif pooling_intensity == 1
                        land_pooled = copy(land)
                    else
                        land_pooled = zeros(convert(Int64, floor(size(land)[1] / pooling_intensity)), convert(Int64, floor(size(land)[1] / pooling_intensity)), size(land)[3])
                    end

                    for third_sim in 1:size(land)[3] #we pool all saved landscapes 
                        if pooling_intensity < 1
                            land_pooled[:, :, third_sim] = inverse_pooling(land[:, :, third_sim], Int(1 / pooling_intensity)) #increasing resolution
                        else
                            land_pooled[:, :, third_sim] = pooling(land[:, :, third_sim], Int(pooling_intensity)) #decreasing resolution
                        end

                    end
                    mat_sum_stat[(y-1)*length(pooling_vec)+x, :] = vcat(Get_summary_stat(land_pooled), b_val, pooling_intensity...)
                end
            else
                mat_sum_stat[(y-1)*length(pooling_vec)+x, :] = vcat(zeros(9), b_val, pooling_intensity...)
            end
        end
        CSV.write("./Data/Spatial_resolution/Trends/Trends_Schneider.csv", Tables.table(mat_sum_stat), writeheader=false)


    elseif model_id == 3




        ##Eby model 

        n_param = 2
        pooling_vec = [0.25 1 / 3 0.5 1]
        param = Get_classical_param_Eby()
        param[2] = 0.92 #same as in the paper
        fraction_cover = [0.8, 0.2]
        size_landscape = 50
        ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)

        p_seq = collect(0:0.01:0.95)

        mat_sum_stat = zeros(length(p_seq) * length(pooling_vec), 11)

        for (y, p_val) in enumerate(p_seq)

            param[1] = p_val #changing conditions

            d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true, n_snapshot=15, n_time_bw_snap=50)

            if sum(land[:, :, size(land)[3]]) != 0

                for (x, pooling_intensity) in enumerate(pooling_vec)

                    if pooling_intensity < 1
                        land_pooled = zeros(Int(1 / pooling_intensity) * size(land)[1], Int(1 / pooling_intensity) * size(land)[2], size(land)[3])
                    elseif pooling_intensity == 1
                        land_pooled = copy(land)
                    else
                        land_pooled = zeros(convert(Int64, floor(size(land)[1] / pooling_intensity)), convert(Int64, floor(size(land)[1] / pooling_intensity)), size(land)[3])
                    end

                    for third_sim in 1:size(land)[3] #we pool all saved landscapes 
                        if pooling_intensity < 1
                            land_pooled[:, :, third_sim] = inverse_pooling(land[:, :, third_sim], Int(1 / pooling_intensity)) #increasing resolution
                        else
                            land_pooled[:, :, third_sim] = pooling(land[:, :, third_sim], Int(pooling_intensity)) #decreasing resolution
                        end

                    end
                    mat_sum_stat[(y-1)*length(pooling_vec)+x, :] = vcat(Get_summary_stat(land_pooled), p_val, pooling_intensity...)
                end
            else
                mat_sum_stat[(y-1)*length(pooling_vec)+x, :] = vcat(zeros(9), p_val, pooling_intensity...)
            end
        end

        CSV.write("./Data/Spatial_resolution/Trends/Trends_Eby.csv", Tables.table(mat_sum_stat), writeheader=false)

    else






        ##Eby model with feedbacks

        n_param = 2
        pooling_vec = [0.25 1 / 3 0.5 1][1]
        param = Get_classical_param_Eby()
        param[2] = 0.4
        fraction_cover = [0.8, 0.2]
        size_landscape = 50
        ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)

        p_seq = collect(0:0.01:0.95)[1:4]

        mat_sum_stat = zeros(length(p_seq) * length(pooling_vec), 11)

        for (y, p_val) in enumerate(p_seq)

            param[1] = p_val #changing conditions

            d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land),
                keep_landscape=true, n_snapshot=15, n_time_bw_snap=50, intensity_feedback=6)

            if sum(land[:, :, size(land)[3]]) != 0

                for (x, pooling_intensity) in enumerate(pooling_vec)

                    if pooling_intensity < 1
                        land_pooled = zeros(Int(1 / pooling_intensity) * size(land)[1], Int(1 / pooling_intensity) * size(land)[2], size(land)[3])
                    elseif pooling_intensity == 1
                        land_pooled = copy(land)
                    else
                        land_pooled = zeros(convert(Int64, floor(size(land)[1] / pooling_intensity)), convert(Int64, floor(size(land)[1] / pooling_intensity)), size(land)[3])
                    end

                    for third_sim in 1:size(land)[3] #we pool all saved landscapes 
                        if pooling_intensity < 1
                            land_pooled[:, :, third_sim] = inverse_pooling(land[:, :, third_sim], Int(1 / pooling_intensity)) #increasing resolution
                        else
                            land_pooled[:, :, third_sim] = pooling(land[:, :, third_sim], Int(pooling_intensity)) #decreasing resolution
                        end

                    end
                    mat_sum_stat[(y-1)*length(pooling_vec)+x, :] = vcat(Get_summary_stat(land_pooled), p_val, pooling_intensity...)
                end
            else
                mat_sum_stat[(y-1)*length(pooling_vec)+x, :] = vcat(zeros(9), p_val, pooling_intensity...)
            end
        end

        CSV.write("./Data/Spatial_resolution/Trends/Trends_Eby_feedbacks.csv", Tables.table(mat_sum_stat), writeheader=false)

    end
end












#endregion
#region, Step 5: Posterior predictive check


post_p_NN = convert(Matrix, CSV.read("../Data/Step11_Inferrence/Post_pred_check/Posterior_pred_check_NN_1.csv",
    DataFrame, delim=";", header=false))
post_q_NN = convert(Matrix, CSV.read("../Data/Step11_Inferrence/Post_pred_check/Posterior_pred_check_NN_2.csv",
    DataFrame, delim=";", header=false))


post_p_rej = convert(Matrix, CSV.read("../Data/Step11_Inferrence/Post_pred_check/Posterior_pred_check_rej_1.csv",
    DataFrame, delim=";", header=false))
post_q_rej = convert(Matrix, CSV.read("../Data/Step11_Inferrence/Post_pred_check/Posterior_pred_check_rej_2.csv",
    DataFrame, delim=";", header=false))


name_pdf = "../Figures/ABC_scale/Posterior_check/Posterior_pred_check_post.pdf"
@rput name_pdf
R"pdf(name_pdf,width=8,height=8)"
@rput post_p_NN
@rput post_q_NN
@rput post_p_rej
@rput post_q_rej


for i in 1:size(post_q_rej)[2]

    R"par(mfrow=c(2,2))"

    #histogram of posterior
    @rput i
    R"hist(post_p_NN[,i])"
    R"hist(post_p_rej[,i])"
    R"hist(post_q_NN[,i])"
    R"hist(post_q_rej[,i])"
end
R"dev.off()"

sampled_sites = convert(Matrix, CSV.read("../Data/Step11_Inferrence/Post_pred_check/Sampled_sites.csv",
    DataFrame, delim=";", header=false))



name_pdf = "../Figures/ABC_scale/Posterior_check/Posterior_pred_check_landscapes.pdf"
@rput name_pdf
R"pdf(name_pdf,width=6,height=15)"

for i in 1:size(post_q_rej)[2]

    R"par(mfrow=c(5,2))"

    #observed vegetation
    land_emp = Get_empirical_site(sampled_sites[i])
    @rput land_emp
    color_land = ["white", "black"]
    @rput color_land
    R"image(land_emp,col=color_land)"
    R"image(land_emp,col=color_land)"

    #NeuralNet
    param = Get_classical_param_Eby()
    param[1] = median(post_p_NN[:, i])
    param[2] = median(post_q_NN[:, i])
    fraction_cover = [0.8, 0.2]
    size_landscape = 75
    ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    d1, land1 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
        keep_landscape=true, burning_phase=1500, intensity_feedback=6)

    #Loclinear
    param = Get_classical_param_Eby()
    param[1] = median(post_p_rej[:, i])
    param[2] = median(post_q_rej[:, i])
    fraction_cover = [0.8, 0.2]
    size_landscape = 75
    ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    d2, land2 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
        keep_landscape=true, burning_phase=1500, intensity_feedback=6)

    for k in 1:4:15
        landscape = copy(land1[:, :, k])
        @rput landscape
        color_land = ["white", "black"]
        @rput color_land
        R"image(landscape,col=color_land)"

        landscape = copy(land2[:, :, k])
        @rput landscape
        color_land = ["white", "black"]
        @rput color_land
        R"image(landscape,col=color_land)"
    end

    # R"par(mfrow=c(5,2))"

    # #observed vegetation
    # land_emp = Get_empirical_site(sampled_sites[i])
    # @rput land_emp
    # color_land = ["white", "black"]
    # @rput color_land
    # R"image(land_emp,col=color_land)"
    # R"image(land_emp,col=color_land)"

    # #NeuralNet
    # param = Get_classical_param_Eby()
    # param[1] = min(post_p_NN[:, i])
    # param[2] = min(post_q_NN[:, i])
    # fraction_cover = [0.8, 0.2]
    # size_landscape = 75
    # ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    # d1, land1 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
    #     keep_landscape=true, burning_phase=1500, intensity_feedback=6)

    # #Loclinear
    # param = Get_classical_param_Eby()
    # param[1] = min(post_p_rej[:, i])
    # param[2] = min(post_q_rej[:, i])
    # fraction_cover = [0.8, 0.2]
    # size_landscape = 75
    # ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    # d2, land2 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
    #     keep_landscape=true, burning_phase=1500, intensity_feedback=6)

    # for k in 1:4:15
    #     landscape = copy(land1[:, :, k])
    #     @rput landscape
    #     color_land = ["white", "black"]
    #     @rput color_land
    #     R"image(landscape,col=color_land)"

    #     landscape = copy(land2[:, :, k])
    #     @rput landscape
    #     color_land = ["white", "black"]
    #     @rput color_land
    #     R"image(landscape,col=color_land)"
    # end


    # R"par(mfrow=c(5,2))"

    # #observed vegetation
    # land_emp = Get_empirical_site(sampled_sites[i])
    # @rput land_emp
    # color_land = ["white", "black"]
    # @rput color_land
    # R"image(land_emp,col=color_land)"
    # R"image(land_emp,col=color_land)"

    # #NeuralNet
    # param = Get_classical_param_Eby()
    # param[1] = max(post_p_NN[:, i])
    # param[2] = max(post_q_NN[:, i])
    # fraction_cover = [0.8, 0.2]
    # size_landscape = 75
    # ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    # d1, land1 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
    #     keep_landscape=true, burning_phase=1500, intensity_feedback=6)

    # #Loclinear
    # param = Get_classical_param_Eby()
    # param[1] = max(post_p_rej[:, i])
    # param[2] = max(post_q_rej[:, i])
    # fraction_cover = [0.8, 0.2]
    # size_landscape = 75
    # ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    # d2, land2 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
    #     keep_landscape=true, burning_phase=1500, intensity_feedback=6)

    # for k in 1:4:15
    #     landscape = copy(land1[:, :, k])
    #     @rput landscape
    #     color_land = ["white", "black"]
    #     @rput color_land
    #     R"image(landscape,col=color_land)"

    #     landscape = copy(land2[:, :, k])
    #     @rput landscape
    #     color_land = ["white", "black"]
    #     @rput color_land
    #     R"image(landscape,col=color_land)"
    # end

    # print(i)
end

R"dev.off()"




#endregion
#region, Step 6: PPC and resampling from posterior parameters






using Distributed



addprocs(15, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end


@everywhere include("./Drylands_ABC_functions.jl")






@everywhere function Performing_pred_post_check(id)

    n_sample = 20

    post_p_NN = convert(Matrix, CSV.read("./Data/Posterior_pred_check/Posterior_pred_check_NN_p.csv",
        DataFrame, delim=";", header=false))
    post_q_NN = convert(Matrix, CSV.read("./Data/Posterior_pred_check/Posterior_pred_check_NN_q.csv",
        DataFrame, delim=";", header=false))


    post_p_rej = convert(Matrix, CSV.read("./Data/Posterior_pred_check/Posterior_pred_check_rej_p.csv",
        DataFrame, delim=";", header=false))
    post_q_rej = convert(Matrix, CSV.read("./Data/Posterior_pred_check/Posterior_pred_check_rej_q.csv",
        DataFrame, delim=";", header=false))

    post_scale = convert(Matrix, CSV.read("./Data/Posterior_pred_check/Posterior_pred_check_rej_scale.csv",
        DataFrame, delim=";", header=false))

    post_pred_NN = zeros(n_sample, 11)
    post_pred_rej = zeros(n_sample, 11)

    i = copy(id)

    scale_mat = post_scale[:, i]
    best_scale = collect(keys(countmap(scale_mat)))[findmax(collect(values(countmap(scale_mat))))[2]]

    for sample_id in 1:n_sample

        #NeuralNet
        param = Get_classical_param_Eby()
        sample_1 = sample(1:50)
        param[1] = post_p_NN[sample_1, i]
        param[2] = post_q_NN[sample_1, i]
        fraction_cover = [0.8, 0.2]
        size_landscape = 75
        ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
        d1, land1 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
            keep_landscape=true, burning_phase=1500, intensity_feedback=6)

        if any([length(findall(land1[:, :, l] .== 1)) for l in 1:size(land1)[3]] .== 0)
            post_pred_NN[sample_id, 3:11] = zeros(9)
        else
            if best_scale > 1
                land_pooled = zeros(Int(best_scale) * size(land1)[1], Int(best_scale) * size(land1)[2], size(land1)[3])
            else
                land_pooled = copy(land1) #no change
            end

            for third_sim in 1:size(land1)[3] #we pool all saved landscapes 
                if best_scale > 1
                    land_pooled[:, :, third_sim] = inverse_pooling(land1[:, :, third_sim], Int(best_scale)) #increasing resolution
                end
            end

            if best_scale > 1
                post_pred_NN[(sample_id), 3:11] = Get_summary_stat(land_pooled, xmin_fit=Int(best_scale)^2)
            else
                post_pred_NN[(sample_id), 3:11] = Get_summary_stat(land_pooled, xmin_fit=1)
            end
        end

        post_pred_NN[(sample_id), 1:2] .= [i, sample_id]

        #rejection
        param = Get_classical_param_Eby()
        sample_2 = sample(1:50)
        param[1] = post_p_rej[sample_2, i]
        param[2] = post_q_rej[sample_2, i]
        fraction_cover = [0.8, 0.2]
        size_landscape = 75
        ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
        d2, land2 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
            keep_landscape=true, burning_phase=1500, intensity_feedback=6)


        if any([length(findall(land2[:, :, l] .== 1)) for l in 1:size(land2)[3]] .== 0)
            post_pred_rej[sample_id, 3:11] = zeros(9)
        else
            if best_scale > 1
                land_pooled = zeros(Int(best_scale) * size(land2)[1], Int(best_scale) * size(land2)[2], size(land2)[3])
            else
                land_pooled = copy(land2) #no change
            end

            for third_sim in 1:size(land2)[3] #we pool all saved landscapes 
                if best_scale > 1
                    land_pooled[:, :, third_sim] = inverse_pooling(land2[:, :, third_sim], Int(best_scale)) #increasing resolution
                end
            end

            if best_scale > 1
                post_pred_rej[sample_id, 3:11] = Get_summary_stat(land_pooled, xmin_fit=Int(best_scale)^2)
            else
                post_pred_rej[sample_id, 3:11] = Get_summary_stat(land_pooled, xmin_fit=1)
            end
        end
        post_pred_rej[sample_id, 1:2] .= [i, sample_id]

    end


    CSV.write("./Data/Posterior_pred_check/Sim_per_site/Site_" * repr(i) * "_NN.csv", Tables.table(post_pred_NN), writeheader=false)
    CSV.write("./Data/Posterior_pred_check/Sim_per_site/Site_" * repr(i) * "_rej.csv", Tables.table(post_pred_rej), writeheader=false)
    print(i)

end

pmap(Performing_pred_post_check, 1:345)







@everywhere function Visual_post_pred_check(id)

    post_p_NN = convert(Matrix, CSV.read("./Data/Posterior_pred_check/Posterior_pred_check_NN_p.csv",
        DataFrame, delim=";", header=false))
    post_q_NN = convert(Matrix, CSV.read("./Data/Posterior_pred_check/Posterior_pred_check_NN_q.csv",
        DataFrame, delim=";", header=false))


    post_p_rej = convert(Matrix, CSV.read("./Data/Posterior_pred_check/Posterior_pred_check_rej_p.csv",
        DataFrame, delim=";", header=false))
    post_q_rej = convert(Matrix, CSV.read("./Data/Posterior_pred_check/Posterior_pred_check_rej_q.csv",
        DataFrame, delim=";", header=false))

    post_scale = convert(Matrix, CSV.read("./Data/Posterior_pred_check/Posterior_pred_check_rej_scale.csv",
        DataFrame, delim=";", header=false))


    name_pdf = "./Data/Posterior_pred_check/Visual/Posterior_pred_check_landscapes" * repr(id) * ".pdf"
    @rput name_pdf
    R"pdf(name_pdf,width=6,height=15)"


    for i in ((id-1)*23+1):(23*id)

        R"par(mfrow=c(5,2))"

        #observed vegetation
        land_emp = Get_empirical_site(i)
        @rput land_emp
        color_land = ["white", "black"]
        @rput color_land
        R"image(land_emp,col=color_land)"
        R"image(land_emp,col=color_land)"

        #NeuralNet
        param = Get_classical_param_Eby()
        param[1] = median(post_p_NN[:, i])
        param[2] = median(post_q_NN[:, i])
        fraction_cover = [0.8, 0.2]
        size_landscape = 75
        ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
        d1, land1 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
            keep_landscape=true, burning_phase=1500, intensity_feedback=6)

        #Loclinear
        param = Get_classical_param_Eby()
        param[1] = median(post_p_rej[:, i])
        param[2] = median(post_q_rej[:, i])
        fraction_cover = [0.8, 0.2]
        size_landscape = 75
        ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
        d2, land2 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
            keep_landscape=true, burning_phase=1500, intensity_feedback=6)

        for k in 1:4:15
            landscape = copy(land1[:, :, k])
            @rput landscape
            color_land = ["white", "black"]
            @rput color_land
            R"image(landscape,col=color_land)"

            landscape = copy(land2[:, :, k])
            @rput landscape
            color_land = ["white", "black"]
            @rput color_land
            R"image(landscape,col=color_land)"
        end

    end

    R"dev.off()"


end

pmap(Visual_post_pred_check, 1:15)

#endregion  
#region, Step 7: Computing the distance to a tipping point and the hysteresis size




using Distributed

addprocs(15, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end

@everywhere include("./Drylands_ABC_functions.jl")


@everywhere function Compute_cover_space(id)

    Param_site = Matrix{Float64}(CSV.read("../Data_new/Inferrence/Posterior_modes_each_sites.csv", DataFrame, header=1, delim=';')[id, :])

    Keeping_data = zeros(2000, 5)

    for traj in ["Degradation", "Restoration"]

        param = Get_classical_param_Eby()
        param[1] = Param_site[id, 3]
        param[2] = Param_site[id, 6]
        Keeping_data[:, 1:2] .= param

        if traj == "Degradation"
            fraction_cover = [0.8, 0.2]
        else
            fraction_cover = [0.1, 0.9]
        end

        size_landscape = 150
        ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)

        index = 1
        while param[1] != 0

            d1, land1 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
                keep_landscape=true, burning_phase=1500, intensity_feedback=6)

            mean_cover = mean([length(findall(land1[:, :, k] .== 1)) / (size(land1)[1] * size(land1)[2]) for k in 1:size(land1)[3]])

            Keeping_data[index, 3] = ifelse(any([length(findall(land1[:, :, k] .== 1)) == 0 for k in 1:size(land1)[3]]), 0, mean_cover)
            Keeping_data[index, 4] = ifelse(traj == "Degradation", 1, 2)

            index += 1
            param[1] -= 1 / 2000 #step size
            if param[1] < 0
                param[1] = 0
            end
        end

    end
    CSV.write("./Data/Inferrence/Dist_tipping_" * repr(id) * ".csv", Tables.table(Keeping_data), writeheader=false)
end


pmap(Compute_cover_space, 1:102)




#endregion
