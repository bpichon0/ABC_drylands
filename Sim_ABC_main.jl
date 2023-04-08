include("./Sim_ABC_functions.jl")

#region, Step 1: Generating the simulations

using Distributed

ncores = 25

addprocs(25, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end


@everywhere include("./Drylands_ABC_functions.jl")


@everywhere function Run_sim_model_EWS(N_sim)

    n_param = 2
    pooling_vec = [0.2 0.25 1 / 3 0.5 1]
    n_keep = 100
    pseudo_param = Matrix{Float64}(CSV.read("./Data/Parameter_sim.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :])
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

    CSV.write("./Data/Simulations/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:40)









#endregion
#region, Step 2: Computing the distance to a tipping point and the hysteresis size


# For each site, we compute the vegetation cover within the modes for p and q

using Distributed

ncores = 25
addprocs(ncores, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end

@everywhere include("./Drylands_ABC_functions.jl")


@everywhere function Distance_tipping(id)

    step_size = 1 / 200

    Param_site = Matrix{Float64}(CSV.read("./Data/Inferrence/Posterior_modes_each_sites.csv", DataFrame, header=1, delim=';')[id, :])
    p_25, p_50, p_75, q_25, q_50, q_75 = Param_site[id, 4], Param_site[id, 3], Param_site[id, 5], Param_site[id, 7], Param_site[id, 6], Param_site[id, 8]
    Keeping_data = zeros(((p_75 - p_25) / step_size) * ((q_75 - q_25) / step_size), 4)
    index = 1
    for traj in ["Degradation", "Restoration"]


        for p_ in p_25:step_size:p_75, q_ in q_25:step_size:p_75

            param = Get_classical_param_Eby()
            param[1] = p_
            param[2] = q_
            Keeping_data[index, 1:2] .= param

            if traj == "Degradation"
                fraction_cover = [0.8, 0.2]
            else
                fraction_cover = [0.2, 0.8]
            end

            size_landscape = 100
            ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)

            d1, land1 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
                keep_landscape=true, burning_phase=1500, intensity_feedback=6)

            mean_cover = mean([length(findall(land1[:, :, k] .== 1)) / (size(land1)[1] * size(land1)[2]) for k in 1:size(land1)[3]])

            Keeping_data[index, 3] = ifelse(any([length(findall(land1[:, :, k] .== 1)) == 0 for k in 1:size(land1)[3]]), 0, mean_cover)
            Keeping_data[index, 4] = ifelse(traj == "Degradation", 1, 2)

            index += 1

        end
    end
    CSV.write("./Data/Inferrence/Dist_tipping_" * repr(id) * ".csv", Tables.table(Keeping_data), writeheader=false)
end


pmap(Distance_tipping, 1:345)




#endregion
