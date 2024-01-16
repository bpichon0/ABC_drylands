include("./Sim_ABC_functions.jl")

#region, Step 1: Generating the simulations

# Code to generate the all simulations from uniform priors 
# for p and q and different spatial resolution
# This can take a while depending on the number of cores used (~2 days with 50 cores CPU)


using Distributed


ncores = 50
addprocs(ncores, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end


@everywhere include("./Drylands_ABC_functions.jl")


@everywhere function Run_sim_model_EWS(N_sim)

    n_param = 2
    n_sim_kept = 1
    pooling_vec = [1 1 / 2 1 / 3 1 / 4 1 / 5]
    n_keep = 200
    pseudo_param = CSV.read("./Data/Parameters.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :]

    param = Get_classical_param_Eby()
    fraction_cover = [0.8, 0.2]
    size_landscape = 100
    ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 11

    summary_stat_table = zeros((length(pooling_vec) + 1) * n_keep * n_sim_kept, size(pseudo_param)[2] + n_metric) #for each site = 

    for n_ in 1:n_keep
        for k in 1:n_param
            summary_stat_table[((n_-1)*length(pooling_vec)*n_sim_kept+1):(n_*length(pooling_vec)*n_sim_kept), k] .= pseudo_param[n_, k]
        end
    end

    for i in 1:size(pseudo_param)[1]
        param[1] = pseudo_param[i, 1]
        param[2] = pseudo_param[i, 2]
        d, land = IBM_Eby_model(time_t=100, param=copy(param), landscape=copy(ini_land),
            keep_landscape=false, n_snapshot=n_sim_kept, n_time_bw_snap=50, intensity_feedback=1)

        if d[100, 1] < 0.04 || d[100, 1] > 0.9 #in case there are missing issues
            summary_stat_table[(1+(i-1)*length(pooling_vec)*n_sim_kept):(i)*length(pooling_vec)*n_sim_kept, (n_param+1):size(summary_stat_table)[2]] .= zeros(length(pooling_vec) * n_sim_kept, n_metric)
        else
            d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land),
                keep_landscape=true, n_snapshot=n_sim_kept, n_time_bw_snap=50, burning_phase=1500, intensity_feedback=1)

            if any([length(findall(land[:, :, l] .== 1)) for l in 1:size(land)[3]] .== 0)
                summary_stat_table[(1+(i-1)*length(pooling_vec)*n_sim_kept):(i)*length(pooling_vec)*n_sim_kept, (n_param+1):size(summary_stat_table)[2]] .= zeros(length(pooling_vec) * n_sim_kept, n_metric)
            else
                for (x, pooling_intensity) in enumerate(pooling_vec) #each spatial resolution
                    for sub_mat in 1:size(land)[3] #for each of the landscapes at asymptotic state

                        if pooling_intensity < 1
                            land_pooled = inverse_pooling(land[:, :, sub_mat], Int(1 / pooling_intensity)) #increasing resolution
                        else #no changes
                            land_pooled = copy(land[:, :, sub_mat])
                        end

                        if pooling_intensity < 1
                            summary_stat_table[((i-1)*length(pooling_vec)*n_sim_kept+(sub_mat-1)*length(pooling_vec)+x), (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled, xmin_fit=Int(1 / pooling_intensity)^2, compute_psd=false)
                        else #no changes
                            summary_stat_table[((i-1)*length(pooling_vec)*n_sim_kept+(sub_mat-1)*length(pooling_vec)+x), (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled, xmin_fit=1)
                        end
                    end
                end
            end
        end
    end
    CSV.write("./Data/Simulations/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:1200) #1200*200 pairs of parameters








#endregion
#region, Step 2: Computing the distance to a tipping point


using Distributed

addprocs(33, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end

@everywhere include("./Drylands_ABC_functions.jl")


@everywhere function Distance_tipping(id)

    step_size = 1 / 200
    n_sample=100
    posteriors=readdlm("./Data/posterior_param.csv", ';')[2:end,2:end]
    post_p = posteriors[:, id]
    post_q = posteriors[:, id+345]
    Keeping_data = zeros(Int((1/step_size))*n_sample, 3)
    index=1
    param=zeros(2)
    println(id)
    for sample_id in 1:n_sample

        if rand(Distributions.Uniform(0, 1)) < .5
          println("")
          GC.gc()
          ccall(:malloc_trim, Cvoid, (Cint,), 0)
          GC.safepoint()
        end
        sample_row=rand(1:length(post_p))
        param[1] = post_p[sample_row]
        param[2] = post_q[sample_row]

        p_to_desert = push!(collect(0:step_size:param[1]),param[1])

        for pcrit_id in 1:length(p_to_desert)

            param[1]=p_to_desert[pcrit_id]

            Keeping_data[index, 1:2] .= param
            fraction_cover = [0.8, 0.2]
            
            size_landscape = 80
            ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)

            d1, land1 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
            keep_landscape=true, burning_phase=2000, intensity_feedback=6)

            mean_cover = mean([length(findall(land1[:, :, k] .== 1)) / (size(land1)[1] * size(land1)[2]) for k in 1:size(land1)[3]])

            Keeping_data[index, 3] = ifelse(any([length(findall(land1[:, :, k] .== 1)) == 0 for k in 1:size(land1)[3]]), 0, mean_cover)
            
            index += 1
            end
    end

    CSV.write("./Data/Prediction/Dist_tipping_" * repr(id) * ".csv", Tables.table(Keeping_data), writeheader=false)
end

sites=readdlm("./Data/Keeping_sites.csv",';',Int64)[:,1]
pmap(Distance_tipping, sites)





#endregion
#region, Step 3: Sensitivity to system site (# of pixels)



addprocs(10, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end


@everywhere include("./Drylands_ABC_functions.jl")

@everywhere function Run_sim_model_EWS(N_sim)

    n_param = 2
    n_sim_kept = 15
    size_vec = [75 125 175 225]
    n_keep = 100
    pseudo_param = CSV.read("./Data/Parameters.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :]

    param = Get_classical_param_Eby()
    fraction_cover = [0.8, 0.2]
    n_metric = 11

    summary_stat_table = zeros((length(size_vec) + 1) * n_keep * n_sim_kept, size(pseudo_param)[2] + n_metric) #for each site = 

    for n_ in 1:n_keep
        for k in 1:n_param
            summary_stat_table[((n_-1)*length(size_vec)*n_sim_kept+1):(n_*length(size_vec)*n_sim_kept), k] .= pseudo_param[n_, k]
        end
    end

    for i in 1:size(pseudo_param)[1]
        param[1] = pseudo_param[i, 1]
        param[2] = pseudo_param[i, 2]
        d, land = IBM_Eby_model(time_t=100, param=copy(param), landscape=copy(Get_initial_lattice_Eby(frac=fraction_cover, size_mat=75)),
            keep_landscape=false, n_snapshot=n_sim_kept, n_time_bw_snap=25, intensity_feedback=6)

        if d[100, 1] < 0.04 || d[100, 1] > 0.9 #in case there are missing issues

            summary_stat_table[(1+(i-1)*length(size_vec)*n_sim_kept):(i)*length(size_vec)*n_sim_kept, (n_param+1):size(summary_stat_table)[2]] .= zeros(length(size_vec) * n_sim_kept, n_metric)

        else

            for (x, size_landscape) in enumerate(size_vec) #each spatial resolution

                d, land = IBM_Eby_model(time_t=100, param=copy(param), landscape=copy(Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)),
                    keep_landscape=true, n_snapshot=n_sim_kept, n_time_bw_snap=25, intensity_feedback=6)

                if any([length(findall(land[:, :, l] .== 1)) for l in 1:size(land)[3]] .== 0)

                    summary_stat_table[(1+(i-1)*length(size_vec)*n_sim_kept):(i)*length(size_vec)*n_sim_kept, (n_param+1):size(summary_stat_table)[2]] .= zeros(length(size_vec) * n_sim_kept, n_metric)

                else

                    for sub_mat in 1:size(land)[3] #for each of the landscapes at asymptotic state

                        land_pooled = copy(land[:, :, sub_mat])

                        summary_stat_table[((i-1)*length(size_vec)*n_sim_kept+(sub_mat-1)*length(size_vec)+x), (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled, xmin_fit=1)

                    end
                end
            end

        end
    end

    CSV.write("./Data/System_size/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:10) #randomly selected parameter sets



#endregion
