include("./Drylands_ABC_functions.jl")



param = Get_classical_param_Eby(p=0.5, q=0.8)
ini_land = Get_initial_lattice_Eby()
d, land = IBM_Eby_model(landscape=copy(ini_land), param=copy(param), time_t=15000)
Plot_dynamics_Eby(d)

#region, Step 1: Simulating of the pseudo-data (from the model)


addprocs(2, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
    using LatinHypercubeSampling, JLD
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

    CSV.write("./Data/Simulations/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end

pmap(Run_sim_model_EWS, 1:1000)







#endregion

#region 2) Getting the summary statistics of the empirical sites

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
#region Spatial resolution


using Distributed



addprocs(1, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random,  LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end


@everywhere include("./Drylands_ABC_functions.jl")


@everywhere function Run_sim_model_EWS(N_sim)

    n_param=2
    pooling_vec=[0.25 1/3 .5 1]
    n_keep=100
    pseudo_param = Matrix{Float64}(CSV.read("./Data/Spatial_resolution/Parameter_Eby.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :])
    param = Get_classical_param_Eby()
    fraction_cover = [0.8, 0.2]
    size_landscape = 50
    ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros((length(pooling_vec)+1)*n_keep, size(pseudo_param)[2] + n_metric) #for each site = 

    for n_ in 1:n_keep
        for k in 1:n_param
            summary_stat_table[((n_-1)*length(pooling_vec)+1):(n_*length(pooling_vec)), k] .= pseudo_param[n_,k]
        end
    end

    for i in 1:size(pseudo_param)[1]
        param[1:2] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land), keep_landscape=false,n_snapshot = 15,n_time_bw_snap = 50)

        if d[200,1]<0.05 || d[200,1]>0.9 #in case there are missing issues
            for x in 1:length(pooling_vec)
                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = zeros(9)
            end
        else    
            d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true,n_snapshot = 15,n_time_bw_snap = 50)

            for (x,pooling_intensity) in enumerate(pooling_vec)
                
                if pooling_intensity < 1
                    land_pooled=zeros(Int(1/pooling_intensity)*size(land)[1],Int(1/pooling_intensity)*size(land)[2],size(land)[3])
                elseif pooling_intensity==1
                    land_pooled=copy(land)
                else
                    land_pooled=zeros(convert(Int64,floor(size(land)[1]/pooling_intensity)),convert(Int64,floor(size(land)[1]/pooling_intensity)),size(land)[3])
                end
                                
                for third_sim in 1:size(land)[3] #we pool all saved landscapes 
                    if pooling_intensity < 1
                        land_pooled[:,:,third_sim] = inverse_pooling(land[:,:,third_sim],Int(1/pooling_intensity)) #increasing resolution
                    else 
                        land_pooled[:,:,third_sim] = pooling(land[:,:,third_sim],Int(pooling_intensity)) #decreasing resolution
                    end

                end
                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled)

            end
        end
    end

    CSV.write("./Data/Spatial_resolution/Simu/Eby/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:250)






@everywhere function Run_sim_model_EWS(N_sim)

    n_param=2
    pooling_vec=[0.25 1/3 .5 1]
    n_keep=100
    pseudo_param = CSV.read("./Data/Spatial_resolution/Parameter_Eby_feedback.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :]

    param = Get_classical_param_Eby()
    fraction_cover = [0.8, 0.2]
    size_landscape = 50
    ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros((length(pooling_vec)+1)*n_keep, size(pseudo_param)[2] + n_metric) #for each site = 

    for n_ in 1:n_keep
        for k in 1:n_param
            summary_stat_table[((n_-1)*length(pooling_vec)+1):(n_*length(pooling_vec)), k] .= pseudo_param[n_,k]
        end
    end


    for i in 1:size(pseudo_param)[1]

        param[1:2] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land), keep_landscape=false,
                    n_snapshot = 15,n_time_bw_snap = 50,intensity_feedback = 6)

        if d[200,1]<0.05 || d[200,1]>0.9 #in case there are missing issues
            for x in 1:length(pooling_vec)
                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = zeros(9)
            end
        else    
            d, land = IBM_Eby_model(time_t=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true,
            n_snapshot = 15,n_time_bw_snap = 50,intensity_feedback = 6)

            for (x,pooling_intensity) in enumerate(pooling_vec)

                if pooling_intensity < 1
                    land_pooled=zeros(Int(1/pooling_intensity)*size(land)[1],Int(1/pooling_intensity)*size(land)[2],size(land)[3])
                elseif pooling_intensity==1
                    land_pooled=copy(land)
                else
                    land_pooled=zeros(convert(Int64,floor(size(land)[1]/pooling_intensity)),convert(Int64,floor(size(land)[1]/pooling_intensity)),size(land)[3])
                end

                for third_sim in 1:size(land)[3] #we pool all saved landscapes 

                    if pooling_intensity < 1
                        land_pooled[:,:,third_sim] = inverse_pooling(land[:,:,third_sim],Int(1/pooling_intensity)) #increasing resolution
                    else 
                        land_pooled[:,:,third_sim] = pooling(land[:,:,third_sim],Int(pooling_intensity)) #decreasing resolution
                    end

                end

                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled)

            end
        end
    end

    CSV.write("./Data/Spatial_resolution/Simu/Eby_feedback/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:250)








@everywhere function Run_sim_model_EWS(N_sim)

    n_param=7
    pooling_vec=[0.25 1/3 .5 1]
    n_keep=100
    pseudo_param = CSV.read("./Data/Spatial_resolution/Parameter_Kefi.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :]

    param = Get_classical_param()
    fraction_cover = [0.8, 0.1,0.1]
    size_landscape = 50
    ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros((length(pooling_vec)+1)*n_keep, size(pseudo_param)[2] + n_metric) #for each site = 

    for n_ in 1:n_keep
        for k in 1:n_param
            summary_stat_table[((n_-1)*length(pooling_vec)+1):(n_*length(pooling_vec)), k] .= pseudo_param[n_,k]
        end
    end

    for i in 1:size(pseudo_param)[1]
        param[1:7] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_drylands(time=200, param=copy(param), landscape=copy(ini_land), keep_landscape=false,n_snapshot = 15,n_time_bw_snap = 50)


        if d[200,2]<0.05 || d[200,2]>0.9 #in case there are missing issues
            for x in 1:length(pooling_vec)
                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = zeros(9)
            end
        else    
            
            d, land = IBM_drylands(time=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true,n_snapshot = 15,n_time_bw_snap = 50)

            for (x,pooling_intensity) in enumerate(pooling_vec)

                if pooling_intensity < 1
                    land_pooled=zeros(Int(1/pooling_intensity)*size(land)[1],Int(1/pooling_intensity)*size(land)[2],size(land)[3])
                elseif pooling_intensity==1
                    land_pooled=copy(land)
                else
                    land_pooled=zeros(convert(Int64,floor(size(land)[1]/pooling_intensity)),convert(Int64,floor(size(land)[1]/pooling_intensity)),size(land)[3])
                end

                
                for third_sim in 1:size(land)[3] #we pool all saved landscapes 

                    if pooling_intensity < 1
                        land_pooled[:,:,third_sim] = inverse_pooling(land[:,:,third_sim],Int(1/pooling_intensity)) #increasing resolution
                    else 
                        land_pooled[:,:,third_sim] = pooling(land[:,:,third_sim],Int(pooling_intensity)) #decreasing resolution
                    end

                end

                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled)

            end
        end
    end

    CSV.write("./Data/Spatial_resolution/Simu/Kefi/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)




end


pmap(Run_sim_model_EWS, 1:250)






@everywhere function Run_sim_model_EWS(N_sim)

    n_param=8
    pooling_vec=[0.25 1/3 .5 1]
    n_keep=100
    pseudo_param = CSV.read("./Data/Spatial_resolution/Parameter_Schneider.csv", DataFrame, header=1, delim=';')[((N_sim-1)*n_keep+1):(N_sim*n_keep), :]

    param = Get_classical_param_Schneider()
    fraction_cover = [0.8, 0.1,0.1]
    size_landscape = 50
    ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros((length(pooling_vec)+1)*n_keep, size(pseudo_param)[2] + n_metric) #for each site = 

    for n_ in 1:n_keep
        for k in 1:n_param
            summary_stat_table[((n_-1)*length(pooling_vec)+1):(n_*length(pooling_vec)), k] .= pseudo_param[n_,k]
        end
    end
    
    for i in 1:size(pseudo_param)[1]
        param[1:8] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_drylands_Schneider(time=200, param=copy(param), landscape=copy(ini_land), keep_landscape=false,n_snapshot = 15,n_time_bw_snap = 50)

        if d[200,2]<0.05 || d[200,2]>0.9 #in case there are missing issues
            for x in 1:length(pooling_vec)
                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = zeros(9)
            end
        else    
            d, land = IBM_drylands_Schneider(time=200, param=copy(param), landscape=copy(ini_land), keep_landscape=true,n_snapshot = 15,n_time_bw_snap = 50)

            for (x,pooling_intensity) in enumerate(pooling_vec)

                if pooling_intensity < 1
                    land_pooled=zeros(Int(1/pooling_intensity)*size(land)[1],Int(1/pooling_intensity)*size(land)[2],size(land)[3])
                elseif pooling_intensity==1
                    land_pooled=copy(land)
                else
                    land_pooled=zeros(convert(Int64,floor(size(land)[1]/pooling_intensity)),convert(Int64,floor(size(land)[1]/pooling_intensity)),size(land)[3])
                end
                
                for third_sim in 1:size(land)[3] #we pool all saved landscapes 

                    if pooling_intensity < 1
                        land_pooled[:,:,third_sim] = inverse_pooling(land[:,:,third_sim],Int(1/pooling_intensity)) #increasing resolution
                    else 
                        land_pooled[:,:,third_sim] = pooling(land[:,:,third_sim],Int(pooling_intensity)) #decreasing resolution
                    end

                end

                summary_stat_table[(i-1)*length(pooling_vec)+x, (n_param+1):size(summary_stat_table)[2]] = Get_summary_stat(land_pooled)

            end
        end
    end

    CSV.write("./Data/Spatial_resolution/Simu/Schneider/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:250)
