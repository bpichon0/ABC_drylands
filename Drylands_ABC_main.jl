include("./Drylands_ABC_functions.jl")

#region, Testing region
param = Get_classical_param()
fraction_cover = [0.8, 0.1, 0.1]
size_landscape = 50
ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)
@time d2, land2 = IBM_drylands(time=1500, param=copy(param),
    landscape=copy(ini_land), keep_landscape=true)

Get_summary_stat(land2)
Plot_dynamics(d2)

#endregion



#region, Step 1: Sensitivity analysis on the model 

# First, naive exploration with every parameter fixed exept one
param_sensitivity = CSV.read("../Data/Step1_sensitivity/sensitivity_1D_param.csv",
    DataFrame, header=1, delim=';')

n_metric = 9
all_data = zeros(size(param_sensitivity)[1], size(param_sensitivity)[2] + n_metric)
all_data[:, 1:7] .= param_sensitivity


param = Get_classical_param()
fraction_cover = [0.8, 0.1, 0.1]
size_landscape = 50
ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)


for i in 1:size(param_sensitivity)[1]
    param[1:7] .= Vector{Float64}(param_sensitivity[i, :])
    d, land = IBM_drylands(time=2000, param=copy(param), landscape=copy(ini_land),
        keep_landscape=true, burning_phase=500, n_time_bw_snap=50, n_snapshot=30)
    all_data[i, 8:size(all_data)[2]] = Get_summary_stat(land)
    #save("../data/Step1_sensitivity/Save_landscape_1D/Sim" * repr(i) * ".jld", "data", land)
    print(repr(i) * "  ")
end

CSV.write("../Data/Step1_sensitivity/All_data_1D.csv", Tables.table(all_data), writeheader=false)



# Second, analysing the pair parameters in the space (p1, p2) 


param_sensitivity = CSV.read("../Data/Step1_sensitivity/sensitivity_2D_param.csv",
    DataFrame, header=1, delim=';')

n_metric = 9
all_data = zeros(size(param_sensitivity)[1], size(param_sensitivity)[2] + n_metric)
all_data[:, 1:7] .= param_sensitivity

param = Get_classical_param()
fraction_cover = [0.8, 0.1, 0.1]
size_landscape = 50
ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)


for i in 1:size(param_sensitivity)[1]

    param[1:7] .= Vector{Float64}(param_sensitivity[i, :])
    d, land = IBM_drylands(time=2000, param=copy(param), landscape=copy(ini_land),
        keep_landscape=true, burning_phase=500, n_time_bw_snap=50, n_snapshot=15)
    all_data[i, 8:size(all_data)[2]] = Get_summary_stat(land)
    #save("../data/Step1_sensitivity/Save_landscape_2D/Sim" * repr(i) * ".jld", "data", land)
    print(repr(i) * "  ")

end

CSV.write("../Data/Step1_sensitivity/All_data_2D.csv", Tables.table(all_data), writeheader=false)


#endregion



#region, Step 2: Simulating pseudo datasets for inference



using Distributed



addprocs(2, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
    using LatinHypercubeSampling, JLD
end


@everywhere include("./Drylands_ABC_functions.jl")




@everywhere function Run_sim_model_EWS(N_sim)

    pseudo_param = CSV.read("../Data/Pseudo_parameters.csv", DataFrame, header=1, delim=';')[((N_sim-1)*100+1):(N_sim*100), :]

    param = Get_classical_param()
    fraction_cover = [0.8, 0.1, 0.1]
    size_landscape = 50
    ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros(100, size(pseudo_param)[2] + n_metric)
    summary_stat_table[:, 1:7] .= pseudo_param

    for i in 1:size(pseudo_param)[1]
        param[1:7] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_drylands(time=2000, param=copy(param), landscape=copy(ini_land), keep_landscape=true)
        summary_stat_table[i, 8:size(summary_stat_table)[2]] = Get_summary_stat(land)

    end

    CSV.write("../Data/Step2_cross_validation/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:500)






#endregion
