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



#region, Step 3: Influence of the number of picture to average   


using Distributed



addprocs(2, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
    using LatinHypercubeSampling, JLD
end


@everywhere include("./Drylands_ABC_functions.jl")

@everywhere function Run_sim_model_EWS(N_sim, number_picture)

    pseudo_param = CSV.read("../Data/Pseudo_parameters.csv", DataFrame, header=1, delim=';')[((N_sim-1)*200+1):(N_sim*200), :]

    param = Get_classical_param()
    fraction_cover = [0.8, 0.1, 0.1]
    size_landscape = 50
    ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros(200, size(pseudo_param)[2] + n_metric)
    summary_stat_table[:, 1:7] .= pseudo_param

    for i in 1:size(pseudo_param)[1]
        param[1:7] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_drylands(time=2000, param=copy(param), landscape=copy(ini_land), keep_landscape=true, n_snapshot=number_picture)
        summary_stat_table[i, 8:size(summary_stat_table)[2]] = Get_summary_stat(land)

    end

    CSV.write("../Data/Step2_cross_validation/" * repr(number_picture) * "_pic" * "/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:250, [3, 5, 15, 25, 35])




#endregion



#region, Step 4: Influence of the number of parameters to infer   


using Distributed



addprocs(2, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
    using LatinHypercubeSampling, JLD
end

@everywhere include("./Drylands_ABC_functions.jl")

@everywhere function Run_sim_model_EWS(N_sim)

    pseudo_param = CSV.read("../Data/Pseudo_param_all_combinations.csv", DataFrame, header=1, delim=';')[((N_sim-1)*200+1):(N_sim*200), :]

    param = Get_classical_param()
    fraction_cover = [0.8, 0.1, 0.1]
    size_landscape = 50
    ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros(200, size(pseudo_param)[2] + n_metric)
    summary_stat_table[:, 1:7] .= pseudo_param

    for i in 1:size(pseudo_param)[1]
        param[1:7] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_drylands(time=2000, param=copy(param), landscape=copy(ini_land), keep_landscape=true, n_snapshot=25)
        summary_stat_table[i, 8:size(summary_stat_table)[2]] = Get_summary_stat(land)

    end

    CSV.write("../Data/Step2_cross_validation/Combination_parameters/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:2800)




#endregion



#region, Step XXXX: Correlation EWS along gradient mortality recruitment   

N_sim = 30
b_seq = collect(range(0.1, 1, length=N_sim))
m_seq = collect(range(0.005, 0.4, length=N_sim))


param = Get_classical_param()
fraction_cover = [0.8, 0.1, 0.1]
size_landscape = 50
ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)
n_metric = 9

summary_stat_table = zeros(N_sim^2, 7 + n_metric)

index = 1
for b_param in b_seq
    for m_param in m_seq

        param[5] = b_param
        param[4] = m_param
        summary_stat_table[index, 1:7] .= param[1:7]
        d, land = IBM_drylands(time=2000, param=copy(param), landscape=copy(ini_land), keep_landscape=true, n_snapshot=25)
        summary_stat_table[index, 8:size(summary_stat_table)[2]] = Get_summary_stat(land)
        index += 1
    end
end

CSV.write("../Data/Correlation_EWS_data.csv", Tables.table(summary_stat_table), writeheader=false)










#endregion



#region, Eby model simulation over the 2D space


addprocs(2, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
    using LatinHypercubeSampling, JLD
end

@everywhere include("./Drylands_ABC_functions.jl")

@everywhere function Run_sim_model_EWS(N_sim)

    pseudo_param = CSV.read("../Data/Eby_model/Pseudo_parameters_Eby.csv", DataFrame, header=1, delim=';')[((N_sim-1)*200+1):(N_sim*200), :]

    param = Get_classical_param()
    fraction_cover = [0.8, 0.1, 0.1]
    size_landscape = 50
    ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)
    n_metric = 9

    summary_stat_table = zeros(200, size(pseudo_param)[2] + n_metric)
    summary_stat_table[:, 1:2] .= pseudo_param

    @inbounds for i in 1:size(pseudo_param)[1]
        param[1:2] .= Vector{Float64}(pseudo_param[i, :])
        d, land = IBM_Eby_model(time_t=2000, param=copy(param), landscape=copy(ini_land), keep_landscape=true, n_snapshot=25)
        summary_stat_table[i, 3:size(summary_stat_table)[2]] = Get_summary_stat(land)

    end

    CSV.write("../Data/Eby_model/Simulation_ABC_number_" * repr(N_sim) * ".csv", Tables.table(summary_stat_table), writeheader=false)
end


pmap(Run_sim_model_EWS, 1:(100000/200))







#endregion

# example of Eby model landscapes
fraction_cover = [0.8, 0.2]
size_landscape = 100
ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)

name_plot = "../Figures/Eby_model_landscapes.pdf"
@rput name_plot

R"pdf(paste0(name_plot),width=6,height=6)"

for i in 1:50
    param = rand(2)
    d, land = IBM_Eby_model(time_t=2000, param=copy(param), landscape=copy(ini_land), keep_landscape=false, n_snapshot=25)
    @rput land
    R"image(as.matrix(land)>0)"
end
R"dev.off()"



#region Analysis of the empirical sites: getting the summarystats

list_land = readdir("../Data/Data_Biocom/landscapes/")

summary_stat = zeros(length(list_land), 9)


index = 1
for file_land in list_land
    print(index)
    summary_stat[index, :] = Get_summary_stat(Matrix{Int64}(CSV.read("../Data/Data_Biocom/landscapes/" * file_land, DataFrame, header=0)))
    index += 1
end

CSV.write("../Data/Data_Biocom/Summary_stats_data.csv", Tables.table(summary_stat), writeheader=false)



#endregion
