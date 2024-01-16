begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end

include("../../Drylands_ABC_functions.jl")


param_list = readdlm("./Parameters_guichard.csv", ';', Float64)
print(size(param_list))

Keeping_data = zeros(size(param_list)[1], 18)
param = Param_Guichard()

for param_id in 1:size(param_list)[1]

    param["d"]=param_list[param_id,1]
    param["a0"]=param_list[param_id,2]
    param["a2"]=param_list[param_id,3]

    Keeping_data[param_id, 1] =param_list[param_id,1]
    Keeping_data[param_id, 2] =param_list[param_id,2]
    Keeping_data[param_id, 3] =param_list[param_id,3] 
    
    ini_land=Get_initial_lattice("guichard", ; size_landscape=100)

    d1, land1 = Run_model(model="guichard", param=copy(param), landscape=copy(ini_land), tmax=3000,
    keep_landscape=false, n_time_bw_snap=50, n_snapshot=25, intensity_feedback=1, burning_phase=2000)

    spatial_stats = Get_summary_stat( convert(Matrix{Int64},land1 .==2), ; log_=true, compute_psd=true, xmin_fit=1)

    Keeping_data[param_id, 4:14] = spatial_stats
    println(param_id)
		
end



CSV.write("./Stats_guichard.csv", Tables.table(Keeping_data), writeheader=false)
