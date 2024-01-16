begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end

include("../../Drylands_ABC_functions.jl")


param_list = readdlm("./Parameters_kefi.csv", ';', Float64)
print(size(param_list))

Keeping_data = zeros(size(param_list)[1], 18)
param = Get_classical_param()

for param_id in 1:size(param_list)[1]

    param[3] = param_list[param_id, 5]
    param[5] = param_list[param_id, 6]
    param[7] = param_list[param_id, 7]
    param[6] = param_list[param_id, 4]
    param[1] = param_list[param_id, 1]
    param[2] = param_list[param_id, 2]
    param[4] = param_list[param_id, 3]
    

    fraction_cover = [0.8, 0.1, 0.1]

    Keeping_data[param_id, 1:7] = param[1:7]

    size_landscape = 100
    ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)

    d1, land1 = IBM_drylands(time=50, param=copy(param), landscape=copy(ini_land),
        keep_landscape=true, burning_phase=2000)

    spatial_stats = Get_summary_stat(land1, ; log_=true, compute_psd=true, xmin_fit=1)

    Keeping_data[param_id, 8:18] = spatial_stats
    println(param_id)
		
end



CSV.write("./Stats_kefi.csv", Tables.table(Keeping_data), writeheader=false)
