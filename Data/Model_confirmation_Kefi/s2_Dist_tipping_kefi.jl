
using Distributed

addprocs(30, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end

@everywhere include("../../Drylands_ABC_functions.jl")




@everywhere function Distance_tipping_kefi(id)

    step_size = 1 / 200

    posteriors = readdlm("./Parameters_kefi.csv", ';', Float64)[id, :]

    Keeping_data = zeros(Int((1 / step_size)), 4)

    p_to_desert = reverse(push!(collect(0:step_size:posteriors[6]),posteriors[6]))
    param = Get_classical_param()

    for pcrit_id in 1:length(p_to_desert)

	if rand(Distributions.Uniform(0, 1)) < .5
          println("")
          GC.gc()
          ccall(:malloc_trim, Cvoid, (Cint,), 0)
          GC.safepoint()
        end


        posteriors[6] = p_to_desert[pcrit_id]

        Keeping_data[pcrit_id, 1:3] = posteriors[5:7]
        fraction_cover = [0.8, 0.1, 0.1]

        param[3] = posteriors[5]
        param[5] = posteriors[6]
        param[7] = posteriors[7]
    	param[6] = posteriors[4]
    	param[1] = posteriors[1]
    	param[2] = posteriors[2]
    	param[4] = posteriors[3]



        size_landscape = 100
        ini_land = Get_initial_lattice(frac=fraction_cover, size_mat=size_landscape)

        d1, land1 = IBM_drylands(time=50, param=copy(param), landscape=copy(ini_land),
            keep_landscape=true, burning_phase=2000)

        mean_cover = mean([length(findall(land1[:, :, k] .== 1)) / (size(land1)[1] * size(land1)[2]) for k in 1:size(land1)[3]])

        Keeping_data[pcrit_id, 4] = ifelse(any([length(findall(land1[:, :, k] .== 1)) == 0 for k in 1:size(land1)[3]]), 0, mean_cover)
	if Keeping_data[pcrit_id, 4]<0.002
	  break
	end

        println(pcrit_id)

    end


    CSV.write("./Dist_kefi/Dist_tipping_" * repr(id) * ".csv", Tables.table(Keeping_data), writeheader=false)
end

pmap(Distance_tipping_kefi, 1:60)
