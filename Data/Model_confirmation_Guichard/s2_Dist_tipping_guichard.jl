
using Distributed

addprocs(60, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end

@everywhere include("../../Drylands_ABC_functions.jl")




@everywhere function Distance_tipping_guichard(id)

    step_size = 1 / 300

    posteriors = readdlm("./Parameters_guichard.csv", ';', Float64)[id, :]

    Keeping_data = zeros(Int((1 / step_size)), 4)

    p_to_desert = posteriors[1]:step_size:1
    param = Param_Guichard()

    for pcrit_id in 1:length(p_to_desert)

	if rand(Distributions.Uniform(0, 1)) < .5
          println("")
          GC.gc()
          ccall(:malloc_trim, Cvoid, (Cint,), 0)
          GC.safepoint()
        end

        posteriors[1] = p_to_desert[pcrit_id]

        Keeping_data[pcrit_id, 1:3] = posteriors[1:3]

        param["d"]=posteriors[1]
        param["a0"]=posteriors[2]
        param["a2"]=posteriors[3]
	
        ini_land=Get_initial_lattice("guichard", ; size_landscape=100)

        d1, land1 = Run_model(model="guichard", param=copy(param), landscape=copy(ini_land), tmax=3000,
           keep_landscape=true, n_time_bw_snap=50, n_snapshot=25, intensity_feedback=1, burning_phase=1500)

        mean_cover = mean([length(findall(land1[:, :, k] .== 2)) / (size(land1)[1] * size(land1)[2]) for k in 1:size(land1)[3]])

        Keeping_data[pcrit_id, 4] = ifelse(any([length(findall(land1[:, :, k] .== 2)) == 0 for k in 1:size(land1)[3]]), 0, mean_cover)
	if Keeping_data[pcrit_id, 4]<0.002
	  break
	end

        println(pcrit_id)

    end

    CSV.write("./Dist_guichard/Dist_tipping_" * repr(id) * ".csv", Tables.table(Keeping_data), writeheader=false)
end

pmap(Distance_tipping_guichard, 1:60)
