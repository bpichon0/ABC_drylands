
using Distributed

addprocs(19, exeflags="--project=$(Base.active_project())")

@everywhere begin
    using StatsBase, RCall, Random, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames
end

@everywhere include("../../Drylands_ABC_functions.jl")


@everywhere function Distance_tipping(id)

    step_size = 1 / 200
    n_sample = 100
    posteriors = readdlm("./param_rej.csv", ';', Float64)[:, 2:end]
    post_p = posteriors[:, id]

    post_q = posteriors[:, id+60]

    Keeping_data = zeros(Int((1 / step_size)) * n_sample, 3)
    index = 1
    param = zeros(2)
    for sample_id in 1:n_sample
	
	if rand(Distributions.Uniform(0, 1)) < 1
          println("")
          GC.gc()
          ccall(:malloc_trim, Cvoid, (Cint,), 0)
          GC.safepoint()
        end
	println(sample_id)
        
	sample_row = rand(1:length(post_p))
        param[1] = post_p[sample_row]
        param[2] = post_q[sample_row]

        p_to_desert = reverse(push!(collect(0:step_size:param[1]),param[1]))

        for pcrit_id in 1:length(p_to_desert)

            param[1] = p_to_desert[pcrit_id]

            Keeping_data[index, 1:2] .= param
            fraction_cover = [0.8, 0.2]
	    
            size_landscape = 100
            ini_land = Get_initial_lattice_Eby(frac=fraction_cover, size_mat=size_landscape)

            d1, land1 = IBM_Eby_model(time_t=50, param=copy(param), landscape=copy(ini_land),
                keep_landscape=true, burning_phase=1500, intensity_feedback=1,n_snapshot=10)

            mean_cover = mean([length(findall(land1[:, :, k] .== 1)) / (size(land1)[1] * size(land1)[2]) for k in 1:size(land1)[3]])

            Keeping_data[index, 3] = ifelse(any([length(findall(land1[:, :, k] .== 1)) == 0 for k in 1:size(land1)[3]]), 0, mean_cover)

            index += 1
        end
    end

    CSV.write("./Dist_Eby/Dist_tipping_Eby" * repr(id) * ".csv", Tables.table(Keeping_data), writeheader=false)
end

print(sites)
pmap(Distance_tipping, sites)

