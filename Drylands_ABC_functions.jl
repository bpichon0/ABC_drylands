using StatsBase, RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings,
    BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames,
    LatinHypercubeSampling, JLD

function Get_classical_param()

    r = 0.05
    d = 0.1
    f = 0.9
    m = 0.1
    b = 1
    c = 0.1
    delta = 0.1
    z = 4
    tau_leap = 0.5


    return vec([
        r,
        d,
        f,
        m,
        b,
        c,
        delta,
        z,
        tau_leap])
end


function Get_initial_lattice(; frac=[0.8, 0.1, 0.1], size_mat=50)

    ini_vec = sample([1, 0, -1], Weights(frac), size_mat * size_mat)

    return reshape(ini_vec, size_mat, size_mat) #reshape by columns
end




function Plot_dynamics(d)

    plot(d[:, 1], d[:, 2], seriescolor=:lightgreen, label="vegetation")
    plot!(d[:, 1], d[:, 3], seriescolor=:orange, label="fertile")
    plot!(d[:, 1], d[:, 4], seriescolor=:grey, label="degraded")
    ylims!((0.0, 1))


end

function Plot_landscape(landscape)

    landscape[findall((landscape .== -1))] .= 0
    colGRAD = cgrad([colorant"#696969", colorant"#ACD87B"])
    heatmap(landscape, yflip=true, fill=true, c=colGRAD)
end


function IBM_drylands(; landscape, param, time, keep_landscape=false, n_snapshot=25, burning_phase=1000, n_time_bw_snap=50)

    if keep_landscape
        all_landscape_snap = zeros(size(landscape)[1], size(landscape)[1], n_snapshot)
        nsave = 1
    end
    r = param[1]
    d = param[2]
    f = param[3]
    m = param[4]
    b = param[5]
    c = param[6]
    delta = param[7]
    z = param[8]
    tau_leap = param[9]


    rules_change = transpose([0 1 -1 0; 1 0 0 -1])

    #Global densities
    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol

    nb_cell = size(landscape)[1]

    #Allocating 
    Rate_landscape = zeros(nb_cell, nb_cell, 4)

    #If we keep all snapshots we determine the minimum time for having n_snapshot after a burning_phase and with n_time_bw_snap time step between each

    if keep_landscape
        time = burning_phase + n_time_bw_snap * n_snapshot
    end
    d2 = zeros(time, 4) #Allocating



    for t in 1:time


        @rput landscape
        R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
        @rget neigh_1

        Rate_landscape[:, :, 1] .= @. (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(b - c * rho_1)
        Rate_landscape[:, :, 1] .= Rate_landscape[:, :, 1] .* (landscape .== 0)

        Rate_landscape[:, :, 2] .= m .* (landscape .== 1)
        Rate_landscape[:, :, 3] .= @.(r + f * neigh_1 / z) .* (landscape .== -1)
        Rate_landscape[:, :, 4] .= d .* (landscape .== 0)

        Rate_landscape[findall(Rate_landscape .< 0)] .= 0 #to avoid problems with propensity

        #calculate propensity

        propensity = [sum(Rate_landscape[:, :, k]) for k in 1:size(Rate_landscape)[3]]

        nb_events = map(x -> rand(Poisson(x)), propensity * tau_leap)

        for event in 1:length(nb_events) #for each type of events
            patches = findall(landscape .== rules_change[event, 1])

            if nb_events[event] != 0 && length(patches) > nb_events[event]
                landscape[wsample(patches, Rate_landscape[patches, event], nb_events[event])] .= rules_change[event, 2]
            end
        end

        rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction vegetation
        rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
        rho_d = 1 - rho_1 - rho_f # fraction degraded

        @views d2[t, :] = [t rho_1 rho_f rho_d]
        Rate_landscape = zeros(nb_cell, nb_cell, 4)

        #keeping the landscapes to average the summary statistics
        if keep_landscape && t > burning_phase && t % ((time - burning_phase) / n_snapshot) == 0
            all_landscape_snap[:, :, nsave] = landscape
            nsave += 1
        end

    end

    if !keep_landscape
        all_landscape_snap = landscape
    end

    return d2, all_landscape_snap
end


function Get_summary_stat(matrix_landscape)
    R"library(spatialwarnings)"
    #The idea maybe is to compute as much as possible of sumamry statistics.

    matrix_landscape[findall(matrix_landscape .< 0)] .= 0

    # vegetation cover
    mean_cover = mean([length(findall(matrix_landscape[:, :, k] .== 1)) / (size(matrix_landscape)[1] * size(matrix_landscape)[2]) for k in 1:size(matrix_landscape)[3]])

    # number of neighbors
    clustering = vec(missings(Float64, size(matrix_landscape)[3], 1))
    nb_neigh = vec(missings(Float64, size(matrix_landscape)[3], 1))
    spatial_skew = vec(missings(Float64, size(matrix_landscape)[3], 1))
    spatial_var = vec(missings(Float64, size(matrix_landscape)[3], 1))
    spatial_corr = vec(missings(Float64, size(matrix_landscape)[3], 1))
    spatial_spectral = vec(missings(Float64, size(matrix_landscape)[3], 1))
    PLR_vec = vec(missings(Float64, size(matrix_landscape)[3], 1))
    alpha_vec = vec(missings(Float64, size(matrix_landscape)[3], 1))

    for k in 1:size(matrix_landscape)[3]

        landscape = matrix_landscape[:, :, k]
        landscape[findall(landscape .<= 0)] .= 0 #make it binary for spatial warnings package
        landscape[findall(landscape .> 0)] .= 1
        @rput landscape

        #vegetation clustering
        R"neighbors_mat = simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
        @rget neighbors_mat
        nb_neigh[k] = mean(neighbors_mat[findall(landscape .== 1)]) #mean number of plant neighbors
        clustering[k] = nb_neigh[k] / (length(findall(landscape .== 1)) / (size(landscape)[1] * size(landscape)[2]))

        #Spatial EWS using spatialwarnings package
        R"landscape=landscape > 0" #make it binary to fit the package input
        R"spatial_ews = generic_sews(landscape,4,moranI_coarse_grain = T)$value"
        @rget spatial_ews

        if sum(landscape) / (size(landscape)[1] * size(landscape)[2]) > 0.01
            spatial_var[k] = spatial_ews[1]
            spatial_skew[k] = spatial_ews[2]
            spatial_corr[k] = spatial_ews[3]
        else
            spatial_var[k] = 0
            spatial_skew[k] = 0
            spatial_corr[k] = 0
        end

        R"spectral_ratio = as.data.frame(spectral_sews(landscape,quiet=T))$value"
        @rget spectral_ratio
        spatial_spectral[k] = spectral_ratio

        R"psd=spatialwarnings::patchdistr_sews(landscape)"
        R"max_patchsize=max(psd$psd_obs)"
        R"PLR=spatialwarnings::raw_plrange(landscape)"
        R"if (nrow(psd$psd_type)==1){ 
            alpha_exp=NA        
        } else {alpha_exp = psd$psd_type$plexpo[which(psd$psd_type$best==T)]}" #i.e., when there is no good fit, return NA
        @rget PLR
        @rget alpha_exp

        PLR_vec[k] = PLR
        alpha_vec[k] = alpha_exp

    end

    mean_nb_neigh = mean(nb_neigh)
    mean_clustering = mean(clustering)
    mean_spatial_skew = mean(spatial_skew)
    mean_spatial_var = mean(spatial_var)
    mean_spatial_corr = mean(spatial_corr)
    mean_spatial_spectral = mean(spatial_spectral)
    mean_PLR = mean(collect(skipmissing(PLR_vec)))
    mean_alpha = mean(collect(skipmissing(alpha_vec)))

    return (vec([mean_cover, mean_nb_neigh, mean_clustering,
        mean_spatial_skew, mean_spatial_var, mean_spatial_corr,
        mean_spatial_spectral, mean_PLR, mean_alpha]))

end
