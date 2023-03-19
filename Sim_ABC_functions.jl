using StatsBase, RCall, Plots, StatsPlots, Random, LaTeXStrings,
    BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames,
    LatinHypercubeSampling, JLD, DelimitedFiles, RData

function Get_classical_param()

    r = 0.0001
    d = 0.2
    f = 0.9
    m = 0.1
    b = 1
    c = 0.3
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


function Plot_dynamics_Eby(d)

    plot(d[:, 1], seriescolor=:lightgreen, label="vegetation")
    ylims!((0.0, 1))

end


function Plot_landscape(landscape)
    landscape[findall((landscape .== -1))] .= 0
    colGRAD = cgrad([colorant"white", colorant"black"])
    heatmap(landscape, yflip=true, fill=true, c=colGRAD)
end

function Plot_landscape_Eby(landscape)
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


function Get_summary_stat(matrix_landscape, ; log_=true, compute_psd=true)
    R"library(spatialwarnings)"
    #The idea maybe is to compute as much as possible of sumamry statistics.

    matrix_landscape[findall(matrix_landscape .< 0)] .= 0

    if length(size(matrix_landscape)) == 3 # multiple landscapes

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
        cv_psd = vec(missings(Float64, size(matrix_landscape)[3], 1))
        median_psd = vec(missings(Float64, size(matrix_landscape)[3], 1))
        mean_psd = vec(missings(Float64, size(matrix_landscape)[3], 1))
        sd_psd = vec(missings(Float64, size(matrix_landscape)[3], 1))
        frac_max_psd = vec(missings(Float64, size(matrix_landscape)[3], 1))


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

            if compute_psd

                R"psd=spatialwarnings::patchdistr_sews(landscape)"
                R"PLR=spatialwarnings::raw_plrange(landscape)"
                R"if (nrow(psd$psd_type)==1){ 
                    alpha_exp=NA        
                } else {alpha_exp = psd$psd_type$plexpo[which(psd$psd_type$best==T)]}" #i.e., when there is no good fit, return NA
                R"psd_obs=psd$psd_obs"
                @rget PLR
                @rget alpha_exp
                @rget psd_obs

                PLR_vec[k] = ifelse(PLR === missing, 0, PLR)
                alpha_vec[k] = ifelse(alpha_exp === missing, 0, alpha_exp)

                cv_psd[k] = std(psd_obs) / mean(psd_obs)
                median_psd[k] = median(psd_obs)
                mean_psd[k] = mean(psd_obs)
                sd_psd[k] = std(psd_obs)
                frac_max_psd[k] = log(maximum(psd_obs) / size(landscape)[1]^2)


            else
                PLR_vec[k] = 0
                alpha_vec[k] = 0
                cv_psd[k] = 0
                median_psd[k] = 0
                mean_psd[k] = 0
                sd_psd[k] = 0
                frac_max_psd[k] = 0

            end

        end

        mean_nb_neigh = mean(nb_neigh)
        mean_clustering = mean(clustering)
        mean_spatial_skew = mean(spatial_skew)
        mean_spatial_var = mean(spatial_var)
        mean_spatial_corr = mean(spatial_corr)
        mean_spatial_spectral = mean(spatial_spectral)
        mean_PLR = mean(collect(skipmissing(PLR_vec)))
        mean_alpha = mean(collect(skipmissing(alpha_vec)))

        mean_cv_psd = mean(cv_psd)
        mean_median_psd = mean(median_psd)
        mean_mean_psd = mean(mean_psd)
        mean_sd_psd = mean(sd_psd)
        mean_frac_max_psd = mean(frac_max_psd)


    else  # only 1 landscape

        mean_cover = sum(matrix_landscape) / (size(matrix_landscape)[1]^2)

        # number of neighbors
        landscape = copy(matrix_landscape)
        landscape[findall(landscape .<= 0)] .= 0 #make it binary for spatial warnings package
        landscape[findall(landscape .> 0)] .= 1
        @rput landscape

        #vegetation clustering
        R"neighbors_mat = simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
        @rget neighbors_mat
        mean_nb_neigh = mean(neighbors_mat[findall(landscape .== 1)]) #mean number of plant neighbors
        mean_clustering = mean_nb_neigh / (length(findall(landscape .== 1)) / (size(landscape)[1] * size(landscape)[2]))

        #Spatial EWS using spatialwarnings package
        R"landscape=landscape > 0" #make it binary to fit the package input
        R"spatial_ews = generic_sews(landscape,4,moranI_coarse_grain = T)$value"
        @rget spatial_ews

        if sum(landscape) / (size(landscape)[1] * size(landscape)[2]) > 0.01
            mean_spatial_var = spatial_ews[1]
            mean_spatial_skew = spatial_ews[2]
            mean_spatial_corr = spatial_ews[3]
        else
            mean_spatial_var = 0
            mean_spatial_skew = 0
            mean_spatial_corr = 0
        end

        R"spectral_ratio = as.data.frame(spectral_sews(landscape,quiet=T))$value"
        @rget spectral_ratio
        mean_spatial_spectral = spectral_ratio

        if compute_psd

            R"psd=spatialwarnings::patchdistr_sews(landscape)"
            R"PLR=spatialwarnings::raw_plrange(landscape)"
            R"if (nrow(psd$psd_type)==1){ 
                alpha_exp=NA        
            } else {alpha_exp = psd$psd_type$plexpo[which(psd$psd_type$best==T)]}" #i.e., when there is no good fit, return NA
            @rget PLR
            @rget alpha_exp
            @rget psd_obs

            mean_PLR = ifelse(PLR === missing, 0, PLR)
            mean_alpha = ifelse(alpha_exp === missing, 0, alpha_exp)

            mean_cv_psd = std(psd_obs) / mean(psd_obs)
            mean_median_psd = median(psd_obs)
            mean_mean_psd = mean(psd_obs)
            mean_sd_psd = std(psd_obs)
            mean_frac_max_psd = log(maximum(psd_obs) / size(landscape)[1]^2)


        else
            mean_PLR = 0
            mean_alpha = 0
            mean_cv_psd = 0
            mean_median_psd = 0
            mean_mean_psd = 0
            mean_sd_psd = 0
            mean_frac_max_psd = 0
        end

    end

    if log_
        mean_spatial_spectral = log(mean_spatial_spectral)
        mean_clustering = log(mean_clustering)
    end


    return (vec([mean_cover, mean_nb_neigh, mean_clustering,
        mean_spatial_skew, mean_spatial_var, mean_spatial_corr,
        mean_spatial_spectral, mean_PLR, mean_alpha,
        mean_cv_psd, mean_median_psd, mean_mean_psd,
        mean_sd_psd, mean_frac_max_psd]))

end



function Get_classical_param_Eby(; p=0.8, q=0.8)
    return vec([
    p q])
end

function Get_initial_lattice_Eby(; frac=[0.9, 0.1], size_mat=50)

    ini_vec = sample([1, 0], Weights(frac), size_mat * size_mat)

    return reshape(ini_vec, size_mat, size_mat) #reshape by columns
end



function Get_neighbor(row, col, N)
    #Accounts for torus landscape (periodic boundaries)

    i, j = copy(row), copy(col)

    if i == 1
        top = N
    else
        top = i - 1
    end


    if i == N
        bottom = 1
    else
        bottom = i + 1
    end

    if j == N
        right = 1
    else
        right = j + 1
    end

    if j == 1
        left = N
    else
        left = j - 1
    end

    coordinate_n = zeros(4, 2)
    coordinate_n[1, :] = [i right]
    coordinate_n[2, :] = [i left]
    coordinate_n[3, :] = [top j]
    coordinate_n[4, :] = [bottom j]

    return convert(Array{Int64}, coordinate_n)
end


function Select_neighbors(neighbors)
    return convert(Array{Int64}, neighbors[shuffle(1:end), :][1, :])
end




function Select_neighbor_pair(coordinate_neighbors, Intensity_feedback)
    return convert(Array{Int64}, coordinate_neighbors[shuffle(1:end), :][1:Intensity_feedback, :])
end


function Get_coordinate_pair(row, col, row_n, col_n, N)
    i, j = copy(row), copy(col)
    i_n, j_n = copy(row_n), copy(col_n)

    #Boundaries condition for focal site (i,j) 

    if i == 1
        top = N
    else
        top = i - 1
    end

    if i == N
        bottom = 1
    else
        bottom = i + 1
    end

    if j == N
        right = 1
    else
        right = j + 1
    end

    if j == 1
        left = N
    else
        left = j - 1
    end

    #Boundaries condition for its neighbor (i_n,j_n)
    if i_n == 1
        topn = N
    else
        topn = i_n - 1
    end

    if i_n == N
        bottomn = 1
    else
        bottomn = i_n + 1
    end

    if j_n == N
        rightn = 1
    else
        rightn = j_n + 1
    end

    if j_n == 1
        leftn = N
    else
        leftn = j_n - 1
    end



    coordinate_n = zeros(8, 2)
    coordinate_n[1, :] = [i right]
    coordinate_n[2, :] = [i left]
    coordinate_n[3, :] = [top j]
    coordinate_n[4, :] = [bottom j]
    coordinate_n[5, :] = [i_n rightn]
    coordinate_n[6, :] = [i_n leftn]
    coordinate_n[7, :] = [topn j_n]
    coordinate_n[8, :] = [bottomn j_n]

    #and filter the ones corresponding to focal site and its selected neighbor
    coordinate_n = coordinate_n[Not(findall(coordinate_n[:, 1] .== i .&& coordinate_n[:, 2] .== j .||
                                            coordinate_n[:, 1] .== i_n .&& coordinate_n[:, 2] .== j_n)), :]

end


function IBM_Eby_model(; landscape, param, time_t, mortality_neigh=false, keep_landscape=false, n_snapshot=25, burning_phase=400, n_time_bw_snap=50, intensity_feedback=1, plot=false)

    if keep_landscape
        all_landscape_snap = zeros(size(landscape)[1], size(landscape)[1], n_snapshot)
        nsave = 1
    end

    p_param = param[1]
    q_param = param[2]
    #If we keep all snapshots we determine the minimum time for having n_snapshot after a burning_phase and with n_time_bw_snap time step between each

    if keep_landscape
        time_t = burning_phase + n_time_bw_snap * n_snapshot
    end
    d2 = zeros(time_t, 2) #Allocating
    d2[1, :] = vec([sum(landscape) / size(landscape)[1]^2, 1 - sum(landscape) / size(landscape)[1]^2])


    if mortality_neigh #context dependant mortality
        g0 = param[3]


        @inbounds for t in 1:time_t

            #for each time step, we perform N**2 (N the dimension of the landscape) iterations to ensure that all sites get change on average 1 time per time step  
            @inbounds for focal_i in eachindex(1:size(landscape)[1]), focal_j in eachindex(1:size(landscape)[1])


                if landscape[focal_i, focal_j] == 1 #if vegetation otherwise do nothing

                    neighbors = Get_neighbor(focal_i, focal_j, size(landscape)[1])
                    neigh = Select_neighbors(neighbors)

                    if landscape[neigh[1], neigh[2]] == 0 #if neighbor is unoccupied
                        rand1 = rand()
                        if rand1 <= p_param #then there is reproduction in a neighbor
                            landscape[neigh[1], neigh[2]] = 1
                        else
                            landscape[focal_i, focal_j] = 0
                        end

                    else #neighbor is occupied
                        rand2 = rand()
                        if rand1 <= q_param #facilitation from the neighbor
                            coord_neigh = Get_coordinate_pair(focal_i, focal_j, neigh[1], neigh[2], size(landscape)[1])
                            neighbors_pair = Select_neighbor_pair(coord_neigh, intensity_feedback) #that is changed to an occupied cell
                            for i in eachindex(neighbors_pair[:, 1])
                                landscape[neighbors_pair[i, 1], neighbors_pair[i, 2]] = 1
                            end
                        elseif rand2 <= 1 - p_param
                            landscape[focal_i, focal_j] = 0
                        end
                    end


                end #end loop for a focal cell

            end #loop on interactions

            @rput landscape
            R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
            @rget neigh_1

            rnum = reshape(rand(length(landscape)), Int64(sqrt(length(landscape))), Int64(sqrt(length(landscape))))
            # one random number between 0 and 1 for each cell

            landscape[findall((landscape .== 1) .& (rnum .<= g0 .* (1 .- neigh_1 / 4)))] .= 0

            d2[t, 1] = sum(landscape) / length(landscape)


            if keep_landscape && t > burning_phase && t % ((time_t - burning_phase) / n_snapshot) == 0 #for each time for which there is a snapshot, we save the landscape
                all_landscape_snap[:, :, nsave] = landscape
                nsave += 1
            end



        end

    else


        @inbounds for t in 1:time_t

            #for each time step, we perform N**2 (N the dimension of the landscape) iterations to ensure that all sites get change on average 1 time per time step  
            @inbounds for focal_i in eachindex(1:size(landscape)[1]), focal_j in eachindex(1:size(landscape)[1])

                if landscape[focal_i, focal_j] == 1 #if vegetation otherwise do nothing

                    neighbors = Get_neighbor(focal_i, focal_j, size(landscape)[1])
                    neigh = Select_neighbors(neighbors)


                    if landscape[neigh[1], neigh[2]] == 0 #if neighbor is unoccupied

                        if rand() <= p_param #then there is reproduction in a neighbor

                            landscape[neigh[1], neigh[2]] = 1

                        else #else focal individual dies

                            landscape[focal_i, focal_j] = 0

                        end

                    else #neighbor is occupied
                        if rand() <= q_param #facilitation from the neighbor

                            coord_neigh = Get_coordinate_pair(focal_i, focal_j, neigh[1], neigh[2], size(landscape)[1])
                            neighbor_pair = Select_neighbor_pair(coord_neigh, intensity_feedback) #that is changed to an occupied cell

                            for i in eachindex(neighbor_pair[:, 1])
                                landscape[neighbor_pair[i, 1], neighbor_pair[i, 2]] = 1
                            end

                        else
                            landscape[focal_i, focal_j] = 0
                        end
                    end
                end #end loop for a focal cell

            end #loop on interactions

            d2[t, 1] = sum(landscape) / length(landscape)


            if keep_landscape && t > burning_phase && t % ((time_t - burning_phase) / n_snapshot) == 0 #for each time for which there is a snapshot, we save the landscape
                all_landscape_snap[:, :, nsave] = landscape
                nsave += 1
            end



        end
    end

    d2[:, 2] = 1 .- d2[:, 1]

    if !keep_landscape
        all_landscape_snap = landscape
    end

    if plot
        if keep_landscape
            display(Plot_landscape(all_landscape_snap[:, :, size(all_landscape_snap)[3]], false))
        else
            display(Plot_landscape(all_landscape_snap, false))
        end
    end



    return d2, all_landscape_snap
end

function Get_empirical_site(id)
    d_biocom = CSV.read("../Data/Data_Biocom/biocom_data.csv", DataFrame, header=1, delim=';')

    return readdlm("../Data/Data_Biocom/landscapes/" * d_biocom[id, 2] * ".txt")

end




# schneider model

function Get_classical_param_Schneider()

    r = 0.0001
    d = 0.2
    f = 0.9
    m = 0.05
    b = 0.3
    c = 0.3
    delta = 0.1
    g0 = 0.2
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
        g0,
        z,
        tau_leap])
end

function IBM_drylands_Schneider(; landscape, param, time, keep_landscape=false, n_snapshot=25, burning_phase=1000, n_time_bw_snap=50)

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
    g0 = param[8]
    z = param[9]
    tau_leap = param[10]


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

        Rate_landscape[:, :, 2] .= @.(m + g0 * (1 - neigh_1 / z))
        Rate_landscape[:, :, 2] .= Rate_landscape[:, :, 2] .* (landscape .== 1)
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




function pooling(mat, submatrix_size) #decreasing resolution

    pooling_matrix = zeros(convert(Int64, floor(size(mat)[1] / submatrix_size)),
        convert(Int64, floor(size(mat)[1] / submatrix_size)))

    for i in 1:size(pooling_matrix)[1]
        for j in 1:size(pooling_matrix)[1]
            start_row = (i - 1) * submatrix_size + 1
            end_row = start_row + submatrix_size - 1
            start_col = (j - 1) * submatrix_size + 1
            end_col = start_col + submatrix_size - 1

            pooling_matrix[i, j] = ifelse(mean(mat[start_row:end_row, start_col:end_col]) > 0.5, 1, 0)
        end
    end
    return pooling_matrix
end




function inverse_pooling(mat, submatrix_size) #increasing resolution
    n = size(mat, 1)
    m = size(mat, 2)
    new_n = submatrix_size * n
    new_m = submatrix_size * m
    new_mat = zeros(new_n, new_m)

    for i in 1:n
        for j in 1:m
            start_row = (i - 1) * submatrix_size + 1
            end_row = start_row + submatrix_size - 1
            start_col = (j - 1) * submatrix_size + 1
            end_col = start_col + submatrix_size - 1

            new_mat[start_row:end_row, start_col:end_col] .= mat[i, j]
        end
    end

    return new_mat
end



function expand_grid(; iters...)
    var_names = collect(keys(iters))
    var_itr = [1:length(x) for x in iters.data]
    var_ix = vcat([collect(x)' for x in Iterators.product(var_itr...)]...)
    out = DataFrame()
    for i = 1:length(var_names)
        out[:, var_names[i]] = collect(iters[i])[var_ix[:, i]]
    end
    return out
end
