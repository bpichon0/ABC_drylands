using StatsBase, RCall, Plots, StatsPlots, Random, LaTeXStrings,
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


function Plot_landscape(landscape, with_fertile)
    if with_fertile
        landscape[findall((landscape .== 0))] .= 1
    else
        landscape[findall((landscape .== 0))] .= -1
    end
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


function Get_summary_stat(matrix_landscape)
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

        R"psd=spatialwarnings::patchdistr_sews(landscape)"
        R"max_patchsize=max(psd$psd_obs)"
        R"PLR=spatialwarnings::raw_plrange(landscape)"
        R"if (nrow(psd$psd_type)==1){ 
            alpha_exp=NA        
        } else {alpha_exp = psd$psd_type$plexpo[which(psd$psd_type$best==T)]}" #i.e., when there is no good fit, return NA
        @rget PLR
        @rget alpha_exp

        mean_PLR = ifelse(PLR === missing, 0, PLR)
        mean_alpha = ifelse(alpha_exp === missing, 0, alpha_exp)

    end



    return (vec([mean_cover, mean_nb_neigh, mean_clustering,
        mean_spatial_skew, mean_spatial_var, mean_spatial_corr,
        mean_spatial_spectral, mean_PLR, mean_alpha]))

end



function Get_classical_param_Eby(; p=0.8, q=0.8)
    return vec([
    p q])
end

function Get_initial_lattice_Eby(; frac=[0.9, 0.1], size_mat=50)

    ini_vec = sample([1, 0], Weights(frac), size_mat * size_mat)

    return reshape(ini_vec, size_mat, size_mat) #reshape by columns
end



function select_neighbor(row, col, N, N_neigh)
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

    test = rand()

    col_neigh = j
    row_neigh = i

    if N_neigh == 4

        if test <= 0.25
            row_neigh = top
        elseif test <= 0.5
            row_neigh = bottom
        elseif test <= 0.75
            col_neigh = left
        else
            col_neigh = right
        end

    else
        if test <= 1 / 8 #top one
            row_neigh = top
        elseif test <= 2 / 8 #bottom one
            row_neigh = bottom
        elseif test <= 3 / 8 #left one 
            col_neigh = left
        elseif test <= 4 / 8 #right one
            col_neigh = right
        elseif test <= 5 / 8 #upper left
            col_neigh = left
            row_neigh = top
        elseif test <= 6 / 8 #bottom left
            col_neigh = left
            row_neigh = bottom
        elseif test <= 7 / 8 #upper right
            col_neigh = right
            row_neigh = top
        else                 # bottom left
            col_neigh = right
            row_neigh = bottom
        end
    end

    return row_neigh, col_neigh

end







function select_neighbor_pair2(coordinate_neighbors, Intensity_feedback)
    return convert(Array{Int64}, coordinate_neighbors[shuffle(1:end), :][1:Intensity_feedback, :])
end


function Get_coordinate(row, col, row_n, col_n, N, N_neigh)
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


    if N_neigh == 4

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

    else

        coordinate_n = zeros(16, 2)
        coordinate_n[1, :] = [i right]
        coordinate_n[2, :] = [i left]
        coordinate_n[3, :] = [top j]
        coordinate_n[4, :] = [bottom j]
        coordinate_n[5, :] = [bottom right]
        coordinate_n[6, :] = [bottom left]
        coordinate_n[7, :] = [top right]
        coordinate_n[8, :] = [top left]
        coordinate_n[9, :] = [i_n rightn]
        coordinate_n[10, :] = [i_n leftn]
        coordinate_n[11, :] = [topn j_n]
        coordinate_n[12, :] = [bottomn j_n]
        coordinate_n[13, :] = [bottomn rightn]
        coordinate_n[14, :] = [bottomn leftn]
        coordinate_n[15, :] = [topn rightn]
        coordinate_n[16, :] = [topn leftn]

        #and filter the ones corresponding to focal site and its selected neighbor
        coordinate_n = coordinate_n[Not(findall(coordinate_n[:, 1] .== i .&& coordinate_n[:, 2] .== j .||
                                                coordinate_n[:, 1] .== i_n .&& coordinate_n[:, 2] .== j_n)), :]
        final_coord = unique(coordinate_n, dims=1) #selecting the uniques sites

    end
end






function IBM_Eby_model(; landscape, param, time_t, keep_landscape=false, n_snapshot=25, burning_phase=400, n_time_bw_snap=50, Type_neighbors=4, intensity_feedback=1)

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

    @inbounds for t in 1:time_t

        #for each time step, we perform N**2 (N the dimension of the landscape) iterations to ensure that all sites get change on average 1 time per time step  
        @inbounds for focal_i in eachindex(1:size(landscape)[1]), focal_j in eachindex(1:size(landscape)[1])

            if landscape[focal_i, focal_j] == 1 #if vegetation otherwise do nothing

                neigh = select_neighbor(focal_i, focal_j, size(landscape)[1], Type_neighbors)

                if landscape[neigh[1], neigh[2]] == 0 #if neighbor is unoccupied
                    if rand() <= p_param #then there is reproduction in a neighbor
                        landscape[neigh[1], neigh[2]] = 1
                    else #else focal individual dies
                        landscape[focal_i, focal_j] = 0
                    end

                else #neighbor is occupied
                    if rand() <= q_param #facilitation from the neighbor
                        coord_neigh = Get_coordinate(focal_i, focal_j, neigh[1], neigh[2], size(landscape)[1], Type_neighbors)
                        neighbors = select_neighbor_pair2(coord_neigh, intensity_feedback) #that is changed to an occupied cell
                        for i in eachindex(neighbors[:, 1])
                            landscape[neighbors[i, 1], neighbors[i, 2]] = 1
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

    d2[:, 2] = 1 .- d2[:, 1]

    if !keep_landscape
        all_landscape_snap = landscape
    end



    return d2, all_landscape_snap
end
