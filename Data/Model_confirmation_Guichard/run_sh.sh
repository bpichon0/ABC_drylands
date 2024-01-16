#! /bin/bash

Rscript s0_param_guichard.R

julia s1_stats_guichard_sim.jl

julia s2_Dist_tipping_guichard.jl

Rscript s3_Inferrence_Guichard.R

julia s4_prediction.jl
