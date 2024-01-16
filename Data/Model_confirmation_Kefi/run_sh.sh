#! /bin/bash

Rscript s0_param_kefi.R

julia s1_stats_Kefi_sim.jl

julia s2_Dist_tipping_kefi.jl

Rscript s3_Inferrence_Kefi.R

julia s4_prediction.jl
