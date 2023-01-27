# Code for running the inference of dryland resilience worldwide

Contact: Benoît Pichon, **benoit.pichon0@gmail.com**

This folder contains all code necessary to replicate the analysis in the main text and in supplementary. 
The file **Dryland_ABC_main.jl** contains the code to generate the simulations of both Kefi et al TPB 2007 model and Eby et al., GEB 2017. It uses the functions in **Dryland_ABC_functions.jl**. 

All inference analyses are made with R using the **ABC_drylands_main.R** script. 
All the steps necessary to *(i)* perform ABC and optimize the ABC method with the simulations (from both Kefi and Eby models), and in the future *(ii)* compare the simulations and empirical data, *(iii)* perform ABC on empirical data, *(iv)* measure site specific resilience for all sites across the globe. 

