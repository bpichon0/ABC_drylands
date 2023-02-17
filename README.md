# Code for running the inference of dryland resilience worldwide

Contact: Benoît Pichon, **benoit.pichon0@gmail.com**

All the code was made on R (v4.1.0) and Julia (1.7.3).

Here is the different steps to reproduce the figures:



## `Simulations`

The file `Sim_ABC_main.jl` allows to run the simulations made for the different models (Kefi et al., TPB, Schneider et al., 2016 TE, Eby et al., GEB 2017). It uses the fonctions in `Sim_ABC_function.jl`. 

## `Analyses`

To run all analyses made, you can run the R file `ABC_Drylands_main.R`.
It allows to perform the inference, as well as the different analyses made to understand the role of image resolution.
