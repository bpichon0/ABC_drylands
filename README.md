# Code for running the inference of dryland resilience worldwide

Contact: Benoît Pichon, **benoit.pichon0@gmail.com**

Prior to analyses, see the preprint: Benoît Pichon, Sophie Donnet, Isabelle Gounand and Sonia Kefi: Learning from vegetation patterns to infer the resilience of drylands, *bioRxiv*.

This repository contains the code used to perform the analyses for both main text and supplementary informations.
All the code was made on R (*v4.1.0*) and in julia (*v1.7.3*).

<p align="center">
    <img src="https://github.com/bpichon0/ABC_drylands/blob/master/Example/Framework_ABC.jpg" width="800">
</p>

## `Installing R & Julia dependancies`


To install the pacakges needed for the analyses and create the folder architecture used to save the data sets, please load the file `Dryland_shift_function.R` using: 

```R
source(ABC_drylands_function.R)
```
In addition, to install julia dependancies go into your working directory, press "]" and enter "activate .". Once activating, all dependancies can be loaded using "] instantiate".


## `Replicating the analyses`

In order to replicate the analyses, both simulations and observed vegetation landscapes are necessary.

- For the data, we provide the data on empirical site in the file `data_sites.csv` within the folder `Data`. This *csv* file contains the summary statistics on the empirical sites as well as the environmental data used in the paper. 

- The simulations are made using julia for computational speed. To generate the simulations, you can run the Step 1 of `Sim_ABC_main.jl` file. We nevertheless, provide a subset of simulations, with $(p, q) \sim \mathcal{U}[0,1]$, in the Data folder to avoid too long computational time (`simulations.csv`).

Next, all the analyses exept from the simulations for prediction are made with the R file `ABC_drylands_main.R`. 
The file is organized in explicit sections (Alt+O to see them). All section starts with the explanation of the code.
 

## `Replicating the Figures`

The file `Make_figs.R` contains the different steps to replicate all analysis in the paper. It is organized into different chunks of code, each doing a separate figure. **To generate a specific figure, run the section id indicaded at the begining of the section of Make_figs.R.**



