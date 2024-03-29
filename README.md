# Code for running the inference of dryland resilience worldwide

> [!NOTE]
> Code to replicate the results of the following preprint: Benoît Pichon, Sophie Donnet, Isabelle Gounand and Sonia Kefi: Estimating distances to tipping points from dryland ecosystem images, *bioRxiv* 2024.
> Contact: Benoît Pichon, **benoit.pichon0@gmail.com**


This repository contains the code used to perform the analyses for both main text and supplementary information.
All the code was made on R (*v4.1.0*) and in julia (*v1.7.3*).

<p align="center">
    <img src="https://github.com/bpichon0/ABC_drylands/blob/master/Example/Framework_ABC.jpg" width="800">
</p>

## `Installing R & Julia dependencies`

> [!IMPORTANT]  
> To install the packages needed for the analyses and create the folder architecture used to save the data sets, please load the file `ABC_drylands_function.R` using: 

```R
source(ABC_drylands_function.R)
```
> [!IMPORTANT]  
> In addition, to install Julia dependencies go into your working directory, press "]" and enter "activate .". Once activating, all dependencies can be loaded using "] instantiate".

## `Replicating the analyses`

People interested in directly replicating the figures can go to "**Replicating the figures**" since we provide all the generated data.

### Step 1: Simulations 

In order to replicate the analyses, we give the spatial statistics computed on the observed and simulated landscapes, as detailed below:

- For the data, we provide the data on empirical site in the file `data_sites.csv` within the folder `Data`. This *csv* file contains the summary statistics on the empirical sites as well as the environmental data used in the paper. 

- The simulations are made using Julia for computational speed. To generate the simulations, you can run the Step 1 of the `Sim_ABC_main.jl` file. Then, the outputs of these simulations are post-processed using the first section of `ABC_drylands_main.R` (Step 0: Merging simulations). We nevertheless provide a subset of simulations (`simulations.csv` file), with $(p, q) \sim \mathcal{U}[0,1]$, in the `Data` folder to avoid too long computational time.
Note that if you want to generate outputs from IBM simulations, we provide also generic functions to generate different IBM models and plot them (<a href="https://github.com/bpichon0/IBM_models">link functions on GitHub<a>)  

### Step 2: Inference part

To infer the parameters and the distance to the tipping point of each observed landscape. Here are the steps: 

- First, run the Step 1.1 and 1.2 in `ABC_drylands_main.R` to infer the parameters and remove sites with bimodal distributions. ***To see the steps press Alt+O***

- Then, Run the Step 2 from `Sim_ABC_main.jl` to estimate the distance to the tipping point.

- Last, run the Step 1.3 in `ABC_drylands_main.R` to merge the outputs from the latter step and save the csv.

### Step 2bis: Methodology around the inference approach

To replicate the analyses made around the inference approach (sensitivity, choice summary statistics etc..), please run the Step 2 in `ABC_drylands_main.R` (chunks 2.1 to 2.5). Also, run the last Step (Step 3) from the `Sim_ABC_main.jl` file to analyse the sensitivity of spatial metrics to system size.
These methodological analyses are nevertheless not needed for the rest of the analyses.

### Step 3: Statistical analyses 

To replicate statistical analyses, please run Step 3 to 6 in the `ABC_drylands_main.R` file (respectively fitting mixed-effects linear models, Running the SEM, comparing predictive power of vegetation cover and/or spatial structure of vegetation, and analysing climatic data).

### Step 4: Comparison to other models 

Last, to replicate the comparison of our approach with other models, we provide a complete framework through a *bash (.sh)* file in `./Data/Model_confirmation_Guichard` and `./Data/Model_confirmation_Kefi`.
The framework draws the parameters, performs the simulations with the two models, infers the parameters, and estimates the distance to the tipping point by our approach and in the two models (mussel-bed and dryland vegetation). 
Then, to postprocess these simulations, you can run Step 7 of the `ABC_drylands_main.R` file.

## `Replicating the figures`

The file `Make_figs.R` plot the figures. As we provide all the data to replicate each figure, everything can be runed without replicating the different analyses. The file is organized in different chunks of code, that can be seen by pressing *Alt+O*.
We give the main figures in the Figure folder but all supplementary figures can be replicated using `Make_figs.R`.

## `Reading the paper`

To have more information about these figures, read the following preprint: Benoît Pichon, Sophie Donnet, Isabelle Gounand and Sonia Kefi: Estimating distances to tipping points from dryland ecosystem images, *bioRxiv* 2024.

<p align="center">
    <img src="https://github.com/bpichon0/ABC_drylands/blob/master/Example/Validation_simulations.jpg" width="800">
</p>
