""" 
This file provides the space to run your GEM model.
Before you start, please read README.txt
"""

# run only once the first time
# include("install_pkgs.jl")

# after closing and re-opening, to have Julia running in the same environment that was originally created, 
# go to the GitHub\ folder and then type "activate ." in the REPO

# load all required packages 
include("functions/Packages.jl")

# load model definition  
include("bd_IGM_model_definition.jl")

# load model configuration
include("bd_IGM_model_config.jl")

# load all functions 
include("functions/GEM_Functions.jl")

# run the GEM simulation
run_sim = GEM_sim(
                  N0, # initial state
                  model_par_vect, # model parameters in a vector
                  gem_const_vect,
                  design_choices, # evolution decisions
                  mappings, # parameter-to-state mapping
                  sim_params,# simulation parameters
                  sim_output, # output containers
                  verbose=false # show time on console
                  ) #

# output: 
# Tuple{Array{Float64, 4}, Array{Float64, 5}, Array{Float64, 5}}

# dataframe for population time series
pop_dat = run_sim.pop_df
CSV.write("pop_time_series.csv", pop_dat)
# 2 dataframes: mean and variance
trait_dat = run_sim.trait_df

# accessing the two trait dataframes
# trait mean dataframe:
trait_dat.median 
CSV.write("trait_mean_time_series.csv", trait_dat.median)

# trait variance dataframeß
trait_dat.var
CSV.write("trait_var_time_series.csv", trait_dat.var)
# =====================================================
# Pop_Plot(pop data, stateID)
Pop_Plot(pop_data = pop_dat, stateID = 2, add_mean=true)

# Trait_Plot(mean, var, stateID, "trait name")
Trait_Plot(mediandf = trait_dat.median, vardf = trait_dat.var, 
           stateID = 1, trait_to_plot = "d_min", add_mean=true)

# Geno_Plot(mean, stateID, "trait name")
Geno_Freq_Plot(freqdf=trait_dat.median, stateID=1, geno_names="g_1")