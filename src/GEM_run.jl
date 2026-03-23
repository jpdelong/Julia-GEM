""" 
This file provides the space to run your GEM model.
Before you start, please read README.txt
"""
# run only once the first time
# include("install_pkgs.jl")

# load all required packages 
include("functions/Packages.jl")

# load model definition  
#include("bdLM_model_definition.jl") 

include("bd_IGM_model_definition.jl") 

#include("RF_model_def_nointerf.jl") 

# load model configuration
include("bd_IGM_model_config.jl")

# load all functions 
include("functions/GEM_Functions.jl")

# run the GEM simulation
using Dates
start_t = now()

run_sim = GEM_sim(
                  N0, # initial state
                  model_par_vect, # model parameters in a vector
                  gem_const_vect,
                  design_choices, # evolution decisions
                  mappings, # parameter-to-state mapping
                  sim_params,# simulation parameters
                  sim_output, # output containers
                  verbose=true # show time on console
                  ) #

end_t = now()




# output: 
# Tuple{Array{Float64, 4}, Array{Float64, 5}, Array{Float64, 5}}

# dataframe for population time series
pop_dat = run_sim.pop_df

# 2 dataframes: mean and variance
trait_dat = run_sim.trait_df

# accessing the two trait dataframes
# trait mean dataframe:
trait_dat.median 

# trait variance dataframe
trait_dat.var

# Pop_Plot(pop data, stateID)
p = Pop_Plot(pop_data=run_sim.pop_df, stateID= 3, add_mean=true)



# Trait_Plot(mean, var, stateID, "trait name")
Trait_Plot(mediandf=trait_dat.median, vardf=trait_dat.var, 
stateID=3, trait_to_plot="a_pn", add_mean=false)


# Geno_Plot(mean, stateID, "trait name")
Geno_Freq_Plot(trait_dat.median, 1, "g_1")

```
- check density dependence forms
- run ODE for model with interference 
 
```