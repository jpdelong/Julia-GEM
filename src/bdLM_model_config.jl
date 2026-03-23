"""
This is the model configuration file. There are 3 parts to this document.
1. Initia state and parameter choices 
2. Instantiate 
"""

# ======================================================================
# function to sample from a lognormal distribution 
using Distributions
include("functions/PickTrait.jl") 

# Set seed for reproducibility
using Random
Random.seed!(42)

# ======================================================================
#                    INITIAL STATE AND PARAMETERS
# ======================================================================

N_init = [ 10 ] # initial abundances (the array takes integer values)

# bd-logistic parameters distribution mu and sigma
# maximum birth
b_max_mu = 4.0
b_max_sigma = 0.0

# minimum death
d_min_mu = 1.0
d_min_sigma = 0.0

# density dependence of birth
b_s_mu = 0.0012
b_s_sigma = 0.0

# density dependence of death
d_s_mu = 1e-5
d_s_sigma = 0.0


# randomly draw one sample from the lognormal distrinution
# (see PickTrait.jl for MU and SIGMA transformations )
b_max = PickTrait(b_max_mu, b_max_sigma)  # max birth
d_min = PickTrait(d_min_mu, d_min_sigma) # min death
b_s = PickTrait(b_s_mu, b_s_sigma) # density dependence of birth
d_s = PickTrait(d_s_mu, d_s_sigma) # density dependence of death

# calculate initial constant 
r_max = b_max-d_min
K = floor(((b_max - d_min)/(b_s + d_s)))

param_vect = [b_max, d_min, b_s, d_s] # parameter vector 
par_names = ["b_max", "d_min", "b_s", "d_s"] # parameter names
no_state = length(N_init) # number of states (do not edit)
no_param = length(param_vect) # number of params (do not edit)

constant_vect = []
# ===================================================================
# mapping arrays 
# nrow = state; ncol = param
state_par_match = [1 1 1 1] # matching parameters to state
# nrow = state, ncol = genotype 
state_geno_match = [1 0 0] # matching genotype to state
geno_names = ["g_1", "g_2", "g_3"] # genotype name
no_columns = no_param + 1 + size(state_geno_match, 2) # do not edit

# ===================================================================
# simulation design choices
GEM_ver = ["ver1"]#, "ver2", "ver3", "ver4"] # number of GEM versions
# nrow: state; ncol  = GEM ver
h2 = [ 0.0 ]#0.2 0.2 0.2 ] # narrow sense heritability
# cv = array{state ID, length(param), GEM ver}
"""
Note: The first stack is for GEM ver 1; typically reserved for "no-evolution". All elements are set to 0.0
In stack 2, set cv value for parameters corresponding to each state. 
You can mirror the dimensions of state_parameter_match matrix defined above. 
1 -> cv value
0 -> n/a for this state
"""
cv = cat([ 0.0 0.0 0.0 0.0], #ver 1
        # [ 0.2 0.0 0.0 0.0], #ver 2
        # [ 0.0 0.2 0.0 0.0], #ver 3
        # [ 0.2 0.2 0.0 0.0], #ver 4
           dims=3) 

# =================================================================== 

# replicate and time
num_rep = 2 # number of replicates
t_max = 20.0 # maximum time 
min_time_step_to_store = 0.1 # time points when data is stored
stand_time = range(0, t_max, step = min_time_step_to_store) # standadized number of time points across all runs (do not edit)
stand_time = collect(stand_time) # do not edit
num_time_steps = length(stand_time) # do not edit

# ===================================================================
# storage containers (do not edit)
pop_stand_out_all = fill(NaN, no_state, num_time_steps, num_rep, length(GEM_ver))
x_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_state,num_rep, length(GEM_ver))
x_var_stand_out_all = fill(NaN, no_columns-1,num_time_steps, no_state,num_rep, length(GEM_ver))


# ======================================================================
#                             INSTANTIATE 
# Only make changes to this block if you make changes to the variable 
# names that are used to instantiate the struct. 
# ======================================================================
""" 1. Instantiate the initial population state """
N0 = InitState(N_init) # Vector{Int}

""" 2. Instantiate ModelParVector """
model_par_vect = ModelParVector(
    param_vect # Vector{Float64}
)

""" 3. Instantiate DesignChoices """
design_choices = DesignChoices(
    h2, # Matrix{Float64}
    cv,  # Matrix{Float64}
    GEM_ver # Vector{String}
)

""" 4. Instantiate SimulationMap """
mappings = SimulationMaps(
    state_par_match, # Matrix{Int}
    state_geno_match, # Matrix{Int}
    par_names, # Vector{String}
    geno_names # Vector{String}
)

"""  5. Instantiate SimulationParameters """
sim_params = SimulationParameter(
    no_state, # Int
    no_param, # Int
    no_columns,  # Int 
    num_time_steps, #Int
    num_rep, # Int
    t_max, # Float64
    min_time_step_to_store # Float64 
    )
""" 6. Instantiate output container """
sim_output = GEMOutput(
    pop_stand_out_all, # Array{Float64, 4}
    x_stand_out_all, # Array{Float64, 5}
    x_var_stand_out_all # Array{Float64, 5}
    )
    
""" 7. Instantiate Constants """
gem_const_vect = GEMConstant(
   constant_vect
   ) 
