"""
This is the model definition file. There are 3 parts to this document.
1. Your model definition.
2. Birth-death terms, and events function.  
3. struct definitons for the model, parameters, and simulation.
"""

# ======================================================================
#           1.   MODEL 
# ======================================================================
"""
Define your model here:
dNdt := the derivative value
N := integer state value for basal resource
C := integer state value for intermediate consumer
P := integer state value for top predator
p := parameters: (b_max, d_min, b_s, a_pc, a_pn, a_cn, h_pc, h_pn, h_cn, w_p, w_c, e_pc, e_pn, e_cn, d_p, d_c)
t := time

ODE: 
P_denominator = 1 + a_pc*h_pc*C + a_pn*h_pn*N + w_p*(P-1)
dPdt = (e_pc*a_pc*C + e_pn*a_pn*N)/P_denominator - d_p*P
dCdt = e_cn*a_cn*N / (1 + a_cn*h_cn*N + w_c*(C-1)) - a_pc*C/P_denominator - d_c*C
dNdt = (b_max-b_s*N)*N - (d_min1+d_s*N)*N - a_cn*N/(1 + a_cn*h_cn*N + w_c*(C-1)) - a_pn*N/P_denominator

# calculate constants:
r_max = b_max-d_min
K = floor(vec((b_max - d_min)/(b_s + d_s))[1])
"""

# ======================================================================
#          2. BIRTH-DEATH FUNCTIONS
#
# Note: Symmetric ODE systems can be generalized using linear algebraic matrix notation, 
# but care must be taken regarding dimensionality of the matrix operations. 
# If your birth or death function is outputting an array corresponding to multiple states, 
# then be extra careful while using the output downstream
# ======================================================================
"""
Define the birth and death terms for each state below.
"""
function Birth_N(b_max::Float64, b_s::Float64, N::Vector{Int})
    b_N = max((b_max - b_s*N[1]),0)*N[1] # logistic density dependent birth
    #=println(max((b_max - b_s*N[1]),0))
    println(b_N)
    println(b_max)
    println(b_s)
    println(N[1])=#
    return b_N
end

function Death_N(d_min::Float64, d_s::Float64, a_cn::Float64, h_cn::Float64, w_c::Float64, a_pc::Float64, h_pc::Float64, a_pn::Float64, h_pn::Float64, w_p::Float64, N::Vector{Int})
    d_NN = (d_min + d_s*N[1])*N[1] # logistic density dependent death
    d_NC = a_cn*N[1]*N[2]/(1 + a_cn*h_cn*N[1] + w_c*(N[2]-1)) # death by C
    d_NP = a_pn*N[1]*N[3]/(1 + a_pn*h_pn*N[1] + a_pc*h_pc*N[2] + w_p*(N[3]-1)) # death by P
    death_N = d_NN + d_NC + d_NP # all the deaths
    return death_N
end

function Birth_C(e_cn::Float64, a_cn::Float64, h_cn::Float64, w_c::Float64, N::Vector{Int})
    b_C = e_cn*a_cn*N[1]*N[2]/(1 + a_cn*h_cn*N[1] + w_c*(N[2]-1)) # birth of C by eating N
    return b_C
end

function Death_C(a_pc::Float64, h_pc::Float64, a_pn::Float64, h_pn::Float64, w_p::Float64, d_c::Float64, N::Vector{Int})
    d_CP = a_pc*N[2]*N[3]/(1 + a_pn*h_pn*N[1] + a_pc*h_pc*N[2] + w_p*(N[3]-1)) # death by P
    d_CC = d_c*N[2] # background death of C
    death_C = d_CC + d_CP # all the deaths
    return death_C
end

function Birth_P(e_pc::Float64, e_pn::Float64, a_pc::Float64, h_pc::Float64, a_pn::Float64, w_p::Float64, N::Vector{Int})
    b_P = (e_pc*a_pc*N[2] + e_pn*a_pn*N[1])*N[3] / (1 + a_pn*h_pn*N[1] + a_pc*h_pc*N[2] + w_p*(N[3]-1)) # birth of P
    return b_P
end

function  Death_P(d_p::Float64, N::Vector{Int})
    death_P = d_p*N[3] # death of P
    return death_P
end

function Event_Terms(param_next::Matrix{Float64}, const_vect::Any, N::Vector{Int})
    b_max = param_next[1,1] # max birth
    d_min = param_next[1,2] # min death
    b_s = param_next[1,3] # density dependence of birth
    d_s = param_next[1,4] # density dependence of death 
    a_cn = param_next[2,5] # space clearance rate of consumers on prey
    h_cn = param_next[2,6] # handling time of consumers on prey
    w_c = param_next[2,7] # interference of consumer
    e_cn = param_next[2,8] # conversion efficiency of consumer
    d_c = param_next[2,9] # death of consumer
    a_pc = param_next[3,10] # space clearance rate of predators on consumers
    a_pn = param_next[3,11] # space clearance rate of predators on prey
    h_pc = param_next[3,12] # handling time of predators on consumers
    h_pn = param_next[3,13] # handling time of predators on prey
    w_p = param_next[3,14] # interference of predators
    e_pc = param_next[3,15] # conversion efficiency of predators on consumers
    e_pn = param_next[3,16] # conversion efficiency of predators on prey
    d_p = param_next[3,17] # death of predators

    bN = Birth_N(b_max, b_s, N)
    dN = Death_N(d_min, d_s, a_cn, h_cn, w_c, a_pc, h_pc, a_pn, h_pn, w_p, N)
    bC = Birth_C(e_cn, a_cn, h_cn, w_c, N)
    dC = Death_C(a_pc, h_pc, a_pn, h_pn, w_p, d_c, N)
    bP = Birth_P(e_pc, e_pn, a_pc, h_pc, a_pn, w_p, N)
    dP = Death_P(d_p, N)

    return bN, dN, bC, dC, bP, dP
end

# ======================================================================
#          3.  STRUCTURES
# only make changes to the following block if you are addint new struct
# ======================================================================
"""
    1. InitState
Initial state struct
"""
mutable struct InitState
    N::Vector{Int} # initial state abundances vector
end

"""
    2. ModelParVector
"""
struct ModelParVector{T}
    param_init::Vector{T} # initial parameter vector
end

"""
    3. DesignChoices{T}
Eco-evo choices
"""

struct DesignChoices
    h2::Array{Float64} # narrow sense heritability
    cv::Array{Float64} # coefficient of variation
    GEM_ver::Vector{String} # GEM version
end

"""
    4. SimulationMapping
"""
struct SimulationMaps
    state_par_match::Matrix{Int} # match parameters to state 
    state_geno_match::Matrix{Int} # match genotype to state
    par_names::Vector{String} # parameter names vector 
    geno_names::Vector{String} # genotype names vector 
end

"""
    5. SimulationParameters
"""
struct SimulationParameter
    no_state::Int # total number of states 
    no_param::Int # total number of parameter 
    no_columns::Int # total number of columns 
    num_time_steps::Int # standardized number of time steps 
    num_rep::Int # number of replication per GEM version 
    t_max::Float64 # total time to run simulation for 
    min_time_step_to_store::Float64 # time points to store data at 
end

"""
    6. GEMOutput: storage for all replicates and versions of GEM
"""
struct GEMOutput
    pop_stand_out_all::Array{Float64, 4} # population data storage 
    x_stand_out_all::Array{Float64, 5} # trait median & genotype frequency data storage 
    x_var_stand_out_all::Array{Float64, 5} # trait variance data storage
end

"""
7. CONSTANTS
"""
struct GEMConstant{T} 
    const_vect::Vector{T}
end
