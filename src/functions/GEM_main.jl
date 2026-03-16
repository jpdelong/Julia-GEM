"""
    1. Function to run one replication 
"""

function run_replicate(init_state::InitState,
                        mod_par_vect::ModelParVector,
                        gem_constants::GEMConstant,
                        dc::DesignChoices,
                        sim_map::SimulationMaps,
                        sim_par::SimulationParameter,
                        j::Int,
                        verbose::Bool)
    @unpack N = init_state
    @unpack param_init = mod_par_vect
    @unpack const_vect = gem_constants
    @unpack h2, cv = dc
    #@unpack state_geno_match, state_par_match, geno_par_match, which_par_quant = sim_map
    @unpack state_geno_match, state_par_match = sim_map
    @unpack t_max, no_state, no_columns, no_param, num_time_steps,min_time_step_to_store  = sim_par    

    # make a copy 
    t = 0.0
    N0 = copy(N)
    N = copy(N)
    params = copy(param_init) 

    # some internal container
    init_comm_mat =  Array{Float64}(fill(NaN, Int(sum(N0)), no_columns))
    pop_slice = Array{Int}(fill(0, no_state, num_time_steps))
    x_slice = fill(NaN, no_columns-1, num_time_steps, no_state)
    x_var_slice = fill(NaN, no_columns-1, num_time_steps, no_state)

    # store details of initial state 
    pop_slice[:,1] .= N0 #first col/time step gets the initial pop 
    
    """ Func Initiate Population """
    x_dist_init = InitiatePop(
        N0, 
        state_geno_match, 
        state_par_match, 
        init_comm_mat,
        params, 
        cv, 
        j)
        
        
    # Store initial state details
    extant_states = findall(N0 .!= 0)
    for ii in extant_states
        #(no_state)
        x_slice[:, 1, ii] = CalcMedian(ii, no_columns, no_param, x_dist_init)
        x_var_slice[1:no_param, 1, ii] = CalcVar(ii, no_param, x_dist_init)
    end
    

    # count up each individual for all states after the first sampling
    for jj = 1:no_state
        x = x_dist_init[:,1] #extract first col
        N[jj] = count(.==(jj), x)
    end
    
    x_dist = x_dist_init
    time_step_index = 2
    time_step = stand_time[time_step_index]
    
    while t < t_max && sum(N) > 0
        """ Func WhoIsNext """
        FindWhoNext = WhoIsNext(x_dist, no_state, no_columns, no_param, N, state_par_match, state_geno_match)
        param_next = FindWhoNext.param_next # FindWhoNext[1]; Using named tuple for cleaner access
        genotype_next = FindWhoNext.genotype_next # FindWhoNext[2]
        whosnext = FindWhoNext.whosnext # FindWhoNext[3]
       # @show param_next
       # @show whosnext

        """ Func Event Terms """
        terms = collect(Event_Terms(param_next, const_vect, N))
        @show terms

        """ Func Pick Event """
        picked_event = PickEvent(terms, no_state)
        c_sum = picked_event.c_sum # PickedEvent[1]
        event = picked_event.event # 1: birth, 2: death
        state = picked_event.state # state
       # @show picked_event

        if event == 1 # birth
            parent_traits = x_dist[Int(whosnext[state]), 2:no_columns] 
            #@show parent_traits   
                 """ Draw New Trait """     
            new_trait = DrawNewTraits(x_dist,parent_traits,h2,no_param,no_columns,state, j)
                   
            new_trait_row = hcat(state, new_trait)
            new_trait_row = hcat(new_trait_row[1],new_trait_row[2][1],new_trait_row[2][2])
            x_dist = vcat(x_dist, new_trait_row) 
        
        elseif event == 2  # death 
            # delete the individual by making a new copy of the matrix
            # without the row 
            x_dist = x_dist[1:size(x_dist, 1) .!= Int(whosnext[state]), :]
        end

        # Update abundances 
        for jj in 1:no_state
            N[jj] = sum(x_dist[:,1].== jj)
        end
        
        # while loop
        while t > time_step
            pop_slice[:,time_step_index] .= N  # assign current values to sliced standard times
            
            extant_states = findall(N .!= 0)
            for ii in extant_states
                x_slice[:,time_step_index,ii] = CalcMedian(ii,no_columns,no_param,x_dist)
                x_var_slice[1:no_param,time_step_index,ii] = CalcVar(ii, no_param, x_dist)
            end

            time_step_index +=  1 # advance to next standardized time
            time_step = stand_time[time_step_index]
        end
        #@show c_sum
        @show N
        # Advance time
        time_advance = exp(-1/c_sum[end])/(c_sum[end])
        if !isnan(time_advance) && time_advance > 0 
            t = t + time_advance
        else
            println("Time advance error. Check cumulative 
            sum. Stopped at time:\nT $t")
            #@show N
            break
        end 
     
        if verbose 
            @show t
        end
        
    end    


    if t > t_max
        println("Simulation reached t_max. Stopped at time:\nT $t")
        # store the last value of the replicate
         pop_slice[:,time_step_index] .= N  # assign current values to sliced standard times
           
            extant_states = findall(N .!= 0)
            for ii in extant_states
                x_slice[:,time_step_index,ii] = CalcMedian(ii,no_columns,no_param,x_dist)
                x_var_slice[1:no_param,time_step_index,ii] = CalcVar(ii, no_param, x_dist)
            end
    end 

        # Return the results for this replicate
    return (pop_time_series=pop_slice, trait_mom1 = x_slice, trait_mom2 = x_var_slice)
end

# ==================================================================

"""
    2. parallel GEM simulation function
"""

function GEM_sim(init_state::InitState,
                        mod_par_vect::ModelParVector,
                        gem_constants::GEMConstant,
                        dc::DesignChoices,
                        sim_map::SimulationMaps,
                        sim_par::SimulationParameter,
                        sim_op::GEMOutput;
                        verbose::Bool)
    
    @unpack N = init_state
    @unpack param_init = mod_par_vect
    @unpack const_vect = gem_constants
    @unpack h2, cv, GEM_ver = dc
#    @unpack state_geno_match, state_par_match, geno_par_match, which_par_quant = sim_map
    @unpack state_geno_match, state_par_match = sim_map
    @unpack num_rep, t_max, no_state, no_columns, no_param, num_time_steps,min_time_step_to_store  = sim_par    
    @unpack pop_stand_out_all, x_stand_out_all, x_var_stand_out_all = sim_op
    
    for j = 1:length(GEM_ver) # loop through the GEM versions

        # some internal containers
        pop_stand = zeros(no_state, num_time_steps, num_rep)
        x_stand = fill(NaN, no_columns - 1, num_time_steps, no_state, num_rep)
        x_var_stand = fill(NaN, no_columns - 1, num_time_steps, no_state, num_rep)
        
        rep_good = fill(false, num_rep)
        Threads.@threads for i = 1:num_rep # loop through the replicates
            @show Threads.threadid()  
            try          
            results = run_replicate(init_state, mod_par_vect, gem_constants, dc, sim_map, 
                sim_par,j, verbose)
            pop_slice = results.pop_time_series
            x_slice = results.trait_mom1
            x_var_slice = results.trait_mom2
            pop_stand[:, :, i] .= pop_slice
            x_stand[:, :, :, i] .= x_slice
            x_var_stand[:, :, :, i] .= x_var_slice
            catch e
                println("Error in replicate $i: $e on thread $(Threads.threadid())")
            end
        end
        
        pop_stand_out_all[:, :, :, j] .= pop_stand
        x_stand_out_all[:, :, :, :, j] .= x_stand
        x_var_stand_out_all[:, :, :, :, j] .= x_var_stand
    end
    
    # turn the multidimensional output arrays into long dataframes
    pop_out = make_pop_df_long(sim_op, sim_par, dc) # population time series dataframe
    trait_out = make_trait_df_long(sim_op, sim_par, dc, sim_map) # trait mean and var time seroes dataframes 
    # trait_out has two dataframes that can be accessed with named tuples: trait_out.mean and trait_out.var 

    return (pop_df = pop_out, trait_df = trait_out)
end

#===================================================================#

