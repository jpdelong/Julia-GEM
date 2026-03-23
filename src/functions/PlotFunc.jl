
function Trait_Plot(; mediandf::DataFrame, vardf::DataFrame, 
    stateID::Int64, trait_to_plot::String, add_mean::Bool)

    #pick the right GEM version number and state ID
    mediandf = mediandf[mediandf.state_ID .== stateID, :]
    vardf = vardf[vardf.state_ID .== stateID, :]
    #pick the trait wanted
    mediantemp = mediandf[:, [:time, :rep, :GEM_ver]]
    mediantemp2 = mediandf[:,trait_to_plot]
    median2plot = hcat(mediantemp, mediantemp2)
    #grab the right var
    vartemp2 = vardf[:,trait_to_plot]
    dftemp = hcat(median2plot, vartemp2, makeunique=true)
    rename!(dftemp, :x1 => :median, :x1_1 => :var)
    upper_bound = dftemp.median .+ dftemp.var
    lower_bound = dftemp.median .- dftemp.var
    df2plot = hcat(dftemp, upper_bound, lower_bound, makeunique = true)
    rename!(df2plot,:x1 => :upper_bound, :x1_1 => :lower_bound )
    
    # Calculate mean and std across replicates for each time point and GEM_ver
    mean_df = combine(groupby(df2plot, [:time, :GEM_ver])) do subdf
        (mean_median = mean(subdf.median),
         std_median = std(subdf.median))
    end
    mean_df.mean_lower = mean_df.mean_median .- mean_df.std_median
    mean_df.mean_upper = mean_df.mean_median .+ mean_df.std_median
    
    median_plot = data(df2plot) * mapping(
        :time,
        :median,
    #color = :rep,
        group = :rep => nonnumeric,
        row = :GEM_ver => nonnumeric
    ) *
        visual(Lines, color = (:grey, 0.5))
    var_plot = data(df2plot) * mapping(
        :time,
        :upper_bound,
        :lower_bound,
        group = :rep => nonnumeric,
        row = :GEM_ver => nonnumeric
    )    * visual(Band, color = :lightgrey, alpha=0.5)
    
    mean_band = data(mean_df) * mapping(
        :time,
        :mean_lower,
        :mean_upper,
        row = :GEM_ver => nonnumeric
    ) * visual(Band, color = (:blue, 0.1))
    
    mean_line = data(mean_df) * mapping(
        :time,
        :mean_median,
        row = :GEM_ver => nonnumeric
    ) * visual(Lines, color = :blue, linewidth = 2)
    
    if add_mean
        tempplot = data(df2plot) * (median_plot + var_plot + mean_band + mean_line)
        finalplot = draw(tempplot, axis=(xlabel="time", ylabel="median trait value"))
    else
        tempplot = data(df2plot) * (median_plot + var_plot)
        finalplot = draw(tempplot, axis=(xlabel="time", ylabel="median trait value"))
    end
return finalplot
end


# =====================================================

function Geno_Freq_Plot(; freqdf::DataFrame, 
                        stateID::Int64, 
                        geno_names::String )

    #pick the right GEM version number and state ID
    freqdf = freqdf[freqdf.state_ID .== stateID, :]
   
    #pick the trait wanted
    freqtemp = freqdf[:, [:time, :rep, :GEM_ver]]
    freqtemp2 = freqdf[:,geno_names]
    freq2plot = hcat(freqtemp, freqtemp2)

    rename!(freq2plot, :x1 => :freq)

    freq_plot = data(freq2plot) * mapping(
        :time,
        :freq,
        #color = :rep,
        group = :rep => nonnumeric,
        row = :GEM_ver => nonnumeric
    ) * visual(Lines, color = :grey)
        

    finalplot = draw(freq_plot, axis=(xlabel="time", ylabel="genotype frequency"))

    return finalplot
end

#====================================================# 
function Pop_Plot(;pop_data::DataFrame, stateID::Int64, add_mean::Bool)
    pop_time = pop_data[pop_data.state_ID .== stateID, :]
   
    rep_col_name = names(pop_time)[4:size(pop_time)[2]]

    for col in rep_col_name
        col_data = pop_time[:, col]
        first_zero = findfirst(==(0.0), col_data)
        if !isnothing(first_zero) && first_zero < length(col_data)
            pop_time[first_zero+1:end, col] .= NaN
        end
    end


    pop_stack = stack(pop_time, 
                    rep_col_name)

    if add_mean  
        # Calculate average and standard deviation across replicates
        pop_time.mean = mean.(eachrow(pop_time[:, rep_col_name]))
        pop_time.std = std.(eachrow(pop_time[:, rep_col_name]))
        pop_time.lower = pop_time.mean .- pop_time.std
        pop_time.upper = pop_time.mean .+ pop_time.std
        # Stack the mean separately
        mean_stack = stack(pop_time, [:mean])
        pop_plot = data(pop_stack) * mapping(
            :time,
            :value,
            group = :variable,
            row = :GEM_ver => nonnumeric 
            ) *
            visual(Lines, color = (:grey, 0.5)) +
            data(pop_time) * mapping(
            :time,
            :lower,
            :upper,
            row = :GEM_ver => nonnumeric
            ) *
            visual(Band, color = (:blue, 0.1)) +
            data(mean_stack) * mapping(
            :time,
            :value,
            row = :GEM_ver => nonnumeric
            ) *
            visual(Lines, color = :blue, linewidth = 2)
        draw(pop_plot, axis=(xlabel="time", ylabel="population abundances"))    
    else
        pop_plot = data(pop_stack) * mapping(
            :time,
            :value,
            #color = :replicate,
            group = :variable,
            row = :GEM_ver => nonnumeric 
            ) *
            visual(Lines, color = :grey)
        draw(pop_plot, axis=(xlabel="time", ylabel="population abundances"))    
    end
end