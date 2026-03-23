# ========================================== #
#		  	  FUNCTION PICK EVENTS           #
# ========================================== #

function PickEvent(terms::Vector{Float64}, no_state::Int)
	
	terms = reshape(terms, 1, length(terms)) #reshape(ele_to_reshape, new_row, new_col)
	terms = transpose(reshape(terms, 2, no_state ))
    
	temp = hcat(1:no_state,terms) # in this temp struc, first col is N, col2 = birth; col3 = death
	state_ids = temp[:, 1]
	rates = temp[:, 2:end]

	nz_rate_pos = findall(rates .!= 0 .&& .!isnan.(rates)) #findall returns a vector of tuples (row, col) where the condition is true
	nz_rate_vals = rates[nz_rate_pos] 
	c_sum = cumsum(nz_rate_vals)  
	pie_slices = c_sum ./ c_sum[end] #generated weighted slices b/w 0-1
	r_num = rand()

	less_than = r_num .< pie_slices  #BitMatrix
	event_index = findfirst(==(1), less_than)

	picked_event = nz_rate_pos[event_index]
	row = picked_event[1]  # state ID
	col = picked_event[2]  # event type (1 = birth, 2 = death) 
	
	return (c_sum = c_sum, event=col, state=row)
	

end


#=
	c_sum = cumsum(terms, dims=2)  
	pie_slices = c_sum ./ c_sum[end] #generated weighted slices b/w 0-1
	r_num = rand()
	
	
	row = -1
	col = -1
	for r = 1:size(event_mat,1) #want to go by row
		col_ind_in_row_r = findfirst(==(1),event_mat[r,:])
		if col_ind_in_row_r != nothing
			row = r
			col = col_ind_in_row_r
			break
		end
	end
=#