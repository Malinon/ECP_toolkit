
function calc_min_max_filtration(ecp, i)
    return (minimum(x->x[i], keys(ecp)), maximum(x->x[i], keys(ecp)))
end


function fill_ecp_deltas_array(ecp, ecp_array_out, step_sizes, min_max_filtrations)
    for filtration in keys(ecp)
        index_filt = Tuple(1 + Int64(trunc((filtration[i] - min_max_filtrations[i][1]) / step_sizes[i])) for i in 1:length(filtration))
        ecp_array_out[index_filt...] += ecp[filtration]
    end
end

function calc_ecp_at_point(ecp_array_out::Array{Int64}, ind)
    acc = 0
    ind_ctr_sum = sum(Tuple(ind))
    # Use inclusion-exclusion princimmple to calculate ECP at ind
    for k in Base.Iterators.product((max(1, ind[i]-1):ind[i] for i in 1:length(ind))...)
        if iseven(ind_ctr_sum - sum(k))
            acc += ecp_array_out[k...]
        else
            acc -= ecp_array_out[k...]
        end
    end
    return acc
end


function agregate_ecp_deltas(ecp_array_out)
    for i in CartesianIndices(ecp_array_out)
        ecp_array_out[i] = calc_ecp_at_point(ecp_array_out, i)
    end
end

function vectorize_ecp(ecp, number_of_steps)
    FILT_PARAM_NUMBER = length(number_of_steps)
    min_max_filtrations = NTuple{FILT_PARAM_NUMBER, Tuple{Float64, Float64}}(calc_min_max_filtration(ecp,  i) for i in 1:FILT_PARAM_NUMBER)
    step_sizes = NTuple{FILT_PARAM_NUMBER, Float64}((min_max_filtrations[i][2] - min_max_filtrations[i][1]) / number_of_steps[i] for i in 1:FILT_PARAM_NUMBER)
    output_dims = Tuple(number_of_steps[i] + 1 for i in 1:FILT_PARAM_NUMBER)
    ecp_array_out = Array{Int64}(undef, output_dims)
    fill_ecp_deltas_array(ecp, ecp_array_out, step_sizes, min_max_filtrations)
    println(ecp_array_out)
    agregate_ecp_deltas(ecp_array_out)
    return ecp_array_out
end