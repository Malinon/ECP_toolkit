
function calc_min_max_filtration(ecp, i)
    return (minimum(x->x[i], keys(ecp)), maximum(x->x[i], keys(ecp)))
end

function fill_ecp_deltas_array(ecp, ecp_array_out, step_sizes, min_max_filtrations)
    for filtration in keys(ecp)
        ecp_array_out[(1 + Int64(trunc((filtration[i] - min_max_filtrations[i][1]) / step_sizes[i])) for i in 1:length(filtration))...] += ecp[filtration]
    end
end

function calc_ecp_at_point(ecp_array_out::Array{Int64}, ind)
    acc = 0 # accumulator for ecp at ind
    ind_ctr_sum = sum(Tuple(ind))
    # Use inclusion-exclusion princimmple to calculate ECP at ind
    for k in Base.Iterators.product((max(1, ind[i]-1):ind[i] for i in 1:length(ind))...)
        curr_ctr_sum = sum(k)
        if iseven(ind_ctr_sum - curr_ctr_sum)
            if curr_ctr_sum == ind_ctr_sum
                acc += ecp_array_out[ind]
            else
                acc -= ecp_array_out[k...]
            end
        else
            acc += ecp_array_out[k...]
        end
    end
    return acc
end

function agregate_ecp_deltas(ecp_array_out)
    for i in CartesianIndices(ecp_array_out)
        ecp_array_out[i] = calc_ecp_at_point(ecp_array_out, i)
    end
end

function replace_zero_with_one(x)
    if x == 0.0
        return 1.0
    else
        return x
    end
end

"""
Vectorize ECP

Input:  ecp - dictionary consisting of filtrations as keys and contributions as values.
        number_of_steps - Tuple containing number of steps for each filtration parameter.
Return: FILT_PARAM_NUMBER - dimensional array containing vectorised ECP
        Shape of the array is (number_of_steps[i] + 1)_{i=1...FILT_PARAM_NUMBAR}
"""
function vectorize_ecp(ecp::Dict{NTuple{PARAM_NUM, T}, Int64}, number_of_steps) where {PARAM_NUM, T}
    # Find minimum and maximum value of filtration
    min_max_filtrations = NTuple{PARAM_NUM, Tuple{T, T}}(calc_min_max_filtration(ecp,  i) for i in 1:PARAM_NUM)
    # Calculate step size for each parameter of filtration
    step_sizes = NTuple{PARAM_NUM, Float64}(replace_zero_with_one(min_max_filtrations[i][2] - min_max_filtrations[i][1]) / number_of_steps[i]
            for i in 1:PARAM_NUM)
    # calculate shape of an output array
    output_dims = Tuple(number_of_steps[i] + 1 for i in 1:PARAM_NUM)
    ecp_array_out = zeros(Int64, output_dims)
    fill_ecp_deltas_array(ecp, ecp_array_out, step_sizes, min_max_filtrations)
    agregate_ecp_deltas(ecp_array_out)
    return ecp_array_out
end
