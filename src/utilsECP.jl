
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

"""
Vectorize ECP

Input:  ecp - dictionary consisting of filtrations as keys and contributions as values.
        number_of_steps - Tuple containing number of steps for each filtration parameter.
Return: FILT_PARAM_NUMBER - dimensional array containing vectorised ECP
        Shape of the array is (number_of_steps[i] + 1)_{i=1...FILT_PARAM_NUMBAR}
"""
function py_vectorize_ecp(ecp::Dict{Any, Any}, number_of_steps)
    for filtration in keys(ecp)
        type_of_filt = typeof(filtration)
        return vectorize_ecp(convert(Dict{type_of_filt, Int64}, ecp), number_of_steps)
    end
end

function read_input_from_perseus_file(file_path)
    file_content = readlines(file_path)
    # Determine dimensions of cubical complex
    dim = parse(Int64, file_content[1])
    sizes = Tuple(parse(Int64, file_content[i]) for i in 2:(dim+1))
    # Determine type of multifiltrations and number of parameters
    first_line_filtration = split(file_content[2 + dim])
    PARAM_NUM = length(first_line_filtration)
    FILT_PARAM_TYPE =  tryparse(Int64, first_line_filtration[1]) == nothing ? Float64 : Int64
    output_array = Array{NTuple{PARAM_NUM, FILT_PARAM_TYPE}}(undef, sizes...)
    # fill array with max dim cells filtration
    file_index = 1 + 1 + dim
    for k in Base.Iterators.product((1:sizes[i] for i in 1:length(sizes))...)
        line_filtration_str = split(file_content[file_index])
        output_array[k...] = NTuple{PARAM_NUM, FILT_PARAM_TYPE}(parse(FILT_PARAM_TYPE, line_filtration_str[i]) for i in 1:PARAM_NUM)
    end

    return output_array
end

function is_less_fil(fil1, fil2)
    for i in 1:length(fil1)
        if fil1[i] >= fil2[i]
            return false
        end
    end
    return true
end

function is_less_eq_fil(fil1, fil2)
    for i in 1:length(fil1)
        if fil1[i] > fil2[i]
            return false
        end
    end
    return true
end

function get_ecp_at_filtration(ecp_vec, filt_value)
    acc = 0
    for i in 1:length(ecp_vec)
        if is_less_fil(ecp_vec[i][1], filt_value)
            acc += ecp_vec[i][2]
        end
    end
    return acc
end


function index_to_filtrations(index, filtration_values)
    return Tuple(filtration_values[i][index[i]] for i in 1:length(index))
end

function norm(ecp_diff::Dict{NTuple{PARAM_NUM, T}, Int64} ) where {PARAM_NUM, T}
    filtratin_values = Array{Array{T}}(undef, PARAM_NUM)
    for i in 1:PARAM_NUM
        filt_vals = Set( filt_tuple[i] for filt_tuple in keys(ecp_diff))
        filtratin_values[i] = collect(filt_vals)
        sort!(filtratin_values[i])
    end

    norm_acc = 0.0
    ecp_diff_vec = collect(ecp_diff)
    for k in Base.Iterators.product((2:length(filtratin_values[i]) for i in 1:length(filtratin_values))...)
        tup_index = Tuple(k)
        area = 1.0
        for i in 1:length(tup_index)
            area *= abs(filtratin_values[i][tup_index[i] - 1] - filtratin_values[i][tup_index[i]])
        end
        norm_acc += area * abs(get_ecp_at_filtration(ecp_diff_vec, index_to_filtrations(k, filtratin_values)))
    end

    return norm_acc
end

function calculate_diff(ecp1::Dict{T, Int64}, ecp2::Dict{T, Int64}) where {T}
    ecp_diff = deepcopy(ecp1)
    for filtration in keys(ecp2)
        ecp_diff[filtration] = get(ecp_diff, filtration, 0) - ecp2[filtration]
    end
    return ecp_diff
end

function prune_ecp!(ecp_diff, fmin, fmax)
    filt_function = contribution -> contribution[2] != 0 && is_less_eq_fil(fmin, contribution[1]) && is_less_eq_fil(contribution[1], fmax) 
    filter!(filt_function, ecp_diff)
end

function distance_ecp(ecp1, ecp2, fmins, fmaxs)
    ecp_diff = calculate_diff(ecp1, ecp2)
    prune_ecp!(ecp_diff, fmins, fmaxs)
    ecp_diff[fmins] = get(ecp_diff, fmins, 0)
    ecp_diff[fmaxs] = get(ecp_diff, fmaxs, 0)
    return norm(ecp_diff)
end
