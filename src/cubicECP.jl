import ImageFiltering

function mins_for_tuples(tuples)
    N = size(tuples[1], 1)
    return NTuple{N, Float64}(minimum(tuples[j][i] for j in 1:length(tuples)) for i in 1:N)
end

function calc_filtration_for_0_dim(i, j, k, image)
    return mins_for_tuples((image[i-1, j-1, k-1],  image[i, j-1, k-1], image[i, j, k-1], image[i-1, j, k-1],
                          image[i-1, j-1, k],  image[i, j-1, k] , image[i, j, k], image[i-1, j, k]))
end

function calc_filtration_for_1_dim_horisontal(i, j, k, image)
    return mins_for_tuples((image[i-1, j-1, k-1],  image[i, j-1, k-1], image[i, j, k-1], image[i-1, j, k-1]))
end

function calc_filtration_for_1_dim_vertical(i, j, k, image)
    return mins_for_tuples((image[i-1, j-1, k-1], image[i, j-1, k-1], image[i-1, j-1, k],  image[i, j-1, k]))
end

function calc_filtration_for_1_dim_depth(i, j, k, image)
    return mins_for_tuples((image[i-1, j-1, k-1], image[i-1, j, k-1], image[i-1, j-1, k], image[i-1, j, k]))
end


"""
Compute cotributions of cells to ECP.
Input: 2-dimensional array. Each entity is tuple describing filtration of high dimensional cell
Return: Dictionary describing contributions at filtration
      keys - filtration values, values - contributions
"""
function compute_contributions_2d(image::Array{NTuple{PARAM_NUM, T}}) where {PARAM_NUM, T}
    # compute contributions of all cells,
    # starting from bottom left
    # uses lowert star filtration
    INFINITY_T = Tuple(typemax(T) for i in 1:PARAM_NUM)
    padded_image = ImageFiltering.padarray(image, ImageFiltering.Fill(INFINITY_T, (1,1), (1,1)))
    contributions = Dict{NTuple{PARAM_NUM, T}, Int64}()
    
    input_size = size(padded_image)
    for i in 1:input_size[1] - 1
        for j in 1:input_size[2] - 1
            # lets track all the contributions
            # from cell i,j
                        
            # itself, 2d cell
            f = padded_image[i,j]
            contributions[f] = get(contributions, f, 0) + 1
                        
            # 0d cell SW
            f = mins_for_tuples((padded_image[i,j], padded_image[i-1,j-1], padded_image[i-1,j], padded_image[i,j-1]))
            contributions[f] = get(contributions, f, 0) + 1
            
            # 1d cell W
            f = mins_for_tuples((padded_image[i,j], padded_image[i,j-1]))
            contributions[f] = get(contributions, f, 0) - 1
            
            # 1d cell S
            f = mins_for_tuples((padded_image[i,j], padded_image[i-1,j]))
            contributions[f] = get(contributions, f, 0) - 1
        end
    end

    delete!(contributions, INFINITY_T)
    
    return contributions
end


"""
Compute cotributions of cells to ECP.
Input: 3-dimensional array. Each entity is tuple describing filtration of high dimensional cell
Return: Dictionary describing contributions at filtration
      keys - filtration values, values - contributions
"""
function compute_contributions_3d(image::Array{NTuple{PARAM_NUM, T}}) where {PARAM_NUM, T<:Number}
    # compute contributions of all cells,
    # starting from bottom left
    # uses lowert star filtration
    INFINITY_T = Tuple(typemax(T) for i in 1:PARAM_NUM)
    padded_image = ImageFiltering.padarray(image, ImageFiltering.Fill(INFINITY_T, (1,1,1), (1,1,1)))
    contributions = Dict{NTuple{PARAM_NUM, T}, Int64}()

    input_size = size(padded_image)
    for i in 1:input_size[1]-1
        for j in 1:input_size[2]-1
            for k in 1:input_size[3]-1
            # lets track all the contributions
            # from cell i,j,k
                        
            # itself, 3d cell
            f = padded_image[i, j, k]
            contributions[f] = get(contributions, f, 0) - 1
                        
            # 0d cell SW
            f = calc_filtration_for_0_dim(i, j, k, padded_image)
            contributions[f] = get(contributions, f, 0)  + 1
            
            # 1d cell i
            f = calc_filtration_for_1_dim_depth(i, j, k, padded_image)
            contributions[f] = get(contributions, f, 0) - 1
            
            # 1d cell j
            f = calc_filtration_for_1_dim_vertical(i, j, k, padded_image)
            contributions[f] = get(contributions, f, 0) - 1

            # 1d cell k
            f = calc_filtration_for_1_dim_horisontal(i, j, k, padded_image)
            contributions[f] = get(contributions, f, 0) - 1

            # 2d cells
            f = mins_for_tuples((padded_image[i-1, j, k], padded_image[i, j, k]))
            contributions[f] = get(contributions, f, 0) + 1

            f = mins_for_tuples((padded_image[i, j-1, k], padded_image[i, j, k]))
            contributions[f] = get(contributions, f, 0) + 1

            f = mins_for_tuples((padded_image[i, j, k-1], padded_image[i, j, k]))
            contributions[f] = get(contributions, f, 0) + 1
            end
        end
    end
    
    delete!(contributions, INFINITY_T)
    
    return contributions
end
