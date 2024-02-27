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


function compute_contributions(image::Array{NTuple{PARAM_NUM, Float64}}) where PARAM_NUM
    # compute contributions of all cells,
    # starting from bottom left
    # uses lowert star filtration
    
    contributions = Dict{NTuple{PARAM_NUM, Float64}, Int64}()
    
    input_size = size(image)
    for i in 2:input_size[1]
        for j in 2:input_size[2]
            for k in 2:input_size[3]
            # lets track all the contributions
            # from cell i,j,k
                        
            # itself, 3d cell
            f = image[i, j, k]
            contributions[f] = get(contributions, f, 0) + 1
                        
            # 0d cell SW
            f = calc_filtration_for_0_dim(i, j, k, image)
            contributions[f] = get(contributions, f, 0)  + 1
            
            # 1d cell i
            f = calc_filtration_for_1_dim_depth(i, j, k, image)
            contributions[f] = get(contributions, f, 0) - 1
            
            # 1d cell j
            f = calc_filtration_for_1_dim_vertical(i, j, k, image)
            contributions[f] = get(contributions, f, 0) - 1

            # 1d cell k
            f = calc_filtration_for_1_dim_horisontal(i, j, k, image)
            contributions[f] = get(contributions, f, 0) - 1
            end
        end
    end
    
    # remove the contributions that are 0
    to_del = []
    for key in keys(contributions)
        if contributions[key] == 0
            push!(to_del, key)
        end
    end
    for key in to_del
        delete!(contributions, key)
    end

    return collect(contributions)
end

const SIZE = 500
img = Array{Tuple{Float64, Float64}}(undef, SIZE, SIZE, SIZE)
for i in 1:SIZE
    for j in 1:SIZE
        for k in 1:SIZE
            img[i, j, k] = (rand(), rand())
        end
    end
end

final = compute_contributions(img)
print(final)