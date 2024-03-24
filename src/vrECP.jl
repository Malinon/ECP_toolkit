using LinearAlgebra
using Statistics

function distance(t1, t2)
    size = length(t1)
    acc = 0.0
    for i in 1:size
        acc += (t1[i] - t2[i])^2
    end
    return sqrt(acc)
end

function create_local_graph(points, id_neigs_of_center_vectex, i, threshold, dbg=false)
    center_vertex = points[i]
    if dbg
        println(id_neigs_of_center_vectex)
    end
    # create the center graph as a list of lists
    # each list corespond to a node and contains its neighbours with distance
    # note, edges are always of the type
    # (i, j) with i<j in the ordering. 
    # This is to save space, we do not need the edge (j, i) 
    # In practice, we are building onlty the upper triangular part
    # of the adjecency matrix
    
    considered_graph = []
    
    mapped_center_vertex = []
    for (j, neigh) in enumerate(id_neigs_of_center_vectex)
        push!(mapped_center_vertex, (j + 1, distance(center_vertex, points[neigh])))
    end
    push!(considered_graph, mapped_center_vertex)
    
    for (j, neigh) in enumerate(id_neigs_of_center_vectex)
        j = j + 1
        if dbg
            println(j, " ", neigh)
        end
        neighbours_of_j = []
        for (z, other_neigh) in enumerate(id_neigs_of_center_vectex[j:end])
            z = z + j
            if dbg
                println("    ", z, " ", other_neigh)
            end
            dist = distance(points[neigh], points[other_neigh])
            if dist ≤ threshold
                push!(neighbours_of_j, (z, dist))
            end
        end
        push!(considered_graph, neighbours_of_j)
    end

    return considered_graph
end

function compute_ECP_single_vertex(considered_graph, vertex_filtrations, max_dimension, dbg=false)
    simplices_in_current_dimension = []
    filtration_of_those_simplices = Vector{Tuple{Float64, Float64}}()
    ECC = Dict{Tuple{Float64, Float64}, Int64}()
    ECC[(0, vertex_filtrations[1])] = 1
    
    if dbg
        println("vertex filtrations : ", vertex_filtrations)
    end
    
    for edge in considered_graph[1]
        push!(simplices_in_current_dimension, [1, edge[1]])
        this_filtration = (edge[2], max(vertex_filtrations[edge[1]], vertex_filtrations[1]))
        push!(filtration_of_those_simplices, this_filtration)
        ECC[this_filtration] = get(ECC, this_filtration, 0) - 1
    end
    
    last_neigh_of_current_vertex = length(considered_graph[1]) + 1
    if dbg
        println("last_neigh_of_current_vertex : ", last_neigh_of_current_vertex)
    end
    
    common_neighs = []
    
    for simplex in simplices_in_current_dimension
        the_other_vertex = simplex[2]
        neighs = Dict{Int64, Float64}()

        for neigh_of_other in considered_graph[the_other_vertex]
            push!(neighs, neigh_of_other[1] => max(neigh_of_other[2], considered_graph[1][the_other_vertex][2]))
        end
        push!(common_neighs, neighs)
    end
    
    neighs_of_vertices = Vector{Set{}}()
    for vertex_list in considered_graph
        push!(neighs_of_vertices, Set([v for v in vertex_list]))
    end

    dimension = 2
    dimm = 1

    while !isempty(common_neighs) && dimension<= max_dimension
        new_filtration_of_those_simplices = Vector{Tuple{Float64, Float64}}()
        new_common_neighs = Vector{Dict{Int64, Float64}}()

        for simplex_iter in eachindex(common_neighs)
            for common_neigh in keys(common_neighs[simplex_iter])
                filtration_of_this_simplex = (max(filtration_of_those_simplices[simplex_iter][1], common_neighs[simplex_iter][common_neigh]),
                max(filtration_of_those_simplices[simplex_iter][2], vertex_filtrations[common_neigh]))
                
                push!(new_filtration_of_those_simplices, filtration_of_this_simplex)
                ECC[filtration_of_this_simplex] = get(ECC, filtration_of_this_simplex, 0) + dimm
                
                neighs_of_new_simplex = Dict{Int64, Float64}()

                for v in neighs_of_vertices[common_neigh]
                    filt_val = get(common_neighs[simplex_iter], v[1], -1.0)
                    if filt_val > 0
                        neighs_of_new_simplex[v[1]] = max(v[2], filt_val)
                    end
                end
                push!(new_common_neighs, neighs_of_new_simplex)
            end
        end
        
        filtration_of_those_simplices = new_filtration_of_those_simplices
        common_neighs = new_common_neighs
        dimension += 1
        dimm = dimm * -1
    end
    
    return ECC      
end

function compute_local_contributions(point_cloud, epsilon, vertex_filtrations, max_dimension, dbg=false)
    ECP_list = Vector{Dict{Tuple{Float64, Float64}, Int64}}(undef, length(point_cloud))
    
    for i in 1:length(point_cloud)
        id_neigs_of_center_vertex = [i]
        for (j, point) in enumerate(point_cloud[i+1:end])
            if distance(point_cloud[i], point) ≤ epsilon
                push!(id_neigs_of_center_vertex, j + i)
            end
        end
        graph_i = create_local_graph(point_cloud, id_neigs_of_center_vertex[2:end], i, epsilon)
        local_ECP = compute_ECP_single_vertex(graph_i, vertex_filtrations[id_neigs_of_center_vertex], max_dimension, dbg)
        ECP_list[i] = local_ECP
    end
    
    total_ECP = Dict{Tuple{Float64, Float64}, Int64}()
    
    for single_ECP in ECP_list
        for (key, value) in single_ECP
            total_ECP[key] = get(total_ECP, key, 0) + value
        end
    end
    
    to_del = []
    for key in keys(total_ECP)
        if total_ECP[key] == 0
            push!(to_del, key)
        end
    end
    for key in to_del
        delete!(total_ECP, key)
    end
    
    return total_ECP
end
