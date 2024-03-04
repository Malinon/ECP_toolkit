import Ecp

using Test

@testset "vectorisation test" begin
    contributions = Dict{Tuple{Int64, Int64, Int64}, Int64}()
    SIZE = 2
    for i in 1:SIZE
        for j in 1:SIZE
            for k in 1:SIZE
                contributions[(i, j, k)] = -1
            end
        end
    end
    vectorisation = Ecp.vectorize_ecp(contributions, (SIZE - 1, SIZE - 1, SIZE - 1))
    for i in 1:SIZE
        for j in 1:SIZE
            for k in 1:SIZE
                @test vectorisation[i, j, k] == - i * j * k
            end
        end
    end
end

@testset "cubic ECP test" begin
    solid_cube = Array{Tuple{Float64, Float64}}(undef, 5, 5, 5)
    for i in eachindex(solid_cube)
        solid_cube[i] = (2.0, 3.0)
    end
    ecp_contributions = Ecp.compute_contributions_3d(solid_cube)
    @test ecp_contributions[(2.0, 3.0)] == 1
    vectorisation = Ecp.vectorize_ecp(ecp_contributions, (4,4))
    for i in 1:5
        for j in 1:5
            @test vectorisation[i, j] == 1
        end
    end


    sphere_to_cube = Array{Tuple{Float64, Float64}}(undef, 3, 3, 3)
    for i in eachindex(sphere_to_cube)
        sphere_to_cube[i] = (2.0, 3.0)
    end
    sphere_to_cube[2, 2, 2] = (6.0, 4.0)
    ecp_contributions = Ecp.compute_contributions_3d(sphere_to_cube)
    @test ecp_contributions[(2.0, 3.0)] == 2
    @test ecp_contributions[(6.0, 4.0)] == -1
    vectorisation = Ecp.vectorize_ecp(ecp_contributions, (4,4))
    for i in 1:5
        for j in 1:5
            if i == 5 & j == 5
                @test vectorisation[i, j] == 1
            else
                @test vectorisation[i, j] == 2
            end
        end
    end

    solid_square = Array{Tuple{Float64, Float64, Float64}}(undef, 5, 5)
    for i in eachindex(solid_square)
        solid_square[i] = (2.0, 3.0, 1.0)
    end
    ecp_contributions = Ecp.compute_contributions_2d(solid_square)
    @test ecp_contributions[(2.0, 3.0, 1.0)] == 1
    vectorisation = Ecp.vectorize_ecp(ecp_contributions, (4,4,4))
    for i in 1:5
        for j in 1:5
            for k in 1:5
                @test vectorisation[i, j, k] == 1
            end
        end
    end


    ring_to_square = Array{Tuple{Float64, Float64}}(undef, 3, 3)
    for i in eachindex(ring_to_square)
        ring_to_square[i] = (2.0, 3.0)
    end
    ring_to_square[2, 2] = (5.0, 6.0)
    ecp_contributions = Ecp.compute_contributions_2d(ring_to_square)
    @test ecp_contributions[(2.0, 3.0)] == 0
    @test ecp_contributions[(5.0, 6.0)] == 1
    vectorisation = Ecp.vectorize_ecp(ecp_contributions, (4,4))
    for i in 1:5
        for j in 1:5
            if i == 5 && j == 5
                @test vectorisation[i, j] == 1
            else
                @test vectorisation[i, j] == 0
            end
        end
    end

end
