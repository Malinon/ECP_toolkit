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

@testset "vectorisation test exponential" begin
    contributions = Dict{Tuple{Int64, Int64, Int64}, Int64}()
    SIZE = 2
    for i in 1:SIZE
        for j in 1:SIZE
            for k in 1:SIZE
                contributions[(i, j, k)] = (2^i) * (3^j) * (5^k)
            end
        end
    end
    vectorisation = Ecp.vectorize_ecp(contributions, (SIZE - 1, SIZE - 1, SIZE - 1))
    @test vectorisation[1, 1, 1] == 30
    @test vectorisation[2, 1, 1] == 90
    @test vectorisation[1, 2, 1] == 120
    @test vectorisation[1, 1, 2] == 180
    @test vectorisation[2, 2, 1] == 30 + 60 + 90 + 180
    @test vectorisation[2, 2, 2] == vectorisation[2, 2, 1] * 6
    @test vectorisation[2, 1, 2] == vectorisation[2, 2, 2] / 4
    @test vectorisation[1, 2, 2] == vectorisation[2, 2, 2] / 3
end

@testset "cubic ECP test" begin
    @testset "solid cube test" begin
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
    end

@testset "sphere to cube test" begin
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
end


@testset "solid square test" begin
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
end


@testset "ring to square test" begin
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

@testset "diagonal snake test 3d" begin
    snake_line = Array{Tuple{Float64, Float64}}(undef, 6, 6, 6)
    for i in eachindex(snake_line)
        snake_line[i] = (7.0, 8.0)
    end
    for i in 1:6
        snake_line[i,i,i] = (i,i)
    end
    ecp_contributions = Ecp.compute_contributions_3d(snake_line)
    @test ecp_contributions[(1.0, 1.0)] == 1
    for i in 2:6
        @test ecp_contributions[(i, i)] == 0
    end
end

@testset "diagonal snake test 2d" begin
    snake_line = Array{Tuple{Float64, Float64, Float64}}(undef, 6, 6)
    for i in eachindex(snake_line)
        snake_line[i] = (7.0, 8.0, 16.0)
    end
    for i in 1:6
        snake_line[i,i] = (i,i,i)
    end
    ecp_contributions = Ecp.compute_contributions_2d(snake_line)
    @test ecp_contributions[(1.0, 1.0, 1.0)] == 1
    for i in 2:6
        @test ecp_contributions[(i, i, i)] == 0
    end
end

@testset "two different holes test 2d" begin
    two_holes = Array{Tuple{Float64, Float64}}(undef, 4, 3)
    for i in eachindex(two_holes)
        two_holes[i] = (1.0, 1.0)
    end
    two_holes[2, 2] = (1.0, 2.0)
    two_holes[3, 2] = (2.0, 1.0)
    ecp_contributions = Ecp.compute_contributions_2d(two_holes)
    @test ecp_contributions[(1.0, 1.0)] == -1
    @test ecp_contributions[(1.0, 2.0)] == 1
    @test ecp_contributions[(2.0, 1.0)] == 1
    vectorisation = Ecp.vectorize_ecp(ecp_contributions, (2,2))
    for i in 1:2
        for j in 1:2
            @test vectorisation[1, 1] == -1
        end
    end
    @test vectorisation[3, 3] == 1
    @test vectorisation[1, 3] == 0
    @test vectorisation[3, 1] == 0
end

@testset "two different holes test 3d" begin
    two_holes = Array{Tuple{Float64, Float64}}(undef, 4, 3, 3)
    for i in eachindex(two_holes)
        two_holes[i] = (1.0, 1.0)
    end
    two_holes[2, 2, 2] = (1.0, 2.0)
    two_holes[3, 2, 2] = (2.0, 1.0)
    ecp_contributions = Ecp.compute_contributions_3d(two_holes)
    @test ecp_contributions[(1.0, 1.0)] == 3
    @test ecp_contributions[(1.0, 2.0)] == -1
    @test ecp_contributions[(2.0, 1.0)] == -1
    vectorisation = Ecp.vectorize_ecp(ecp_contributions, (2,2))
    for i in 1:2
        for j in 1:2
            @test vectorisation[i, j] == 3
        end
    end
    @test vectorisation[3, 3] == 1
    @test vectorisation[1, 3] == 2
    @test vectorisation[3, 1] == 2
end

end
