import Ecp
using Test


@testset "vr ECP test" begin
    @testset "square test" begin
        #square test
        square_poins = [(0,0), (1,0), (1,1), (0,1)]
        time_filtraion = [1.0, 2.0 ,3.0 , 4.0]

        ecp_contributions = Ecp.compute_local_contributions(square_poins, 2.0, time_filtraion, 2, false)
        println(ecp_contributions)
        @test ecp_contributions[(0.0, 1.0)] == 1
        @test ecp_contributions[(0.0, 2.0)] == 1
        @test ecp_contributions[(0.0, 3.0)] == 1
        @test ecp_contributions[(0.0, 4.0)] == 1
        @test ecp_contributions[(1.0, 2.0)] == -1
        @test ecp_contributions[(1.0, 3.0)] == -1
        @test ecp_contributions[(1.0, 4.0)] == -2
        @test ecp_contributions[(sqrt(2.0), 4.0)] == 2
    end
    @testset "line test" begin
        # line test
        line_poins = [(0,0), (1,0), (2,0)]
        time_filtraion = [1.0, 2.0 ,3.0]

        ecp_contributions = Ecp.compute_local_contributions(line_poins, 2.0, time_filtraion, 2, false)
        @test ecp_contributions[(0.0, 1.0)] == 1
        @test ecp_contributions[(0.0, 2.0)] == 1
        @test ecp_contributions[(0.0, 3.0)] == 1
        @test ecp_contributions[(1.0, 2.0)] == -1
        @test ecp_contributions[(1.0, 3.0)] == -1
    end
end
