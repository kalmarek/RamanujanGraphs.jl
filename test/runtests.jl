using Test

using LinearAlgebra
using RamanujanGraphs
using LightGraphs

@testset "GL₂" begin
    @testset "PGL₂" begin
        @test_throws AssertionError PGL₂{5}(-4, 2 , 6, 2)
        m = PGL₂{5}(-4, 2 , 6, 4)
        @test m == PGL₂{5}(1, 2, 1, 4)
        @test det(m) == 2
        @test invmod(det(m), 5) == 3
        @test m*m == PGL₂{5}(3, 0, 0, 3)
        @test m^2 == PGL₂{5}(1, 0, 0, 1)
        @test inv(m) * m == PGL₂{5}(1, 0, 0, 1)

        m = PGL₂{5}([0 3; 1 4])
        @test m[1] == 0
        @test m[2] == 1
        @test m[3] == 3
        @test m[4] == 4

        @test m[1,1] == 0
        @test m[1,2] == 3
        @test m[2,1] == 1
        @test m[2,2] == 4
    end

    @testset "PSL₂" begin
        @test_throws ArgumentError PSL₂{5}(-4, 2 , 6, 4)
        m = PSL₂{5}(-4, 2 , 6, 3)
        @test det(m) == 1
        @test m^3 == PSL₂{5}([3 0; 0 3])
        @test m^3 == PSL₂{5}([1 0; 0 1])
        @test PGL₂{5}(m) isa PGL₂{5}
        @test m == PSL₂{5}(PGL₂{5}(m))
        @test inv(m) == m^2
    end
end

@testset "quadruples" begin
    @test length(RamanujanGraphs.quadruples_4k_plus1(5)) == 6
    @test length(RamanujanGraphs.quadruples_4k_plus1(13)) == 14
    @test length(RamanujanGraphs.quadruples_4k_plus1(17)) == 18

    @test length(RamanujanGraphs.quadruples_4k_plus3(7)) == 8
    @test length(RamanujanGraphs.quadruples_4k_plus3(11)) == 12
    @test length(RamanujanGraphs.quadruples_4k_plus3(19)) == 20
end

@testset "lpsgenerators" begin

    @test lps_generators(3, 5) isa Vector{PGL₂{5}}
    @test lps_generators(7, 5) isa Vector{PGL₂{5}}
    @test lps_generators(11, 5) isa Vector{PSL₂{5}}

    @test length(lps_generators(3, 5)) == 4
    @test length(lps_generators(7, 5)) == 8
    @test length(lps_generators(11, 5)) == 12


    for p in [3,7,11]
        S = lps_generators(p, 5);
        E, sizes = RamanujanGraphs.generate_balls(S, radius=5);
        @test all(isequal(p+1), (length(unique(g*s for s in S)) for g in E))
    end
    # S = lps_generators(7, 5);
    # E, sizes = generate_balls(S, radius=4);
    # @test all(isequal(8), (length(unique(g*s for s in S)) for g in E))
    #
    # S = lps_generators(11, 5);
    # E, sizes = generate_balls(S, radius=4);
    # @test all(isequal(12), (length(unique(g*s for s in S)) for g in E))
end


@testset "Cayley graphs" begin
    let (p, q) = (3, 5)
        S = lps_generators(p, q)
        c, vertices, vlabels = cayley_graph(S, radius=6)
        @test length(vertices) == RamanujanGraphs.order(eltype(vertices))
        @test all(isequal(p+1), degree(c))
        @test all(v in keys(vlabels) for v in vertices)
        @test all(vlabels[vertices[i]] == i for (g,i) in vlabels)
    end
end
