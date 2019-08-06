using Test

using LinearAlgebra
using RamanujanGraphs
using LightGraphs

@testset "GL₂" begin
    @testset "PGL₂" begin
        @test_throws ArgumentError PGL₂{5}(-4, 2 , 6, 2)
        m = PGL₂{5}(-4, 2 , 6, 4)
        @test m == PGL₂{5}(1, 2, 1, 4)
        @test det(m) == 2
        @test det(m) isa RamanujanGraphs.IntMod{5}
        @test inv(det(m)) == 3
        @test m*m == PGL₂{5}(3, 0, 0, 3)
        @test m^2 == PGL₂{5}(1, 0, 0, 1)
        @test inv(m) * m == PGL₂{5}(1, 0, 0, 1)

        # testing correct indexing
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

    for p in [3,5,7,11,13,17,19,23]
        @test length(RamanujanGraphs.quadruples(p)) == p+1
    end
end

@testset "lpsgenerators" begin

    @test lps_generators(13, 5) isa Vector{PGL₂{5}}
    @test lps_generators(17, 13) isa Vector{PSL₂{13}}

    for p in [5, 13, 17, 29, 37]
        for q in [13, 17, 29, 37]
            if p ≠ q
                @test length(lps_generators(p, q)) == p+1
            end
        end
    end

    for p in [5, 13, 17, 29, 37]
        for q in [13, 17]
            if p ≠ q
                S = lps_generators(p, q);
                E, sizes = RamanujanGraphs.generate_balls(S, radius=2p)
                @test sizes[end] == RamanujanGraphs.order(eltype(S))
                @test all(isequal(p+1), (length(unique(g*s for s in S)) for g in E))
            end
        end
    end
end


@testset "Cayley" begin
    let (p, q) = (13, 5)
        S = lps_generators(p, q)

        c, vs, vls, els = cayley_graph(S, 6)
        @test length(vs) == RamanujanGraphs.order(eltype(vs))
        @test all(isequal(p+1), LightGraphs.degree(c))
        @test all(v in keys(vls) for v in vs)
        @test all(vls[vs[i]] == i for i in eachindex(vs))

        c, vs, vls, els = cayley_graph(S)
        @test length(vs) == RamanujanGraphs.order(eltype(vs))
        @test all(isequal(p+1), LightGraphs.degree(c))
        @test all(v in keys(vls) for v in vs)
        @test all(vls[vs[i]] == i for i in eachindex(vs))
    end

    let (p, q) = (37, 29)
        S = lps_generators(p, q)

        c, vs, vls, els = cayley_graph(S, 6)
        @test length(vs) == RamanujanGraphs.order(eltype(vs))
        @test all(isequal(p+1), LightGraphs.degree(c))
        @test all(v in keys(vls) for v in vs)
        @test all(vls[vs[i]] == i for i in eachindex(vs))

        c, vs, vls, els = cayley_graph(S)
        @test length(vs) == RamanujanGraphs.order(eltype(vs))
        @test all(isequal(p+1), LightGraphs.degree(c))
        @test all(v in keys(vls) for v in vs)
        @test all(vls[vs[i]] == i for i in eachindex(vs))
    end

    c, vs, vls, els = LPS(13, 5)
    @test length(vs) == RamanujanGraphs.order(eltype(vs))
    @test nv(c) == length(vs)
end
