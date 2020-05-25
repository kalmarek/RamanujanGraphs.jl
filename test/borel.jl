@testset "Borel subgroup" begin
    q = 31
    @test Borel{SL₂{q}}() isa Borel
    @test Borel(SL₂{q}) isa Borel

    B = Borel(SL₂{q})

    @test length(B) == q * (q - 1)
    @test first(B) == one(SL₂{q})
    @test collect(B) isa Vector{<:SL₂{q}}

    @test length(Set(collect(B))) == length(B)
    @test all(RamanujanGraphs.isupper, B)

    @test all(b == prod(bruhat(b)) for b in B)

    q = 31

    SL2q = let
        a = SL₂{q}(1, 0, 1, 1)
        b = SL₂{q}(1, 1, 0, 1)

        E, sizes = RamanujanGraphs.generate_balls(
            [a, b, inv(a), inv(b)],
            radius = 15,
        )
        @assert sizes[end] == order(SL₂{q})
        E
    end

    @test all(g == prod(bruhat(g)) for g in SL2q)
end
