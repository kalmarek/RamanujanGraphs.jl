@testset "Principal Representations" begin
    q = 5

    SL2q = let
        a = SL₂{q}(1, 0, 1, 1)
        b = SL₂{q}(1, 1, 0, 1)

        E, sizes = RamanujanGraphs.generate_balls([a,b, inv(a), inv(b)], radius=20);
        @assert sizes[end] == order(SL₂{q})
        E
    end

    # characteristic(::Type{<:AbstractGL₂{q}}) where q = q

    coset_reps = let α = RamanujanGraphs.generator(first(SL2q)[1])
        borel_coset_representative(u::GF{q}) = SL₂{q}(0, -1, 1, -u)
        borel_coset_representative(q::Int) = one(SL₂{q})
        reps = [
            [borel_coset_representative(α^i) for i in 2:2:q-1];
            [borel_coset_representative(α^i) for i in 1:2:q-2];
            [borel_coset_representative(0*α), borel_coset_representative(q)];
        ]
        @assert length(unique(reps)) == length(reps)
        reps
    end


    Borel_cd = CosetDecomposition(coset_reps, inv.(coset_reps), Borel(SL₂{q}))
    α = RamanujanGraphs.generator(one(GF{q}))
    ϱ = PrincipalRepresentation(α=>1im, Borel_cd)

    @testset "specific basis" begin

        @test isone(det(ϱ(SL2q[2], RamanujanGraphs.Unipotent)))

        @test isone(ϱ(SL2q[1]))
        @test !isone(ϱ(SL2q[2]))
        @test !isone(ϱ(SL2q[3]))

        @test isone(det(ϱ(SL₂{q}(3^2, 0, 0, invmod(3^2, q)), Diagonal)))

        @test isone(det(ϱ(SL₂{q}(0, -1, 1, 0), RamanujanGraphs.Weyl)))

        @test isone(ϱ(SL2q[1]))
        @test ϱ(SL2q[2]) == ϱ(SL2q[2], RamanujanGraphs.Unipotent)
        @test isone(det(ϱ(SL2q[3])))
        @test isone(ϱ(SL2q[2])*ϱ(inv(SL2q[2])))
        @test isone(ϱ(SL2q[3])*ϱ(inv(SL2q[3])))

        @test ϱ(SL2q[2]*SL2q[3]) == ϱ(SL2q[2])*ϱ(SL2q[3])

        for g in SL2q
            @test all(ϱ(g*h) == ϱ(g)*ϱ(h) for h in SL2q)
        end
    end

    for val in [1, im, -1, -im]
        ϱ = PrincipalRepresentation(α=>val,
            CosetDecomposition(SL2q, Borel(SL₂{q})))

        @testset "arbirtrary basis: $val" begin

            @test isone(det(ϱ(SL2q[2], RamanujanGraphs.Unipotent)))

            @test isone(ϱ(SL2q[1]))
            @test !isone(ϱ(SL2q[2]))
            @test !isone(ϱ(SL2q[3]))

            @test isone(det(ϱ(SL₂{q}(3^2, 0, 0, invmod(3^2, q)), Diagonal)))

            @test isone(det(ϱ(SL₂{q}(0, -1, 1, 0), RamanujanGraphs.Weyl)))

            @test isone(ϱ(SL2q[1]))
            @test ϱ(SL2q[2]) == ϱ(SL2q[2], RamanujanGraphs.Unipotent)
            @test isone(det(ϱ(SL2q[3])))
            @test isone(ϱ(SL2q[2])*ϱ(inv(SL2q[2])))
            @test isone(ϱ(SL2q[3])*ϱ(inv(SL2q[3])))

            @test ϱ(SL2q[2]*SL2q[3]) == ϱ(SL2q[2])*ϱ(SL2q[3])

            for g in SL2q
                @test all(ϱ(g*h) == ϱ(g)*ϱ(h) for h in SL2q)
            end
        end
    end
end
