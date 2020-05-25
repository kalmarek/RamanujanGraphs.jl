@testset "Principal Series" begin
    q = 5

    SL2q = let
        a = SL₂{q}(1, 0, 1, 1)
        b = SL₂{q}(1, 1, 0, 1)

        E, sizes = RamanujanGraphs.generate_balls([a,b, inv(a), inv(b)], radius=6);
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



    Borel_cd = RamanujanGraphs.CosetDecomposition(coset_reps, inv.(coset_reps), Borel(SL₂{q}))
    α = RamanujanGraphs.generator(one(GF{q}))
    ϱ = PrincipalRepr(α=>1im, Borel_cd)
    @info ϱ

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

        @test all(ϱ(g*h) == ϱ(g)*ϱ(h) for g in SL2q for h in SL2q)
    end

    for val in [1, im, -1, -im]
        ϱ = PrincipalRepr(α=>val,
            RamanujanGraphs.CosetDecomposition(SL2q, Borel(SL₂{q})))

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

            @test all(ϱ(g*h) == ϱ(g)*ϱ(h) for g in SL2q for h in SL2q)
        end
    end
end

@testset "Discrete Series" begin

    q = 7

    SL2q = let
        a,b = let
            a = SL₂{q}(1, 0, 1, 1)
            b = SL₂{q}(1, 1, 0, 1)
            a, b
        end

        E, sizes = RamanujanGraphs.generate_balls([a,b, inv(a), inv(b)], radius=7);
        @assert sizes[end] == order(SL₂{q})
        E
    end

    ζ, ξ = let k = 2
        ζ, ξ = exp(2π*im/q), exp((2π*im)/((q+1)÷2))^(k*(q-1)÷4)
        # ζ, ξ = exp(2π*im/q), 1.0im
    end

    α = RamanujanGraphs.generator(GF{q}(0))
    β = generator(QuadraticExt(α))

    h = RamanujanGraphs.DiscreteRepr(one(GF{q})=>ζ, β=>ξ)
    @info h


    @testset "DiscreteRepr" begin

        tol = 1e-10

        @test all(isapprox.(0.0,
            [h.indecomposable[β^i] - h.indecomposable[β]^i for i in 1:q^2-1], atol=tol))
        @test all(isapprox.(0.0,
            [h.decomposable[i*α] - h.decomposable[α]^i for i in 1:q], atol=tol))

        Id = Matrix{Float64}(I, RamanujanGraphs.degree(h), RamanujanGraphs.degree(h));
        @test all(isapprox.(h(SL2q[1], RamanujanGraphs.Unipotent), Id, atol=tol))
        @test all(isapprox.(h(SL2q[1], Diagonal), Id, atol=tol))

        @test all(isapprox.(h(SL2q[1]), Id, atol=tol))

        w = h(SL2q[3], RamanujanGraphs.Weyl)

        @test all(isapprox.(w^4, Id, atol=tol))
        @test all(isapprox.(h(SL2q[2]^2), h(SL2q[2])^2, atol=tol))
        @test all(isapprox.(h(SL2q[3]^2), h(SL2q[3])^2, atol=tol))

        @test all(isapprox.(h(SL2q[3]*SL2q[2]), h(SL2q[3])*h(SL2q[2]), atol=tol))

        @test all(all(isapprox.(h(a*b), h(a)*h(b), atol=tol))
            for a in SL2q[rand(1:length(SL2q), 50)]
                for b in SL2q[rand(1:length(SL2q), 50)])
    end
end
