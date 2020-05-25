@testset "Coset decompositions" begin

    import RamanujanGraphs: CosetDecomposition, right_action

    q = 31

    SL2q = let
        a = SL₂{q}(1, 0, 1, 1)
        b = SL₂{q}(1, 1, 0, 1)

        E, sizes = RamanujanGraphs.generate_balls(
            [a, b, inv(a), inv(b)],
            radius = 15,
        )
        @test sizes[end] == order(SL₂{q})
        E
    end

    @test all(b == prod(bruhat(b)) for b in SL2q)

    coset_reps = let α = RamanujanGraphs.generator(first(SL2q)[1])
        borel_coset_representative(u::GF{q}) = SL₂{q}(0, -1, 1, -u)
        borel_coset_representative(q::Int) = one(SL₂{q})
        reps = [
            [borel_coset_representative(α^i) for i = 2:2:q-1]
            [borel_coset_representative(α^i) for i = 1:2:q-2]
            [borel_coset_representative(0 * α), borel_coset_representative(q)]
        ]
        @test length(unique(reps)) == length(reps)
        reps
    end

    @testset "perm_repr" begin

        function perm_repr(
            x::T,
            coset_representatives,
            trivial_coset,
        ) where {T<:RamanujanGraphs.AbstractGL₂}

            perm = zeros(Int, length(coset_representatives))
            inv_coset_representatives = inv.(coset_representatives)

            for (i, c) in enumerate(coset_representatives)
                a = c * x
                for (j, inv_d) in enumerate(inv_coset_representatives)
                    if a * inv_d ∈ trivial_coset
                        perm[i] = j
                        break
                    end
                end
            end
            @test isperm(perm)
            return perm
        end

        let (id, a, b) = SL2q[1:3], Borel = Borel(SL₂{q})
            perm_id = perm_repr(id, coset_reps, Borel)
            perm_a = perm_repr(a, coset_reps, Borel)
            perm_b = perm_repr(b, coset_reps, Borel)

            l = length(coset_reps)

            @test perm_repr(id, coset_reps, Borel) == 1:l
            @test !any(iszero, perm_a)
            @test !any(iszero, perm_b)
            @test length(unique(perm_a)) == length(perm_a)
            @test length(unique(perm_b)) == length(perm_b)
            @test sort(perm_a) == sort(perm_b) == 1:l

            @test perm_repr(a^2, coset_reps, Borel) ==
                  [perm_a[perm_a[i]] for i = 1:l]

            @test perm_repr(b^2, coset_reps, Borel) ==
                  [perm_b[perm_b[i]] for i = 1:l]

            @test [perm_a[i] for i in perm_repr(inv(a), coset_reps, Borel)] == 1:l

            @test [perm_b[i] for i in perm_repr(inv(b), coset_reps, Borel)] == 1:l

            @test perm_repr(a * b, coset_reps, Borel) ==
                  [perm_b[perm_a[i]] for i = 1:l]
        end
    end

    @testset "B\\SL₂{q} (left cosets)" begin

        id, a, b = SL2q[1:3]
        Bcosets = CosetDecomposition(
            coset_reps,
            inv.(coset_reps),
            Borel(eltype(SL2q)),
        )

        perm_id = right_action(id, Bcosets)
        perm_a = right_action(a, Bcosets)
        perm_b = right_action(b, Bcosets)

        # @show perm_a
        # @show perm_b

        @test a ∈ Bcosets.trivial_coset
        @test b ∉ Bcosets.trivial_coset

        @test right_action(id, Bcosets) == 1:length(Bcosets)

        if !isone(a)
            @test perm_a[end] == perm_a[end]
            @test perm_a !== 1:length(Bcosets)
        end

        let cd = Bcosets
            @test perm_b != 1:length(cd)
            @test !any(iszero, perm_a)
            @test !any(iszero, perm_b)
            @test length(unique(perm_a)) == length(perm_a)
            @test length(unique(perm_b)) == length(perm_b)
            @test sort(perm_a) == sort(perm_b) == 1:length(cd)

            @test right_action(a^2, cd) ==
                  [perm_a[perm_a[i]] for i = 1:length(cd)]

            @test right_action(b^2, cd) ==
                  [perm_b[perm_b[i]] for i = 1:length(cd)]

            @test [perm_a[i] for i in right_action(inv(a), cd)] == 1:length(cd)

            @test [perm_b[i] for i in right_action(inv(b), cd)] == 1:length(cd)

            @test right_action(a * b, cd) ==
                  [perm_b[perm_a[i]] for i = 1:length(cd)]
        end

        let cd = Bcosets,
            triv = findfirst(i -> cd[i] ∈ cd.trivial_coset, 1:length(cd))

            @test all(cd[i] * cd[-i] ∈ cd.trivial_coset for i = 1:length(cd))
            g = prod(rand(collect(cd.trivial_coset), 4))

            @test all(
                right_action(g, cd)[triv] == triv for g in cd.trivial_coset
            )

            l = length(cd)
            @test isone([
                cd[i] * cd[-j] ∈ cd.trivial_coset for i = 1:l, j = 1:l
            ])

            @test all(right_action(cd[i], cd)[triv] == i for i = 1:length(cd))

            h = prod(rand(cd.representatives, 4))
            @test sum(cd[i] * h ∈ cd.trivial_coset for i = 1:length(cd)) == 1

            perm = right_action(h, cd)
            @test all(
                cd[i] * h * cd[-j] ∈ cd.trivial_coset
                for (i, j) in enumerate(perm)
            )

            @test isperm(perm)
            perm
        end
    end

    @testset "CosetDecomposition" begin
        id, a, b = SL2q[1:3]
        Bcosets = CosetDecomposition(SL2q, Borel(SL₂{q}))

        perm_id = right_action(id, Bcosets)
        perm_a = right_action(a, Bcosets)
        perm_b = right_action(b, Bcosets)

        # @show perm_a
        # @show perm_b

        @test a ∈ Bcosets.trivial_coset
        @test b ∉ Bcosets.trivial_coset

        @test right_action(id, Bcosets) == 1:length(Bcosets)

        @test right_action(id, Bcosets) == 1:length(Bcosets)

        if !isone(a)
            @test perm_a[end] == perm_a[end]
            @test perm_a !== 1:length(Bcosets)
        end

        let cd = Bcosets
            @test perm_b != 1:length(cd)
            @test !any(iszero, perm_a)
            @test !any(iszero, perm_b)
            @test length(unique(perm_a)) == length(perm_a)
            @test length(unique(perm_b)) == length(perm_b)
            @test sort(perm_a) == sort(perm_b) == 1:length(cd)

            @test right_action(a^2, cd) ==
                  [perm_a[perm_a[i]] for i = 1:length(cd)]

            @test right_action(b^2, cd) ==
                  [perm_b[perm_b[i]] for i = 1:length(cd)]

            @test [perm_a[i] for i in right_action(inv(a), cd)] == 1:length(cd)

            @test [perm_b[i] for i in right_action(inv(b), cd)] == 1:length(cd)

            @test right_action(a * b, cd) ==
                  [perm_b[perm_a[i]] for i = 1:length(cd)]
        end

        let cd = Bcosets,
            triv = findfirst(i -> cd[i] ∈ cd.trivial_coset, 1:length(cd))

            @test all(cd[i] * cd[-i] ∈ cd.trivial_coset for i = 1:length(cd))
            g = prod(rand(collect(cd.trivial_coset), 4))

            @test all(
                right_action(g, cd)[triv] == triv for g in cd.trivial_coset
            )

            l = length(cd)
            @test isone([
                cd[i] * cd[-j] ∈ cd.trivial_coset for i = 1:l, j = 1:l
            ])

            @test all(right_action(cd[i], cd)[triv] == i for i = 1:length(cd))

            h = prod(rand(cd.representatives, 4))
            @test sum(cd[i] * h ∈ cd.trivial_coset for i = 1:length(cd)) == 1

            perm = right_action(h, cd)
            @test all(
                cd[i] * h * cd[-j] ∈ cd.trivial_coset
                for (i, j) in enumerate(perm)
            )

            @test isperm(perm)
            perm
        end
    end

    @testset "cross-checking CosetDecomposition" begin
        Borel_cd =
            CosetDecomposition(coset_reps, inv.(coset_reps), Borel(SL₂{q}))
        cd = CosetDecomposition(SL2q, Borel(SL₂{q}))

        @test all(cd[i] * cd[-i] ∈ cd.trivial_coset for i = 1:length(cd))
        g = prod(rand(SL2q, 4))

        @test sum(cd[i] * g ∈ cd.trivial_coset for i = 1:length(cd)) == 1

        l = length(cd)
        m = zeros(Int, l, l)

        for i = 1:l
            for j = 1:l
                if cd[i] * cd[-j] ∈ cd.trivial_coset
                    m[i, j] = 1
                end
            end
        end

        @test isone(m)

        perm = zeros(Int, l)

        for i = 1:l
            for j = 1:l
                if cd[i] * Borel_cd[-j] ∈ Borel_cd.trivial_coset
                    perm[i] = j
                    # break
                end
            end
        end
        @test isperm(perm)

        Borel_reps = Borel_cd.representatives
        reps = cd.representatives
        for i = 1:length(reps)
            @test reps[i] * inv(Borel_reps[perm[i]]) ∈ Borel_cd.trivial_coset
        end
    end
end
