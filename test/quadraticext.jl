@testset "Quadratic Extension of GF{q}" begin
    q = 11
    k = generator(GF{q}(1))
    @test !issqrt(k)

    @test QuadraticExt{k}(one(k), zero(k)) isa QuadraticExt
    @test QuadraticExt(k) isa QuadraticExt
    L = QuadraticExt{k}(one(k), zero(k))
    @test RamanujanGraphs.oneim(L) isa QuadraticExt
    @test reim(RamanujanGraphs.oneim(L)) == (0, 1)
    imL = RamanujanGraphs.oneim(L)
    @test imL == QuadraticExt{k}(zero(k), one(k))

    @testset "arithmetic" begin
        @test inv(L) isa QuadraticExt
        @test -L isa QuadraticExt

        for (a, b) in ((L, L), (L, k))
            @test a + b isa QuadraticExt
            @test b + a isa QuadraticExt
            @test a * b isa QuadraticExt
            @test b * a isa QuadraticExt

            @test a / b isa QuadraticExt
            @test b / a isa QuadraticExt
        end
    end

    C = RamanujanGraphs.Units(L)
    @test length(C) == q + 1
    @test eltype(C) == QuadraticExt{k,typeof(k)}
    @test collect(C) isa Vector{<:QuadraticExt}
    cC = collect(C)

    @test last(cC) == QuadraticExt{k,GF{q}}((-1, 0))
    @test all(isone, norm.(C))
    @test all(
        ==(-RamanujanGraphs.sqroot(L)),
        norm.(RamanujanGraphs.oneim(L) .* C),
    )

    l = generator(L)
    @test length(unique!([l^i for i = 1:q^2-1])) == q^2 - 1
    # frobenius map
    @test all(isreal, (l^i)^(q + 1) for i = 1:q^2-1)
    # trace map
    @test all(isreal, (l^i) + (l^i)^q for i = 1:q^2-1)
end
