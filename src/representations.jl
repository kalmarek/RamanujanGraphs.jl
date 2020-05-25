struct Unipotent end
struct Weyl end

abstract type AbstractRepresentation{T} end

####
# Principal series

struct PrincipalRepr{K,T,GL<:AbstractGL₂} <: AbstractRepresentation{T}
    character::Dict{K,T}
    borel_cd::CosetDecomposition{GL,Borel{GL}}

    function PrincipalRepr(
        α_ψα::Pair{G,T},
        boreldec::CosetDecomposition{GL,Borel{GL}},
    ) where {G,T,q,GL<:AbstractGL₂{q}}
        @assert length(boreldec) == q+1
        α, ψα = α_ψα
        new{G,T,GL}(Dict(α^j => ψα^j for j = 1:q), boreldec)
    end
end

degree(ϱ::PrincipalRepr{GF{q},T, SL₂{q}}) where {q, T} = q+1

function Base.show(io::IO, ϱ::PrincipalRepr{GF{q}, T, GL}) where {q, T, GL}
    α = first(first(ϱ.character))
    println(io, "Principal series representation of $GL")
    print(io, "\tcharacter of 𝔽$(subscriptify(q))ˣ: ", α, " → ", ϱ.character[α])
end

function (ϱ::PrincipalRepr{GF{q},T,SL₂{q}})(m::SL₂{q}) where {q,T}
    # for now only for SL₂
    u, w, D, U = bruhat(m)
    isone(w) && return ϱ(D, Diagonal) * ϱ(U, Unipotent)
    return ϱ(u, Unipotent) * ϱ(w, Weyl) * ϱ(D, Diagonal) * ϱ(U, Unipotent)
end

function (ϱ::PrincipalRepr{GF{q},T,SL₂{q}})(U::SL₂{q},::Type{Unipotent}) where {q, T}

    𝟙 = one(last(first(ϱ.character)))
    ϱU = fill(zero(𝟙), degree(ϱ), degree(ϱ))

    for (i, pi) in zip(1:degree(ϱ), right_action(U, ϱ.borel_cd))
        ϱU[i, pi] = 𝟙
    end

    return ϱU
end

function (ϱ::PrincipalRepr{GF{q},T,SL₂{q}})(D::SL₂{q},::Type{Diagonal}) where {q, T}
    a = D[1]
    ψa = ϱ.character[a]
    ψa_inv = inv(ψa)
    ϱD = fill(zero(ψa), degree(ϱ), degree(ϱ))

    perm_repr = right_action(D, ϱ.borel_cd)

    for (i, pi) in zip(1:degree(ϱ), perm_repr)
        if ϱ.borel_cd[i] ∈ ϱ.borel_cd.trivial_coset
            ϱD[i, pi] = ψa
        else
            ϱD[i, pi] = ψa_inv
        end
    end

    return ϱD
end

function (ϱ::PrincipalRepr{GF{q},T})(w::SL₂{q}, ::Type{Weyl}) where {q, T}

    𝟙 = one(last(first(ϱ.character)))
    ϱw = fill(zero(𝟙), degree(ϱ), degree(ϱ))

    perm_repr = right_action(w, ϱ.borel_cd)

    for (i, pi) in zip(1:degree(ϱ), perm_repr)
        if ϱ.borel_cd[i] ∈ ϱ.borel_cd.trivial_coset
            ϱw[i, pi] = ϱ.character[GF{q}(-1)] # ψ(-1)
        elseif w * ϱ.borel_cd[-i] ∈ ϱ.borel_cd.trivial_coset
            ϱw[i, pi] = 𝟙 # ψ(1)
        else
            repr = ϱ.borel_cd[i]
            # [ c    0 ][ 1 -a/c ][ a b ] =  [ 0    1 ]
            # [ 0 -1/c ][ 0    1 ][ c d ]    [-1 -d/c ]
            c, d = repr[2], repr[4]
            # we deal with the trivial coset above
            @assert !iszero(c)
            # we deal with the coset of w above
            @assert !iszero(d)
            u = -d / c
            ϱw[i, pi] = ϱ.character[inv(u)]
        end
    end
    return ϱw
end

####
# Discrete series

struct DiscreteRepr{K,T,GL,ε} <: AbstractRepresentation{T}
    decomposable::Dict{K,T}
    indecomposable::Dict{QuadraticExt{ε, K},T}
    basis::Dict{K,Int}
    function DiscreteRepr(
        α_χα::Pair{GF{q},T},
        β_νβ::Pair{QuadraticExt{ε, GF{q}},T},
        basis = Dict(generator(first(α_χα))^i => i for i = 1:q-1),
    ) where {q,T,ε}
        @assert length(basis) == q - 1 "Basis should be Kˣ"
        @assert !isone(last(α_χα)) "α_χα should be a non-trivial character of K⁺"
        @assert isone(first(α_χα)) "α_χα should be a non-trivial character of K⁺"
        decomposable = Dict(i * first(α_χα) => last(α_χα)^i for i = 0:q-1)
        indecomposable = Dict(
            generator(first(β_νβ))^i => last(β_νβ)^i for i = 1:q^2-1)

        @assert any(!isone, (indecomposable[u] for u in RamanujanGraphs.Units(first(β_νβ)))) "β_νβ should be a non-trivial character of Lˣ"
        return new{GF{q},T,SL₂{q}, ε}(decomposable, indecomposable, basis)
    end
end

degree(ϱ::DiscreteRepr{GF{q},T,SL₂{q}}) where {q,T} = q - 1
basis(ϱ::DiscreteRepr) = ϱ.basis

function Base.show(io::IO, ϱ::DiscreteRepr{GF{q},T,GL}) where {q,T,GL}
    χ = ϱ.decomposable
    ν = ϱ.indecomposable

    q_idx = subscriptify(q)

    β = generator(first(first(ν)))

    println(io, "Discrete series representation of $GL")
    println(io, "\tcharacter of 𝔽$(q_idx): ", GF{q}(1), " → ", χ[GF{q}(1)])
    println(io, "\tcharacter of 𝔽$(q_idx)(√$(sqroot(β))ˣ: ", β, " → ", ν[β])
end

function (ϱ::DiscreteRepr{GF{q},T,SL₂{q}})(m::SL₂{q}) where {q,T}
    # for now only for SL₂
    B = Bruhat(m)
    isupper(m) && return ϱ(B.D, Diagonal) * ϱ(B.U, Unipotent)
    return ϱ(B.u, Unipotent) * ϱ(B.w, Weyl) * ϱ(B.D, Diagonal) * ϱ(B.U, Unipotent)
end

function (ϱ::DiscreteRepr{GF{q},T,SL₂{q}})(U::SL₂{q},::Type{Unipotent}) where {q,T}
    b = U[3]
    χ = ϱ.decomposable
    ϱU = fill(zero(last(first(χ))), degree(ϱ), degree(ϱ))

    for (x, i) in basis(ϱ)
        ϱU[i, i] = χ[x*b]
    end

    return ϱU
end

function (ϱ::DiscreteRepr{GF{q},T,SL₂{q}, ε})(D::SL₂{q}, ::Type{Diagonal}) where {q,T,ε}
    a² = D[1]^2
    d = D[4]

    χ = ϱ.decomposable
    ν = ϱ.indecomposable

    ϱD = fill(zero(last(first(χ))), degree(ϱ), degree(ϱ))

    for (x, i) in basis(ϱ)
        ϱD[i, basis(ϱ)[a²*x]] = ν[QuadraticExt{ε}(d, zero(d))]
    end

    return ϱD
end

function _j(z::GF{q}, χ, ν) where q
    tmp = first(first(ν))
    w = _elt_of_norm(tmp, z)

    res = zero(last(first(ν)))
    for u in Units(tmp)
        t = w*u
        res += χ[real(t+t^q)]*ν[t]
    end
    return res/oftype(res, q)

    # L = keys(ν)
    # return sum(χ[real(t+t^q)]*ν[t] for t in L if norm(t) == z)//q
end

function (ϱ::DiscreteRepr{GF{q},T,SL₂{q}, ε})(w::SL₂{q}, ::Type{Weyl}) where {q, T, ε}

    χ = ϱ.decomposable
    ν = ϱ.indecomposable

    ϱw = fill(zero(last(first(χ))), degree(ϱ), degree(ϱ))

    for (y, j) in basis(ϱ)
        νy¯¹ = ν[QuadraticExt{ε}(inv(y), zero(y))]
        for (x, i) in basis(ϱ)
            ϱw[i, j] += -νy¯¹*_j(x*y, χ, ν)
        end
    end

    return ϱw
end
