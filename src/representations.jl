struct Unipotent end
struct Weyl end

abstract type AbstractRepresentation{T} end

struct PrincipalRepresentation{G,T,GL} <: AbstractRepresentation{T}
    # ψ(α), the value of a character ψ:K → ℂ
    # for a given generator of K
    character_value::Pair{G,T}
    borel_cd::CosetDecomposition{GL,Borel{GL}}
end

function Base.show(io::IO, ϱ::PrincipalRepresentation{GF{q}, T, GL}) where {q, T, GL}
    println(io, "Principal series representation of $GL")
    print(io, " · associated character of 𝔽_q: ", ϱ.character_value)
end

function (ϱ::PrincipalRepresentation{GF{q},T,SL₂{q}})(m::SL₂{q}) where {q,T}
    # for now only for SL₂
    u, w, D, U = bruhat(m)
    isone(w) && return ϱ(D, Diagonal) * ϱ(U, Unipotent)
    return ϱ(u, Unipotent) * ϱ(w, Weyl) * ϱ(D, Diagonal) * ϱ(U, Unipotent)
end

function (ϱ::PrincipalRepresentation{GF{q},T,SL₂{q}})(
    U::SL₂{q},
    ::Type{Unipotent},
) where {T,q}

    α, ψα = ϱ.character_value
    ϱU = fill(zero(ψα), q + 1, q + 1)

    for (i, pi) in zip(1:length(ϱ.borel_cd), right_action(U, ϱ.borel_cd))
        ϱU[pi, i] = one(ψα)
    end

    return ϱU
end

function (ϱ::PrincipalRepresentation{GF{q},T,SL₂{q}})(
    D::SL₂{q},
    ::Type{Diagonal},
) where {T,q}
    a = D[1]
    α, ψα = ϱ.character_value
    ψα_inv = inv(ψα)
    ϱD = fill(zero(ψα), q + 1, q + 1)

    j = something(findfirst(j -> α^j == a, 1:q), 0)

    perm_repr = right_action(D, ϱ.borel_cd)

    for (i, pi) in zip(1:length(ϱ.borel_cd), perm_repr)
        if ϱ.borel_cd[i] ∈ ϱ.borel_cd.trivial_coset
            ϱD[i, pi] = ψα^j
        else
            ϱD[i, pi] = ψα_inv^j
        end
    end

    return ϱD
end

function (ϱ::PrincipalRepresentation{GF{q},T})(
    w::SL₂{q},
    ::Type{Weyl},
) where {T,q}
    α, ψα = ϱ.character_value
    ϱw = fill(zero(ψα), q + 1, q + 1)

    exps = Dict(α^j => j for j = 1:q)

    perm_repr = right_action(w, ϱ.borel_cd)

    for (i, pi) in zip(1:length(ϱ.borel_cd), perm_repr)
        if ϱ.borel_cd[i] ∈ ϱ.borel_cd.trivial_coset
            ϱw[i, pi] = ψα^exps[-one(α)] # ψ(-1)
        elseif w * ϱ.borel_cd[-i] ∈ ϱ.borel_cd.trivial_coset
            ϱw[i, pi] = one(ψα)
        else
            repr = ϱ.borel_cd[i]
            # [ c    0 ][ 1 -a/c ] [ a b ]   [ 0    1 ]
            # [ 0 -1/c ][ 0    1 ] [ c d ]   [-1 -d/c ]
            c, d = repr[2], repr[4]
            # we deal with the trivial coset above
            @assert !iszero(c) #
            # we deal with the coset of w above
            @assert !iszero(d) #
            u = -d / c
            ϱw[i, pi] = ψα^exps[-inv(u)]
        end
    end
    return ϱw
end
