struct Unipotent end
struct Weyl end

abstract type AbstractRepresentation{T} end

struct PrincipalRepresentation{G,T,GL<:AbstractGL₂} <: AbstractRepresentation{T}
    # ψ(α), the value of a character ψ:K → ℂ
    # for a given generator of K
    character::Dict{G,T}
    borel_cd::CosetDecomposition{GL,Borel{GL}}

    function PrincipalRepresentation(α_ψα::Pair{G, T},  boreldec::CosetDecomposition{GL, Borel{GL}}) where
        {G, T, q, GL<:AbstractGL₂{q}}
        α, ψα = α_ψα
        new{G, T, GL}(Dict(α^j=>ψα^j for j in 1:q), boreldec)
    end
end

function Base.show(io::IO, ϱ::PrincipalRepresentation{GF{q}, T, GL}) where {q, T, GL}
    println(io, "Principal series representation of $GL")
    print(io, " · associated character of 𝔽_q: ", ϱ.character)
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

    𝟏 = one(last(first(ϱ.character)))
    ϱU = fill(zero(𝟏), q + 1, q + 1)

    for (i, pi) in zip(1:length(ϱ.borel_cd), right_action(U, ϱ.borel_cd))
        ϱU[pi, i] = 𝟏
    end

    return ϱU
end

function (ϱ::PrincipalRepresentation{GF{q},T,SL₂{q}})(
    D::SL₂{q},
    ::Type{Diagonal},
) where {T,q}
    a = D[1]
    ψa = ϱ.character[a]
    ψa_inv = inv(ψa)
    ϱD = fill(zero(ψa), q + 1, q + 1)

    perm_repr = right_action(D, ϱ.borel_cd)

    for (i, pi) in zip(1:length(ϱ.borel_cd), perm_repr)
        if ϱ.borel_cd[i] ∈ ϱ.borel_cd.trivial_coset
            ϱD[i, pi] = ψa
        else
            ϱD[i, pi] = ψa_inv
        end
    end

    return ϱD
end

function (ϱ::PrincipalRepresentation{GF{q},T})(
        w::SL₂{q}, ::Type{Weyl}) where {T,q}

    a, ψa = first(ϱ.character)
    ϱw = fill(zero(ψa), q + 1, q + 1)

    perm_repr = right_action(w, ϱ.borel_cd)

    for (i, pi) in zip(1:length(ϱ.borel_cd), perm_repr)
        if ϱ.borel_cd[i] ∈ ϱ.borel_cd.trivial_coset
            ϱw[i, pi] = ϱ.character[-one(a)] # ψ(-1)
        elseif w * ϱ.borel_cd[-i] ∈ ϱ.borel_cd.trivial_coset
            ϱw[i, pi] = one(ϱw[i, pi])
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
            ϱw[i, pi] = ϱ.character[-inv(u)]
        end
    end
    return ϱw
end
