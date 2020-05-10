abstract type AbstractGL₂{q} <: AbstractMatrix{GF} end

Base.size(::AbstractGL₂) = (2,2)
Base.length(::AbstractGL₂) = 4

Base.one(::Type{T}) where T <: AbstractGL₂ = T(1, 0, 0, 1)
Base.one(::T) where T <: AbstractGL₂ = one(T)
Base.similar(::T) where T <: AbstractGL₂ = one(T)

Base.IndexStyle(::AbstractGL₂) = IndexLinear()

function Base.getindex(m::AbstractGL₂, i::Integer)
    i == 1 && return m.a
    i == 2 && return m.c # julia assumes column wise storage for IndexLinear
    i == 3 && return m.b
    i == 4 && return m.d
    throw(BoundsError(m, i))
end

function Base.setindex!(m::AbstractGL₂, v, i::Integer)
    @boundscheck 1 ≤ i ≤ 4 || throw(BoundsError(m, i))
    if i == 1
        m.a = v
    elseif i == 2
        m.c = v
    elseif i == 3
        m.b = v
    elseif i == 4
        m.d = v
    end
    return v
end

ints(m::AbstractGL₂) = int(m[1]), int(m[2]), int(m[3]), int(m[4])

function Base.:(==)(m::T, n::T) where T <: AbstractGL₂
    m = normalform!(m)
    n = normalform!(n)
    return ints(m) == ints(n)
end

function Base.hash(m::T, h::UInt) where {q, T<:AbstractGL₂{q}}
    m = normalform!(m)
    a,c,b,d = ints(m)
    val = d + q*(b + q*(c + q*a)) # q-adic expression of m
    return hash(T, hash(val, h))
end

function Base.:(*)(m::T, n::T) where T <: AbstractGL₂
    a,c,b,d = ints(m)
    A,C,B,D = ints(n)
    return T(a*A + b*C, c*A + d*C, a*B + b*D, c*B + d*D)
end

function mul!(x::Number, m::T) where T<:AbstractGL₂
    m[1] = x*int(m[1])
    m[2] = x*int(m[2])
    m[3] = x*int(m[3])
    m[4] = x*int(m[4])
    return m
end

LinearAlgebra.det(m::AbstractGL₂{q}) where q = ((a,c,b,d) = ints(m); GF{q}(a*d-c*b))

function Base.inv(m::T) where {q, T <: AbstractGL₂{q}}
    a,c,b,d = ints(m)
    p = a*d - b*c
    p¯¹ = invmod(p, q)

    return T(p¯¹*d, -p¯¹*c+q, -p¯¹*b+q, p¯¹*a)
end

isupper(m::AbstractGL₂) = iszero(m[2])

"""
    Bruhat(m::AbstractGL₂)
Return Bruhat decomposition of `m` consisting of four matrices: `u, w, D, U`,
where `m == u * w * D * U` and
* `u` and `U` are upper-triangular, with `1`s on the diagonal,
* `w` is either identity matrix, or `w = [0 1; -1 0]` corresponds to the Weyl group generator,
* `D` is diagonal.


If `b = Bruhat(m)` one can access those matrices through property query:
> `b.u, b.w, b.D, b.U`

In mathematical terms these matrices correspond to Bruhat decomposition of GL₂ as
double cosets of the Borel subgroup of upper-triangular matrices:

    `G = BWB = ∐_w BwB`

where `w` ranges over the Weyl group. In the case of `GL₂` we can write

    `GL₂ = DU ⊔ UwDU`

where `D` is the subgroup of diagonal matrices, `U` of unipotent ones
(thus `DU = B`) and `w² = -Id₂`.
"""
struct Bruhat{T<:AbstractGL₂}
    matrix::T
end

function Base.getproperty(bru::Bruhat{T}, S::Symbol) where T
    m = getfield(bru, :matrix)
    a,c,b,d = m[1], m[2], m[3], m[4]
    if S ∈ (:u, :w, :U, :D)
        if isupper(getfield(bru, :matrix)) # istrivial_weylcoset(bru)
            S === :u && return one(T)
            S === :w && return one(T)
            S === :U && return T(1, 0, b*inv(a), 1)
            S === :D && return T(a, 0, 0, d)
        else
            S === :u && return T(1, 0, a*inv(c), 1)
            S === :w && return T(0, -1, 1, 0)
            S === :D && return T(-c, 0, 0, -(a*d - b*c)*inv(c))
            S === :U && return T(1, 0, d*inv(c), 1)
        end
    else
        return getfield(bru, S)
    end
end

bruhat(m::AbstractGL₂) = (b = Bruhat(m); (b.u, b.w, b.D, b.U))

############################################
# GL₂{q}

mutable struct GL₂{q} <: AbstractGL₂{q}
    a::GF{q}
    c::GF{q}
    b::GF{q}
    d::GF{q}

    function GL₂{q}(a,c,b,d) where q
        @assert q isa Integer
        q > 1 || error(ArgumentError("$q (the modulus) must be > 1"))
        m = new{q}(a,c,b,d)
        @assert !iszero(det(m)) "Singular Matrix in GL₂{$q}: $m"
        return m
    end

    GL₂{q}(m::AbstractMatrix) where q = GL₂{q}(m[1,1], m[2,1], m[1,2], m[2,2])
end

isnormal(m::GL₂) = true

normalform!(m::GL₂) = m

order(::Type{GL₂{q}}) where q = (q^2 - 1)*(q^2 - q)


############################################
# SL₂{q}

mutable struct SL₂{q} <: AbstractGL₂{q}
    a::GF{q}
    c::GF{q}
    b::GF{q}
    d::GF{q}

    function SL₂{q}(a,c,b,d) where q
        @assert q isa Integer
        q > 1 || error(ArgumentError("$q (the modulus) must be > 1"))
        m = new{q}(a,c,b,d)
        @assert isone(det(m)) "Matrix of determinant ≠ 1 in SL₂{$q}: $m"
        return m
    end

    SL₂{q}(m::AbstractMatrix) where q = SL₂{q}(m[1,1], m[2,1], m[1,2], m[2,2])
end

isnormal(m::SL₂) = true

normalform!(m::SL₂) = m

order(::Type{SL₂{q}}) where q = q^3 - q


############################################
# PGL₂{q}

mutable struct PGL₂{q} <: AbstractGL₂{q}
    a::GF{q}
    c::GF{q}
    b::GF{q}
    d::GF{q}

    function PGL₂{q}(a,c,b,d) where q
        @assert q isa Integer
        q > 1 || error(ArgumentError("$q (the modulus) must be > 1"))
        m = new{q}(a,c,b,d)
        m = normalform!(m)
        @assert !iszero(det(m)) "Singular Matrix in PGL₂{$q}: $m"
        return m
    end

    PGL₂{q}(m::AbstractMatrix) where q = PGL₂{q}(m[1,1], m[2,1], m[1,2], m[2,2])
end

isnormal(m::PGL₂) = isone(m[1]) || (iszero(m[1]) && isone(m[2]))

function normalform!(m::PGL₂)
    isnormal(m) && return m
    if !iszero(m[1])
        a = inv(m[1])
    else
        a = inv(m[2])
    end
    m = mul!(int(a), m)
    return m
end

order(::Type{PGL₂{q}}) where q = q^3 - q


############################################
# PSL₂{q}

mutable struct PSL₂{q} <: AbstractGL₂{q}
    a::GF{q}
    c::GF{q}
    b::GF{q}
    d::GF{q}

    function PSL₂{q}(a,c,b,d) where q
        @assert q isa Integer
        q > 1 || error(ArgumentError("$q (the modulus) must be > 1"))
        m = new{q}(a,c,b,d)
        m = normalform!(m)
        @assert isone(det(m)) "Matrix of determinant ≠ 1 in PSL₂{$q}: $m"
        return m
    end

    PSL₂{q}(m::AbstractMatrix) where q = PSL₂{q}(m[1,1], m[2,1], m[1,2], m[2,2])
end

function isnormal(m::PSL₂{q}) where q
    isone(det(m)) || return false
    if 0 < m[1] <= div(q-1,2)
        return true
    elseif iszero(m[1]) && 0 < m[2] <= div(q-1,2)
        return true
    else
        return false
    end
end

function normalform!(m::PSL₂{q}) where q
    isnormal(m) && return m
    p = det(m)
    if isone(p)
        xinv = 1
    else
        x = sqrt(p)
        xinv = int(inv(x))
    end

    elt = ifelse(iszero(m[1]), m[2], m[1])

    if xinv*elt > div(q-1, 2)
        xinv = q - xinv
    end

    m = mul!(xinv, m)

    return m
end

order(::Type{PSL₂{q}}) where q = div(q^3 - q, 2)
