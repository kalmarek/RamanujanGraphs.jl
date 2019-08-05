struct IntMod{q} <: Number
    x::Int
    function IntMod{q}(x) where q
        k = x % q;
        return (k >= 0 ? new{q}(k) : new{q}(k+q))
    end
end

Base.:(==)(n::IntMod{q}, m::IntMod{q}) where q = n.x == m.x
Base.hash(n::IntMod{q}, h::UInt) where q = hash(IntMod{q}, hash(n.x, h))

Base.:+(n::IntMod{q}, m::IntMod{q}) where q = IntMod{q}(n.x + m.x)
Base.:-(n::IntMod{q}, m::IntMod{q}) where q = IntMod{q}(n.x - m.x)
Base.:*(n::IntMod{q}, m::IntMod{q}) where q = IntMod{q}(n.x * m.x)

Base.:-(n::IntMod{q}) where q = IntMod{q}(q-n.x)
Base.inv(n::IntMod{q}) where q = IntMod{q}(invmod(n.x, q))

Base.zero(::Type{IntMod{q}}) where q = IntMod{q}(0)
Base.one(::Type{IntMod{q}}) where q = IntMod{q}(1)

Base.promote_rule(::Type{IntMod{q}}, ::Type{I}) where {q, I<:Integer} = IntMod{q}

Base.show(io::IO, n::IntMod) = print(io, n.x)

Int(n::IntMod) = n.x
isless(n::IntMod, y::Integer) = n.x < y


function Base.sqrt(n::IntMod{q}) where q
    l = legendresymbol(Int(n), q)
    l == 0 && return zero(IntMod{q})
    l == -1 && throw(ArgumentError("$(n.x) is not a square modulo $q"))
    for i in 1:q
        i^2 % q == n.x && return IntMod{q}(i)
    end
end

abstract type GL₂{q} <: AbstractMatrix{IntMod} end

Base.size(::GL₂) = (2,2)
Base.length(::GL₂) = 4

Base.one(::Type{T}) where T <: GL₂ = T(1, 0, 0, 1)
Base.one(::T) where T <: GL₂ = one(T)

Base.IndexStyle(::GL₂) = IndexLinear()

function Base.getindex(m::GL₂, i::Integer)
    i == 1 && return m.a
    i == 2 && return m.c # julia assumes column wise storage for IndexLinear
    i == 3 && return m.b
    i == 4 && return m.d
end

function Base.setindex!(m::GL₂{q}, v, i::Integer) where q
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

function Base.:(==)(m::T, n::T) where T <: GL₂
    m = normalform!(m)
    n = normalform!(n)
    return all(m[i] == n[i] for i in eachindex(m))
end

function Base.hash(m::T, h::UInt) where {q, T<:GL₂{q}}
    return UInt(h)
    m = normalform!(m)
    return hash(T, hash(q, hash(m[1], hash(m[2], hash(m[3], hash(m[4], h))))))
end

function Base.:(*)(m::T, n::T) where T <: GL₂
    a = m[1,1]*n[1,1] + m[1,2]*n[2,1]
    b = m[1,1]*n[1,2] + m[1,2]*n[2,2]
    c = m[2,1]*n[1,1] + m[2,2]*n[2,1]
    d = m[2,1]*n[1,2] + m[2,2]*n[2,2]
    return T(a,c,b,d)
end

LinearAlgebra.det(m::GL₂{q}) where q = m[1,1]*m[2,2] - m[2,1]*m[1,2]

function Base.inv(m::T) where {q, T <: GL₂{q}}
    p = det(m)
    p == 0 && throw(ArgumentError("Element is not invertible!\n $m"))
    p¯¹ = inv(p)

    return T(p¯¹*m.d, -p¯¹*m.c, -p¯¹*m.b, p¯¹*m.a)
end

############################################
#
# PGL₂{q}
#
############################################

mutable struct PGL₂{q} <: GL₂{q}
    a::IntMod{q}
    c::IntMod{q}
    b::IntMod{q}
    d::IntMod{q}

    function PGL₂{q}(a,c,b,d) where q
        @assert q isa Integer
        q > 1 || error("q (the modulus) must be > 1")
        m = new{q}(a,c,b,d)
        m = normalform!(m)
        @assert det(m) ≠ 0
        return m
    end

    PGL₂{q}(m::AbstractMatrix) where q = PGL₂{q}(m[1,1], m[2,1], m[1,2], m[2,2])
end

isnormal(m::PGL₂) = m[1] == 1 || (m[1] == 0 && m[2] == 1)

function normalform!(m::PGL₂)
    isnormal(m) && return m
    if m[1] ≠ 0
        a = inv(m[1])
    elseif m[2] ≠ 0
        a = inv(m[2])
    else
        error("The first column of $(typeof(m)) matrix must be non-zero! $m")
    end
    for i in eachindex(m)
        m[i] = a*m[i]
    end
    return m
end

order(::Type{PGL₂{q}}) where q = q^3 - q

############################################
#
# PSL₂{q}
#
############################################
mutable struct PSL₂{q} <: GL₂{q}
    a::IntMod{q}
    c::IntMod{q}
    b::IntMod{q}
    d::IntMod{q}

    function PSL₂{q}(a,c,b,d) where q
        @assert q isa Integer
        q > 1 || error("$q (the modulus) must be > 1")
        m = new{q}(a,c,b,d)
        m = normalform!(m)
        @assert det(m) == 1 "m = $m, det(m) = $(det(m))"
        @assert m[1] <= div(q-1,2) "m = $m is not in normal form!"
        return m
    end

    PSL₂{q}(m::AbstractMatrix) where q = PSL₂{q}(m[1,1], m[2,1], m[1,2], m[2,2])
end

function isnormal(m::PSL₂{q}) where q
    det(m) == 1 || return false
    if 0 < m[1] <= div(q-1,2)
        return true
    elseif m[1] == 0 && 0 < m[2] <= div(q-1,2)
        return true
    else
        return false
    end
end

function normalform!(m::PSL₂{q}) where q
    isnormal(m) && return m
    old_m = deepcopy(m)
    p = det(m)

    x = sqrt(p)
    xinv = inv(x)

    elt = (m[1] ≠ 0 ? m[1] : m[2])

    if xinv*elt > div(q-1, 2)
        xinv *= -1
    end

    for i in eachindex(m)
        m[i] = xinv*m[i]
    end
    return m
end

order(::Type{PSL₂{q}}) where q = div(q^3 - q, 2)
