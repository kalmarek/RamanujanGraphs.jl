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

abstract type GL₂{q} <: AbstractMatrix{Int} end

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

function Base.:(==)(m::T, n::T) where T <: GL₂
    return all(normalform!(m)[i] == normalform!(n)[i] for i in eachindex(m))
end

function Base.hash(m::PT, h::UInt) where {q, PT<:GL₂{q}}
    m = normalform!(m)
    return hash(PT, hash(q, hash(m.a, hash(m.b, hash(m.c, hash(m.d, h))))))
end

function Base.:(*)(m::T, n::T) where T <: GL₂
    a = m[1,1]*n[1,1] + m[1,2]*n[2,1]
    b = m[1,1]*n[1,2] + m[1,2]*n[2,2]
    c = m[2,1]*n[1,1] + m[2,2]*n[2,1]
    d = m[2,1]*n[1,2] + m[2,2]*n[2,2]
    return normalform!(T(a,c,b,d))
end

LinearAlgebra.det(m::GL₂{q}) where q = modd(m[1,1]*m[2,2] - m[2,1]*m[1,2], q)

function Base.inv(m::T) where {q, T <: GL₂{q}}
    D = modd(det(m), q)
    D == 0 && throw(ArgumentError("Element is not invertible!\n $m"))
    D¯¹ = invmod(D, q)

    return T(D¯¹*m.d, -D¯¹*m.c, -D¯¹*m.b, D¯¹*m.a)
end

############################################
#
# PGL₂{q}
#
############################################

mutable struct PGL₂{q} <: GL₂{q}
    a::Int
    c::Int
    b::Int
    d::Int

    function PGL₂{q}(a,c,b,d) where q
        @assert q isa Integer
        q > 1 || error("q (the modulus) must be > 1")
        a,c,b,d = modd(a,q), modd(c,q), modd(b,q), modd(d, q)
        m = new{q}(a,c,b,d)
        m = normalform!(m)
        @assert det(m) ≠ 0
        return m
    end

    PGL₂{q}(m::AbstractMatrix{<:Integer}) where q = PGL₂{q}(m[1,1], m[2,1], m[1,2], m[2,2])
end

isnormal_pgl2(a,c,b,d) = a == 1 || (a == 0 && c == 1)
isnormal(m::PGL₂) = isnormal_pgl2(m...)

function normalform_pgl2(q, a,c,b,d)
    isnormal_pgl2(a,c,b,d) && return (a,c,b,d)
    if a ≠ 0
        a¯¹ = invmod(a, q)
        return modd.((1, a¯¹*c, a¯¹*b, a¯¹*d), q)
    elseif a == 0
        c¯¹ = invmod(c, q)
        return modd.((0,     1, c¯¹*b, c¯¹*d), q)
    else
        error("The first element of quadruple should be non-negative")
    end
end

function normalform!(m::PGL₂{q}) where q
    isnormal(m) && return m
    m.a, m.c, m.b, m.d = normalform_pgl2(q, m...)
    return m
end

order(::Type{PGL₂{q}}) where q = q^3 - q

############################################
#
# PSL₂{q}
#
############################################
mutable struct PSL₂{q} <: GL₂{q}
    a::Int
    c::Int
    b::Int
    d::Int

    function PSL₂{q}(a,c,b,d) where q
        @assert q isa Integer
        q > 1 || error("$q (the modulus) must be > 1")
        a,c,b,d = modd(a,q), modd(c,q), modd(b,q), modd(d, q)
        m = new{q}(a,c,b,d)
        m = normalform!(m)
        @assert det(m) == 1
        return m
    end

    PSL₂{q}(m::AbstractMatrix{<:Integer}) where q = PSL₂{q}(m[1,1], m[2,1], m[1,2], m[2,2])
end

isnormal_psl2(q, a,c,b,d) = modd(a*d - c*b, q) == 1 && a <= div(q-1,2)
isnormal(m::PSL₂{q}) where q = isnormal_psl2(q, m...)

function normalform_psl2(q, a,c,b,d)
    isnormal_psl2(q, a,c,b,d) && return (a,c,b,d)
    D = a*d - c*b
    sqrtD = sqrtmod(D, q)
    sqrtD¯¹ = invmod(sqrtD, q)
    if modd(sqrtD¯¹*a, q) > div(q-1,2)
        sqrtD¯¹ *= -1
    end
    return modd.((sqrtD¯¹*a, sqrtD¯¹*c, sqrtD¯¹*b, sqrtD¯¹*d), q)
end

function normalform!(m::PSL₂{q}) where q
    isnormal(m) && return m
    m.a, m.c, m.b, m.d = normalform_psl2(q, m...)
    return m
end

order(::Type{PSL₂{q}}) where q = div(q^3 - q, 2)
