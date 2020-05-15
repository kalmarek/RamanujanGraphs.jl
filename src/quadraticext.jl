struct QuadraticExt{ε,T} <: Number
    reim::Tuple{T,T}
end
QuadraticExt(ε::T) where {T} = QuadraticExt{ε, T}((zero(ε),zero(ε)))
QuadraticExt{ε}(re::T, im::T) where {ε,T} = QuadraticExt{ε, T}((re, im))
QuadraticExt{ε, T}(re::T) where {ε,T} = QuadraticExt{ε, T}((re, zero(re)))

Base.reim(x::QuadraticExt) = x.reim
Base.real(x::QuadraticExt) = first(reim(x))
Base.imag(x::QuadraticExt) = last(reim(x))
Base.isreal(x::QuadraticExt) = iszero(imag(x))

sqroot(x::QuadraticExt{ε}) where {ε} = ε

Base.show(io::IO, x::QuadraticExt) =
    print(io, real(x), " + ", imag(x), "√", sqroot(x))

Base.:(==)(x::QuadraticExt{ε}, y::QuadraticExt{ε}) where {ε} = reim(x) == reim(y)

Base.hash(x::QuadraticExt, h::UInt) =
    hash(reim(x), hash(sqroot(x), hash(QuadraticExt, h)))

Base.:+(x::QuadraticExt{ε}, y::QuadraticExt{ε}) where {ε} =
    QuadraticExt{ε}((reim(x) .+ reim(y))...)

Base.:-(x::QuadraticExt{ε}, y::QuadraticExt{ε}) where {ε} =
    QuadraticExt{ε}((reim(x) .- reim(y))...)

Base.:-(x::QuadraticExt{ε}) where {ε} = QuadraticExt{ε}((.-reim(x))...)

function Base.:*(x::QuadraticExt{ε}, y::QuadraticExt{ε}) where {ε}
    a, b = reim(x)
    c, d = reim(y)
    return QuadraticExt{ε}(a*c + ε*b*d, a*d + b*c)
end

function Base.inv(x::QuadraticExt{ε}) where {ε}
    (a, b) = reim(x)
    if iszero(a)
        c = a
        d = inv(b * ε)
    else
        Nx = norm(x)
        c = a * inv(Nx)
        d = -(Nx * b * inv(a))
    end
    return QuadraticExt{ε}(c, d)
end

Base.:/(x::QuadraticExt{ε}, y::QuadraticExt{ε}) where ε = x*inv(y)

Base.:+(x::QuadraticExt, a::Number) =
    QuadraticExt{sqroot(x)}(real(x) + a, imag(x))
Base.:+(a::Number, x::QuadraticExt) = x + a
Base.:*(x::QuadraticExt, a::Number) = QuadraticExt{sqroot(x)}((a .* reim(x))...)
Base.:*(a::Number, x::QuadraticExt) = x * a
Base.:/(x::QuadraticExt, a::Number) = QuadraticExt{sqroot(x)}((reim(x) ./ a)...)

Base.promote_rule(::Type{T}, ::Type{QuadraticExt{ε,T}}) where {ε, T} =
    QuadraticExt{ε,T}

LinearAlgebra.norm(x::QuadraticExt) = real(x)^2 - sqroot(x) * imag(x)^2

Base.conj(x::QuadraticExt) = QuadraticExt{sqroot(x)}(real(a), -imag(b))

Base.zero(x::QuadraticExt) = QuadraticExt{sqroot(x)}(zero.(reim(x))...)

Base.one(x::QuadraticExt) =
    (v = real(x); QuadraticExt{sqroot(x)}(one(v), zero(v)))

oneim(x::QuadraticExt) =
    (v = real(x); QuadraticExt{sqroot(x)}(zero(v), one(v)))

Base.iszero(x::QuadraticExt) = iszero(real(x)) && iszero(imag(x))
Base.isone(x::QuadraticExt) = isone(real(x)) && iszero(imag(x))

struct Units{ε,T} end
Units(x::QuadraticExt{ε,T}) where {ε,T<:GF} = Units{ε,T}()

Base.eltype(u::Units{ε,T}) where {ε,T} = QuadraticExt{ε,T}
Base.length(u::Units{ε,T}) where {ε,T} = order(T) + 1

function Base.iterate(u::Units{ε,T}, state=(zero(T), 1)) where {ε,T}
    t, count = state
    if count > order(T)
        return QuadraticExt{ε}(-one(t), zero(t)), nothing
    end
    v = inv(one(t) - ε * t^2)
    x = (one(t) + ε * t^2) * v
    y = (t + t) * v

    return QuadraticExt{ε}(x, y), (t + one(t), count + 1)
end

Base.iterate(u::Units{ε,T}, ::Nothing) where {ε,T} = nothing

function generator(x::QuadraticExt{ε,<:GF}) where {ε}
    y = _elt_of_norm(x, ε)
    q = Int(int(ε))

    for u in Units(x)
        w = y*u
        for i in 1:q^2-2
            isone(w^i) && break
        end
        return w
    end
    return zero(x) # never returned to keep compiler happy
end

function _elt_of_norm(x::QuadraticExt{ε,T}, n::T) where {ε,T<:GF}
    iszero(n) && return zero(n)
    isone(n) && return one(n)
    for d in T
        c² = n + ε*d^2
        if RamanujanGraphs.issqrt(c²)
            y = QuadraticExt{ε}(sqrt(c²), d)
            @assert norm(y) == n
            return y
        end
    end
    return zero(x) # never returned to keep compiler happy
end
