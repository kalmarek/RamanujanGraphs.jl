struct IntMod{q} <: Number
    value::Int16

    function IntMod{q}(n) where q
        @assert q > 1
        k = (-q<n<q ? n : n%q)
        k = ifelse(k >= 0, k, k + q)
        return new{q}(k)
    end
end

int(n::IntMod) = n.value

Base.:(==)(n::IntMod{q}, m::IntMod{q}) where q = int(n) == int(m)
# hash(RamanujanGraphs.IntMod) == 0x04fd9e474909f8bf
Base.hash(n::IntMod{q}, h::UInt) where q = xor(0x04fd9e474909f8bf, hash(q, hash(int(n), h)))

Base.:+(n::IntMod{q}, m::IntMod{q}) where q = IntMod{q}(int(n) + int(m))
Base.:-(n::IntMod{q}, m::IntMod{q}) where q = IntMod{q}(int(n) - int(m))
Base.:*(n::IntMod{q}, m::IntMod{q}) where q = IntMod{q}(int(n) * int(m))

Base.:-(n::IntMod{q}) where q = IntMod{q}(q - int(n))
Base.inv(n::IntMod{q}) where q = IntMod{q}(invmod(int(n), q))

function Base.:^(n::IntMod{q}, i::Integer) where q
   i < 0 && return inv(n)^-i
   return IntMod{q}(powermod(int(n), i, q))
end

Base.zero(::Type{IntMod{q}}) where q = IntMod{q}(0)
Base.one(::Type{IntMod{q}}) where q = IntMod{q}(1)
Base.iszero(n::IntMod) = int(n) == 0
Base.isone(n::IntMod) = int(n) == 1

Base.promote_rule(::Type{IntMod{q}}, ::Type{I}) where {q, I<:Integer} = IntMod{q}

# taken from ValidatedNumerics, under under the MIT "Expat" License:
# https://github.com/JuliaIntervals/ValidatedNumerics.jl/blob/master/LICENSE.md
function subscriptify(n::Integer)
    subscript_0 = Int(0x2080) # Char(0x2080) -> subscript 0
    @assert 0 <= n <= 9
    return Char(subscript_0 + n)
end

Base.show(io::IO, n::IntMod{q}) where q =
    print(io, "$(int(n))"*join(subscriptify.(reverse(digits(q))), ""))

import Base: <, <=
for ord in [:<, :(<=)]
    @eval begin
        $ord(n::IntMod{q}, m::IntMod{q}) where q = $ord(int(n),int(m))
        $ord(n::IntMod, y::Number) = $ord(int(n), y)
        $ord(y::Number, n::IntMod) = $ord(y, int(n))
    end
end

function legendresymbol(n, q)
    iszero(mod(n, q)) && return zero(n)
    isone(powermod(n, (q-1)÷2, q)) && return one(n)
    return -one(n)
end

Base.sqrt(n::IntMod{q}) where q = IntMod{q}(sqrtmod(int(n), q))

function sqrtmod(n::Integer, q::Integer)
    l = legendresymbol(n, q)
    l == 0 && return zero(n)
    l == -1 && throw(DomainError(n, "$n is not a square modulo $q"))
    for i in 1:q # bruteforce loop
        y = powermod(i, 2, q)
        y == n && return oftype(n, i)
    end
    return zero(n) # never hit, to keep compiler happy
end
