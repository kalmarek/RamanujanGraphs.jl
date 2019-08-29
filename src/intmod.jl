struct IntMod{q} <: Number
    value::Int16

    function IntMod{q}(n) where q
        @assert q > 1
        k = n % q
        k = ifelse(k >=0, k, k+q)
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

Base.zero(::Type{IntMod{q}}) where q = IntMod{q}(0)
Base.one(::Type{IntMod{q}}) where q = IntMod{q}(1)
Base.iszero(n::IntMod) = int(n) == 0
Base.isone(n::IntMod) = int(n) == 1

Base.promote_rule(::Type{IntMod{q}}, ::Type{I}) where {q, I<:Integer} = IntMod{q}


Combinatorics.legendresymbol(n::IntMod{q}) where q = legendresymbol(int(n), q)

function Base.sqrt(n::IntMod{q}) where q
    l = legendresymbol(n)
    l == 0 && return zero(IntMod{q})
    l == -1 && throw(DomainError(n, "$(int(n)) is not a square modulo $q"))
    for i in 1:q
        i^2 % q == int(n) && return IntMod{q}(i)
    end
    return zero(n) # never hit, to keep compiler happy
end
