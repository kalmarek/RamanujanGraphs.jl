struct Borel{T<:AbstractGL₂{q} where {q}} end

Borel(x::AbstractGL₂) = Borel{typeof(x)}()
Borel(::Type{T}) where {T} = Borel{T}()

Base.eltype(::Borel{T}) where {T} = T
Base.in(x::T, ::Borel{T}) where {T} = isupper(x)

Base.length(::Borel{SL₂{q}}) where {q} = q * (q - 1)

Base.iterate(B::Borel{SL₂{q}}) where {q} = one(SL₂{q}), (α = 1, invα = 1, u = 0)

function Base.iterate(B::Borel{SL₂{q}}, s) where {q}
    s.α == q - 1 && s.u == q - 1 && return nothing

    if s.u == q - 1
        u = 0
        α = s.α + 1
        invα = invmod(α, q)
    else
        u = s.u + 1
        α = s.α
        invα = s.invα
    end

    return SL₂{q}(α, 0, u, invα), (α = α, invα = invα, u = u)
end
