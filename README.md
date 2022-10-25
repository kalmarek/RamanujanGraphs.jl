[![Build Status](https://github.com/kalmarek/RamanujanGraphs.jl/workflows/CI/badge.svg)](https://github.com/kalmarek/RamanujanGraphs.jl/actions)
![ci](https://github.com/kalmarek/RamanujanGraphs.jl/workflows/ci/badge.svg)
[![codecov](https://codecov.io/gh/kalmarek/RamanujanGraphs.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/kalmarek/RamanujanGraphs.jl)

Following paper
> _Ramanujan graphs_ by Lubotzky, A., Phillips, R. & Sarnak, P. Combinatorica (1988) 8: 261. https://doi.org/10.1007/BF02126799

this package implements function `lps(p::Integer, q::Integer)` for different primes `p`,`q` congruent to `1 modulo 4` returning the appropriate Cayley graph.

A basic syntax is as follows:

```julia
using RamanujanGraphs
G, verts, vlabels, elabels = lps(p, q)
```

where
 * `G` is `(p+1)`-regular graph with
   - `q³ - q` vertices if `p` is not a square modulo `q` (Cayley graph of `PGL₂(q)`)
   - `(q³ - q)/2` vertices if `p` is a square modulo `q` (Cayley graph of `PSL₂(q)`)
 * `verts` is a plain array of vertices (=group elements)
 * `vlabels` is a labelling dictionary for vertices: group element pointing to its vertex in the graph
 * `elabels` is a dictionary for edges:
   - a tuple of integers `(src, dst)` points to a generator `g` of the group iff `verts[src]^-1*verts[dst] == g`, i.e. if we travel from `src` to `dst` by multiplying `src` by `g` on the right. Note that only one of `(src, dst)` and `(dst, src)` is stored.

Timings:

```julia
julia> using RamanujanGraphs, RamanujanGraphs.Graphs

julia> let (p,q) = (13,61)
           lps(p, q);
           @time G, verts, vlabels, elabels = lps(p, q);
           @assert nv(G) == RamanujanGraphs.order(eltype(verts))
           @info "Cayley graph of $(eltype(verts)):" degree=length(neighbors(G,1)) size=nv(G)
       end
 1.701829 seconds (2.05 M allocations: 272.833 MiB)
┌ Info: Cayley Graph of PSL₂{61}:
│   degree = 14
└   size = 113460

julia> let (p,q) = (13,73)
           lps(p, q);
           @time G, verts, vlabels, elabels = lps(p, q);
           @assert nv(G) == RamanujanGraphs.order(eltype(verts))
           @info "Cayley graph of $(eltype(verts)):" degree=length(neighbors(G,1)) size=nv(G)
       end
 6.400727 seconds (6.68 M allocations: 655.549 MiB)
┌ Info: Cayley Graph of PGL₂{73}:
│   degree = 14
└   size = 388944

```

# Irreducible representations for _SL₂(p)_

### Principal Series

These representations are associated to the induced representations of _B(p)_,
the _Borel subgroup_ (of upper triangular matrices) of _SL₂(p)_.
All representations of the Borel subgroup come from the representations of the
torus inside (i.e. diagonal matrices), hence are _1_-dimensional.

Therefore to define a matrix representation of _SL₂(p)_ one needs to specify:
 * a complex character of 𝔽ₚ (finite field of _p_ elements)
 * an explicit set of representatives of _SL₂(p)/B(p)_.

In code this can be specified by

```julia
p = 109 # our choice of a prime
ζ = root_of_unity((p-1)÷2, ...) # ζ is (p-1)÷2 -th root of unity
# two particular generators of SL₂(109):
a = SL₂{p}([0 1; 108 11])
b = SL₂{p}([57 2; 52 42])

S = [a, b, inv(a), inv(b)] # symmetric generating set
SL2p, _ = RamanujanGraphs.generate_balls(S, radius = 21)

Borel_cosets = RamanujanGraphs.CosetDecomposition(SL2p, Borel(SL₂{p}))
# the generator of 𝔽ₚˣ
α = RamanujanGraphs.generator(RamanujanGraphs.GF{p}(0))

ν₅ = let k = 5 # k runs from 0 to (p-1)÷4, or (p-3)÷4 depending on p (mod 4)
  νₖ = PrincipalRepr(
      α => ζ^k, # character sending α ↦ ζᵏ
      Borel_cosets
    )
end

```

### Discrete Series

These representations are associated with the action of _SL₂(p)_ (or in more
generality of _GL₂(p)_) on ℂ[𝔽ₚˣ], the vector space of complex valued functions
on 𝔽ₚˣ. There are however multiple choices how to encode such action.

Let _L_ = 𝔽ₚ(√_α_) be the unique quadratic extension of 𝔽ₚ by a square of a
generator _α_ of 𝔽ₚˣ. Comples characters of _Lˣ_ can be separated into
_decomposable_ (the ones that take constant 1 value on the unique cyclic
subgroup of order _(p+1)_ in _Lˣ_) and _nondecomposable_. Each _nondecomposable_
character corresponds to a representation of _SL₂(p)_ in discrete series.

To define matrix representatives one needs to specify
* _χ_:𝔽ₚ⁺ → ℂ, a complex, non-trivial character of the _additive group_ of 𝔽ₚ
* _ν_:_Lˣ_ → ℂ, a complex indecomposable character of _Lˣ_
* a basis for ℂ[𝔽ₚ].

Continuing the snippet above we can write

```julia
α = RamanujanGraphs.generator(RamanujanGraphs.GF{p}(0)) # a generator of 𝔽ₚˣ
β = RamanujanGraphs.generator_min(QuadraticExt(α))
# a generator of _Lˣ_ of minimal "Euclidean norm"

ζₚ = root_of_unity(p, ...)
ζ = root_of_unity(p+1, ...)

ϱ₁₇ = let k = 17 # k runs from 1 to (p-1)÷4 or (p+1)÷4 depending on p (mod 4)
    DiscreteRepr(
    RamanujanGraphs.GF{p}(1) => ζₚ, # character of the additive group of 𝔽ₚ
    β => ζ^k, # character of the multiplicative group of _L_
    basis = [α^i for i in 1:p-1] # our choice for basis: the dual of
)
```

A priori ζ needs to be a complex _(p²-1)_-th root of unity, however one can show
that a reduction to _(p+1)_-th Cyclotomic field is possible.
