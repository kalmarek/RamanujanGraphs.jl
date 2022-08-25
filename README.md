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
