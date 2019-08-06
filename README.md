Following paper
> _Ramanujan graphs_ by Lubotzky, A., Phillips, R. & Sarnak, P. Combinatorica (1988) 8: 261. https://doi.org/10.1007/BF02126799

this package implements function `LPS(p::Integer, q::Integer)` for different primes `p`,`q` congruent to `1 modulo 4` returning the appropriate Cayley graph.

A basic syntax is as follows:

```julia
using RamanujanGraphs
G, verts, vlabels, elabels = LPS(p, q)
```

where
 * `G` is `(p+1)`-regular graph with
   - `q³ - q` vertices if `p` is not a square modulo `q` (Cayley graph of `PGL₂(q)`)
   - `(q³ - q)/2` vertices if `p` is a square modulo `q` (Cayley graph of `PSL₂(q)`)
 * `verts` is a plain array of vertices (=group elements)
 * `vlabels` is a labelling dictionary for vertices: group element pointing to its vertex in the graph
 * elabels is a dictionary for edges:
   - a tuple of integers `(src, dst)` points to a GENERATOR `g` of the group iff `verts[src]^-1*verts[dst] == g`, i.e. if we travel from `src` to `dst` by multiplying `src` by `g` on the right.

Timings:

```julia

julia> using RamanujanGraphs

julia> using RamanujanGraphs.LightGraphs

julia> @time G, verts, vlabels, elabels = LPS(13, 61); @info "Cayley Graph of $eltype(verts):" degree=length(neighbors(G,1)) size=nv(G)
  2.597834 seconds (2.05 M allocations: 320.664 MiB, 20.40% gc time)
┌ Info: Cayley Graph of eltype(verts):
│   degree = 14
└   size = 113460

julia> @time G, verts, vlabels, elabels = LPS(13, 61); @info "Cayley Graph of $eltype(verts):" degree=length(neighbors(G,1)) size=nv(G)
  2.507784 seconds (2.05 M allocations: 320.664 MiB, 19.57% gc time)
┌ Info: Cayley Graph of eltype(verts):
│   degree = 14
└   size = 113460

julia> @time G, verts, vlabels, elabels = LPS(13, 73); @info "Cayley Graph of $(eltype(verts)):" degree=length(neighbors(G,1)) size=nv(G)
 10.958349 seconds (9.68 M allocations: 963.480 MiB, 18.57% gc time)
┌ Info: Cayley Graph of PGL₂{73}:
│   degree = 14
└   size = 388944

julia> @time G, verts, vlabels, elabels = LPS(13, 73); @info "Cayley Graph of $(eltype(verts)):" degree=length(neighbors(G,1)) size=nv(G)
  9.801834 seconds (6.68 M allocations: 811.452 MiB, 21.35% gc time)
┌ Info: Cayley Graph of PGL₂{73}:
│   degree = 14
└   size = 388944
```
