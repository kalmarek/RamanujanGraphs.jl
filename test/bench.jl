using RamanujanGraphs
using LightGraphs

let (p,q) = (13,61)
    LPS(p, q);
    @time G, verts, vlabels, elabels = LPS(p, q);
    @assert nv(G) == RamanujanGraphs.order(eltype(verts))
    @info "Cayley Graph of $(eltype(verts)):" degree=length(neighbors(G,1)) size=nv(G)
end

let (p,q) = (13,73)
    LPS(p, q);
    @time G, verts, vlabels, elabels = LPS(p, q);
    @assert nv(G) == RamanujanGraphs.order(eltype(verts))
    @info "Cayley Graph of $(eltype(verts)):" degree=length(neighbors(G,1)) size=nv(G)
end

# using Pkg
# Pkg.test("RamanujanGraphs")
