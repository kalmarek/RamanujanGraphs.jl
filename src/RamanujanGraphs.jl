module RamanujanGraphs

using Primes
using LinearAlgebra
using LightGraphs

export PGL₂, PSL₂, lps_generators, cayley_graph, lps

include("cayley.jl")
include("intmod.jl")
include("gl2.jl")
include("lps_generators.jl")

function lps(p::Integer, q::Integer)
    S = lps_generators(p,q)
    G, verts, vlabels, elabels = cayley_graph(order(eltype(S)), S)
    @assert order(eltype(S)) == length(verts)
    @assert all(isequal(p+1), degree(G))
    return G, verts, vlabels, elabels
end

function lps(p::Integer, q::Integer, radius::Integer)
    S = lps_generators(p,q)
    G, verts, vlabels, elabels = cayley_graph(S, radius=radius)
    radius < diameter_ub(p,q) && @warn "Radius given is smaller than its upper bound, cayley graph might be not complete!"
    return G, verts, labels, elabels
end

# from Lubotzky-Phillips-Sarnak
function diameter_ub(p::Integer, q::Integer)
    n = order(GLtype(p,q))
    return floor(Int, 2log(p, n) + 2log(p,2) + 1)
end

end # module
