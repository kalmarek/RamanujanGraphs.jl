module RamanujanGraphs

using Primes
using LinearAlgebra
using Combinatorics
using LightGraphs

export PGL₂, PSL₂, lps_generators, cayley_graph, LPS

include("cayley.jl")
include("intmod.jl")
include("gl2.jl")
include("lps_generators.jl")

cayley_graph(S::AbstractVector{T}) where T<:GL₂ = cayley_graph(order(T), S)

function LPS(p::Integer, q::Integer, radius::Integer)
    S = lps_generators(p,q)
    G, verts, vlabels, elabels = cayley_graph(S, radius=radius)
    # @assert order(eltype(verts)) == length(verts)
    # @assert all(isequal(p+1), degree(G))
    return G, verts, labels, elabels
end

function LPS(p::Integer, q::Integer)
    S = lps_generators(p,q)
    G, verts, vlabels, elabels = cayley_graph(S)
    @assert order(eltype(S)) == length(verts)
    @assert all(isequal(p+1), degree(G))
    return G, verts, vlabels, elabels
end

# from Lubotzky-Phillips-Sarnak
function diameter_ub(p::Integer, q::Integer)
    n = order(GLtype(p,q))
    return floor(Int, 2log(p, n) + 2log(p,2) + 1)
end

end # module
