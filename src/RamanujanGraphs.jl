module RamanujanGraphs

using Primes
using LinearAlgebra
using LightGraphs

export lps_generators, cayley_graph, lps
export GF, GL₂, SL₂, PGL₂, PSL₂, Borel, Bruhat, bruhat
export PrincipalRepresentation, CosetDecomposition, right_action, order

include("cayley.jl")
include("gf.jl")
include("gl2.jl")
include("borel.jl")
include("coset_decomposition.jl")
include("representations.jl")

include("lps_generators.jl")

end # module
