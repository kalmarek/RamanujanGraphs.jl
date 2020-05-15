module RamanujanGraphs

using Primes
using LinearAlgebra
using LightGraphs

export lps_generators, cayley_graph, lps
export GF, GL₂, SL₂, PGL₂, PSL₂, Borel, Bruhat, bruhat
export issqrt, generator
export PrincipalRepr, CosetDecomposition, right_action, order

include("gf.jl")
include("gl2.jl")
include("borel.jl")
include("coset_decomposition.jl")
include("representations.jl")

include("cayley.jl")
include("lps_generators.jl")

end # module
