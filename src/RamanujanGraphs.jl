module RamanujanGraphs

using Primes
using LinearAlgebra
using LightGraphs

export lps_generators, cayley_graph, lps
export GL₂, SL₂, PGL₂, PSL₂, Borel, Bruhat, bruhat

include("cayley.jl")
include("intmod.jl")
include("gl2.jl")

include("lps_generators.jl")

end # module
