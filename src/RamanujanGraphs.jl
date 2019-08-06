module RamanujanGraphs

using Primes
using LinearAlgebra
using Combinatorics
using LightGraphs

export PGL₂, PSL₂, lps_generators, cayley_graph

include("gl2.jl")
include("lps_generators.jl")

function generate_balls(S::AbstractVector{T}, Id::T=one(T);
        radius=2) where T
    sizes = Int[]
    B = [Id]
    B_Set = Set([Id])

    for i in 1:radius
        new_elts = T[]
        for (g,h) in Base.product(B,S)
            elt = g*h
            g*h in B_Set && continue
            push!(new_elts, elt)
            push!(B_Set, elt)
        end

        if length(new_elts) == 0
            @info "Given radius = $radius, but the group already saturated at" radius=i-1 size=sizes[end]
            break
        end
        B = append!(B, new_elts)
        push!(sizes, length(B))
    end
    return B, sizes
end

function reduced_generating_set(S::AbstractVector)
    @assert all((!isone).(S))
    rS = Set{eltype(S)}()
    for s in S
        s^2 == s && push!(rS, s)
        inv(s) in rS && continue
        push!(rS, s)
    end
    return collect(rS)
end

function cayley_graph(S::AbstractVector{T}; radius::Integer=3) where T
    @assert all((!isone).(S))
    # @assert all(inv(s) in S for s in S)

    vertices, sizes = generate_balls(unique([S; inv.(S)]), radius=radius)
    @info "Generated balls of radii" sizes

    cayley = SimpleGraph(length(vertices))
    vlabels = Dict(g=>idx for (idx, g) in enumerate(vertices))

    elabels = Dict{Tuple{Int, Int}, T}()

    for (idx,g) in enumerate(vertices)
        for s in S
            src = vlabels[g]
            dst = vlabels[g*s] # assuming it exists
            add_edge!(cayley, src, dst)
        end
    end

    return cayley, vertices, vlabels
end

end # module
