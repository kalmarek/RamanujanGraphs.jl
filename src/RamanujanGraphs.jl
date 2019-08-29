module RamanujanGraphs

using Primes
using LinearAlgebra
using Combinatorics
using LightGraphs

export PGL₂, PSL₂, lps_generators, cayley_graph, LPS

include("intmod.jl")
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

function diameter_ub(p::Integer, q::Integer)
    n = order(GLtype(p,q))
    return floor(Int, 2log(p, n) + 2log(p,2) + 1)
end

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

function cayley_graph(S::AbstractVector{T}, radius::Integer=10) where T
    @assert all((!isone).(S))
    S = reduced_generating_set(S)

    verts, sizes = generate_balls(unique([S; inv.(S)]), radius=radius)
    @show sizes

    cayley = SimpleGraph(length(verts))
    vlabels = Dict(g=>idx for (idx, g) in enumerate(verts))

    elabels = Dict{Tuple{Int, Int}, T}()

    for g in verts
        for s in S # S is reduced so that the graph is simple
            src = vlabels[g]
            dst = vlabels[g*s] # assuming it exists
            add_edge!(cayley, src, dst)
            elabels[(src, dst)] = s
        end
    end

    all(isequal(2length(S)), degree(cayley)) || @warn(
        "The degree is not constant = $(2length(S)). A truncated part of the graph is returned.")

    return cayley, verts, vlabels, elabels
end

cayley_graph(S::AbstractVector{T}) where T<:GL₂ = cayley_graph(order(T), S)

function cayley_graph(order::Integer, S::AbstractVector{T}) where T
    @assert all((!isone).(S))
    @assert all(inv(s) in S for s in S)

    sizes = Int[]
    cayley = SimpleGraph(1)
    verts = [one(T)]
    vlabels = Dict(one(T)=>1)
    elabels = Dict{Tuple{Int, Int}, T}()

    new_elts = Vector{T}(undef, length(S))
    nverts = 1

    for g in verts
        for i in 1:length(S)
            new_elts[i] = g*S[i]
        end

        for (s, elt) in zip(S, new_elts)
            if !haskey(vlabels, elt)
                add_vertex!(cayley)
                nverts +=1
                vlabels[elt] = nverts
                push!(verts, elt)
            end
            src, dst = vlabels[g], vlabels[elt]
            add_edge!(cayley, src, dst)
            elabels[(src, dst)] = s
        end
        nverts == order && break
    end

    verts_missing_edges = findall(!isequal(length(S)), degree(cayley))

    for v in verts_missing_edges
        nbhood = neighbors(cayley, v)
        seen_edges = T[]
        for n in nbhood
            if haskey(elabels, (v,n))
                push!(seen_edges, elabels[(v,n)])
            elseif haskey(elabels, (n, v))
                push!(seen_edges, inv(elabels[(n,v)]))
            end
        end

        K = setdiff(S, seen_edges)

        g = verts[v]
        src = vlabels[g]
        for s in K
            dst = vlabels[g*s]
            add_edge!(cayley, src, dst)
            elabels[(src, dst)] = s
        end
    end

    return cayley, verts, vlabels, elabels
end

end # module
