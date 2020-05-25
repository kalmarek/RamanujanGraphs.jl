struct CosetDecomposition{T,TC}
    representatives::Vector{T}
    inv_representatives::Vector{T}
    trivial_coset::TC
end

Base.length(cd::CosetDecomposition) = length(cd.representatives)
Base.eltype(cd::CosetDecomposition{T}) where {T} = T

function Base.getindex(cd::CosetDecomposition, i::Integer)
    @boundscheck 0 < abs(i) <= length(cd) || throw(BoundsError(cd, i))
    i < 0 && return cd.inv_representatives[-i]
    return cd.representatives[i]
end
"""
    CosetDecomposition(G, H)
Decomposes `G` into __right__ cosets of `H`, usually written `H\\G`.
The multiplication is given by `Hg * Hk = H(g*h)` with `Hg == Hk` if and only if
`g*h^-1 ∈ H`.

The constructor uses:
* iteration over `G` hence `G` must implement iterator protocol,
* membership testing for `H`, i.e. `Base.in(g, H)` must be overloaded.

For convenience `cd::CosetDecomposition` implements `getindex(cd, i::Integer)`,
so that if `i` is a positive integer
 * `cd[i]` returns the representative of `i`-th coset
 * `cd[-i]` returns the inverse of `c[i]` (stored, not evaluated).
"""
function CosetDecomposition(group, subgroup)
    ordg = length(group)
    ordh = length(subgroup)
    k, r = divrem(ordg, ordh)
    @assert r == 0

    coset_reps = [one(first(group))]
    coset_reps_inv = [one(first(group))]
    sizehint!(coset_reps, k)
    sizehint!(coset_reps_inv, k)

    i = 1
    cosets_found = 1

    for g in group
        any(g * coset_reps_inv[i] ∈ subgroup for i = 1:cosets_found) && continue
        cosets_found += 1
        push!(coset_reps, g)
        push!(coset_reps_inv, inv(g))
        cosets_found == k && break
    end

    @assert cosets_found == k
    @assert length(Set(coset_reps)) == length(coset_reps)
    return CosetDecomposition(coset_reps, coset_reps_inv, subgroup)
end

"""
    right_action(x, cd::CosetDecomposition)
Return the permutation `p::Vector{Int}` representing the __right__ action
of `x` on the set of cosets of `cd`, that is
> `p[i] == j` iff `H*cd[i]*g == H*cd[j]`, or `cd[i]*g*cd[-j] ∈ H`.

There is an in-place version (`right_action!`) as well.
"""
function right_action(x::AbstractGL₂, cd::CosetDecomposition)
    perm = zeros(Int, length(cd))
    perm = right_action!(perm, x, cd)
    @assert all(!iszero, perm)
    return perm
end

function right_action!(
    perm::AbstractVector{<:Integer},
    g::AbstractGL₂,
    cosets::CosetDecomposition,
)

    # i → j when
    # Hc_i*g = Hc_j, i.e. c_i*g*c_j^-1 ∈ B

    for i = 1:length(cosets)
        cig = cosets[i] * g
        for j = 1:length(cosets)
            if cig * cosets[-j] ∈ cosets.trivial_coset
                perm[i] = j
                break
            end
        end
    end
    return perm
end
