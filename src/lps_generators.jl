function quadruples_4k_plus1(p::Integer)
    @assert p % 4 == 1
    N = floor(Int, sqrt(p))
    N = iseven(N) ? N : N+1
    quads = NTuple{4, Int}[]
    for a in 1:2:N
        s1 = a^2
        for b in -N:2:N
            s2 = s1 + b^2
            for c in -N:2:N
                s3 = s2 + c^2
                for d in -N:2:N
                    s4 = s3 + d^2
                    s4 == p && push!(quads, (a,b,c,d))
                end
            end
        end
    end
    return quads
end

function quadruples_4k_plus3(p::Integer)
    @assert p % 4 == 3
    N = floor(Int, sqrt(p))
    N = iseven(N) ? N+1 : N
    quads = NTuple{4, Int}[]
    for a in 0:2:N
        s1 = a^2
        b_range = (a == 0 ? StepRange(1,2,N) : StepRange(-N,2,N))
        for b in b_range
            s2 = s1 + b^2
            for c in -N:2:N
                s3 = s2 + c^2
                for d in -N:2:N
                    s4 = s3 + d^2
                    s4 == p && push!(quads, (a,b,c,d))
                end
            end
        end
    end
    return quads
end

function hyperboloid_solution(q::Integer, b=1)
    for x in 0:q
        x² = x^2
        for y in 0:q
            y² = y^2
            iszero((x² + y² + b) % q) && return x, y
        end
    end
    throw("a solution should always exist in finite field!")
end

generator(a,b,c,d, x,y) = [a + b*x + d*y -b*y + c + d*x
			     -b*y - c + d*x a - b*x- d*y]

function lps_generators(p::Integer, q::Integer)
    x,y = hyperboloid_solution(q)
    @assert (x^2 + y^2 + 1) %q == 0

    if p % 4 == 1
        mats = [generator(a,b,c,d,x,y) for (a,b,c,d) in quadruples_4k_plus1(p)]
    elseif p % 4 == 3
        mats = [generator(a,b,c,d,x,y) for (a,b,c,d) in quadruples_4k_plus3(p)]
    end

    legendre = legendresymbol(p,q)
    if legendre == -1
        S = PGL₂{q}.(mats)
    elseif legendre == 1
        S = PSL₂{q}.(mats)
    end

	S = unique([S; inv.(S)])
	@assert all(inv(s) in S for s in S)
	@assert all((!isone).(S))
    return S
end
