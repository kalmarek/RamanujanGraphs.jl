function quadruples_4k_plus1(p::Integer)
    @assert p % 4 == 1
	@assert isprime(p)
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
	@assert isprime(p)
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

function quadruples(p::Integer)
	@assert p>0
	@assert isprime(p)
	if p % 4 == 1
		return quadruples_4k_plus1(p)
	elseif p % 4 == 3
		return quadruples_4k_plus3(p)
	end
end

generator(a₀,a₁,a₂,a₃,i) = [a₀ + i*a₁  a₂ + i*a₃;
			     			-a₂ + i*a₃ a₀ - i*a₁]

function PGLtype(p::Integer, q::Integer)
	legendre = legendresymbol(p,q)

    if legendre == -1
        return PGL₂{q}
    elseif legendre == 1
        return PSL₂{q}
	else
		throw("legendresymbol(p,q) = $legendre")
    end
end

function lps_generators(p::Integer, q::Integer)
	@assert p > 2
	@assert q > 2
	@assert p ≠ q
	@assert isprime(p)
	@assert isprime(q)
	@assert p % 4 == 1
	@assert q % 4 == 1

	i = sqrt(IntMod{q}(q-1))

	mats = [generator(a₀,a₁,a₂,a₃,i) for (a₀,a₁,a₂,a₃) in quadruples(p)]

	GL_t = PGLtype(p,q)
	S = GL_t.(mats)

	S = unique([S; inv.(S)])
	@assert all(inv(s) in S for s in S)
	@assert all((!isone).(S))
    return S
end
