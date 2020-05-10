export logbinomial, krawtchouk, upperbound_ltfailure

"""Return log(binomial(n, k))."""
function logbinomial(n::Integer, k::Integer)::Float64
    k <= n || return -Inf
    k != n || return zero(Float64)
    if k > (n - k)
        return logfactorial(n) - logfactorial(n-k) - logfactorial(k)
    else
        return logfactorial(n) - logfactorial(k) - logfactorial(n-k)
    end
end

"""Evaluate the krawtchouk polynomial."""
function krawtchouk(ξ; ν::Integer, ς::Integer, q::Integer)
    0 <= ς <= ν || throw(ArgumentError("ς must be in [0, ν]"))
    ξ >= 0 || throw(DomainError(ξ, "ξ must be non-negative"))
    rv = zero(Float64)
    for i in 0:ς
        v1 = (ς-i) * log(q-1)
        v2 = logbinomial(ξ, i)
        if isinf(v2) continue end
        v3 = logbinomial(ν-ξ, ς-i)
        if isinf(v3) continue end
        v = v1+v2+v3
        if iszero(rem(i, 2))
            rv += exp(v)
        else
            rv -= exp(v)
        end
    end
    return rv
end

"""
    upperbound_ltfailure(γ::Real; k::Integer, q::Integer=2, ds, ps)::Float64

Return an upper bound on the decoding failure probability of LT codes
under optimal erasure decoding.

Source: Theorem 1, Analysis of LT Codes over Finite Fields under
Optimal Erasure Decoding, by Birgit Schotsch, Giuliano Garrammone and
Peter Vary, published in IEEE Communications Letters vol. 9, no. 17, 2013.

# Arguments

- `γ::Real`: Inverse reception code rate, i.e., the number of received
  symbols divided by the number of source symbols.

- `k::Integer`: Number of source symbols.

- `q::Integer=2`: Field size.

- `ds`, `ps`: Represents the PDF of the degree distribution, with
  `ps[i]` being the probability of a symbol having degree `ds[i]`.

"""
function upperbound_ltfailure(γ::Real; k::Integer, q::Integer=2, ds, ps)::Float64
    2 <= q || throw(DomainError(q, "q must be at least 2, which represents the binary field."))
    0 < k || throw(DomainError(k, "k must be positive."))
    length(ds) == length(ps) || throw(ArgumentError("ds and ps must be of equal length."))
    maximum(ds) <= k || raise(ArgumentError("Degrees larger than k are impossible."))
    sum(ps) ≈ 1 || raise(ArgumentError("Expected sum(ps)=1, but got $(sum(ps))."))
    rv = 0.0
    if γ < 1 return rv end
    for w in 1:k
        v1 = logbinomial(k, w)
        if isinf(v1) continue end
        v2 = (w-1) * log(q-1)
        v3 = 0.0 # inner term
        for (d, p) in zip(ds, ps)
            v = Float64(p)
            v *= krawtchouk(w; ν=k, ς=d, q=q)
            v /= krawtchouk(0; ν=k, ς=d, q=q)
            v3 += v
        end
        v3 *= (q-1)/q
        v3 += 1/q
        v3 = k*γ * log(v3)

        rv += exp(v1+v2+v3)
    end
    return max(min(rv, 1.0), 0.0)
end
