export logbinomial, krawtchouk, upperbound_ltfailure

"""return log(binomial(n, k))"""
function logbinomial(n::Integer, k::Integer)::Float64
    k <= n || return -Inf
    k != n || return zero(Float64)
    if k > (n - k)
        rv = logfactorial(n, n-k) - logfactorial(k)
    else
        rv = logfactorial(n, k) - logfactorial(n-k)
    end
    return rv
end

"""evaluate the krawtchouk polynomial"""
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

Return an upper bound on the decoding failure probability of LT codes
under optimal erasure decoding (Theorem 1 of Schotsch 2013).

γ: Inverse reception code rate, i.e., the number of received symbols
divided by the number of source symbols.

k: Number of source symbols.

q: Field size.

ds, ps: Represents the PDF of the degree distribution, with ps[i]
being the probability of a symbol having degree ds[i]. ps must sum to
1.

"""
function upperbound_ltfailure(γ::Real; k::Integer, q::Integer, ds, ps)::Float64
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

### lower bound - old, not sure if this is correct ###

"""inner term of ltfailure_lower_reference"""
function ltfailure_lower_reference_inner(i::Int; K::Int, ϵ::Number, Ω::Distribution{Univariate, Discrete})
    r = 0.0
    for d in 1:K
        println("d=$d, K=$K, f.K=$(Ω.K)")
        r += pdf(Ω, d) * binomial(K-i, d) / binomial(K, d)
    end
    r = r^(K*(1+ϵ))
    return r
end

"""

reference implementation of ltfailure_lower. use only for testing
ltfailure_lower.

"""
function ltfailure_lower_reference(K::Int, ϵ::Number, Ω::Distribution{Univariate, Discrete})
    println(Ω)
    r = 0.0
    for i in 1:K
        r += (-1)^(i+1)*binomial(K, i) * ltfailure_lower_reference_inner(
            i, K=K, ϵ=ϵ, Ω=Ω,
        )
    end
    return r
end

"""inner term of ltfailure_lower"""
function ltfailure_lower_inner(i::Int; K::Int, ϵ::Number, Ω::Distribution{Univariate, Discrete}) :: Float64
    r = zero(Float64)
    for d in 1:K
        if K-i < d || K < d
            continue
        end
        r += exp(log(pdf(Ω, d)) + logbinomial(K-i, d) - logbinomial(K, d))
        if isnan(r) # can't remember why I added this
            println("r=$r")
            println(pdf(Omega,d))
            println(log(pdf(Omega, d)))
            println(logbinomial(k-i, d))
            println(logbinomial(k, d))
            error("foo")
        end
    end
    return log(r) * (K*(1+ϵ))
end

"""
    ltfailure_lower(K::Int, ϵ::Number, Ω::Distribution{Univariate, Discrete})

Lower-bound the decoding failure probability of LT codes with K input
symbols, degree distribution Ω, and relative reception overhead ϵ
(i.e., ϵ=0.2 is a 20% reception overhead).

"""
function ltfailure_lower(ϵ::Real; K::Integer, Ω::Distribution{Univariate, Discrete})
    r = zero(Float64)
    tiny = eps(Float64)
    for i in 1:div(K, 2)
        v = zero(r)
        v += exp(logbinomial(K, 2i-1) + ltfailure_lower_inner(2i-1; K=K, ϵ=ϵ, Ω=Ω))
        v -= exp(logbinomial(K, 2i) + ltfailure_lower_inner(2i, K=K, ϵ=ϵ, Ω=Ω))

        # the terms become smaller. stop after reaching the smallest
        # normally represented float.
        if v < tiny
            break
        end
        r += v
    end
    return r
end

function plot_lower(K=100, M=K-2, δ=1e-2)
    Ω = Soliton(K, M, δ)
    ϵs = range(0, 0.1, length=10)
    # fs_ref = ltfailure_lower_reference(K, 0.2, Ω)
    fs = ltfailure_lower.(ϵs, K=K, Ω=Ω)
    # plt.semilogy(ϵs, fs_ref, "-s")
    plt.semilogy(ϵs, fs, "-o")
    plt.grid()
    plt.xlim(0, 0.1)
    plt.ylim(1e-6, 1)
    return fs
end
