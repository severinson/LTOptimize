using PyPlot

"""evaluate the krawtchouk polynomial"""
function krawtchouk(ξ; ν::Integer, ς::Integer, q::Integer)
    0 <= ς <= ν || throw(ArgumentError("ς must be in [0, ν]"))
    ξ >= 0 || throw(DomainError(ξ, "ξ must be non-negative"))
    rv = zero(Float64)
    for i in 0:ς
        v1 = (ς-i) * log(q-1)
        # println("v1=$v1")
        v2 = logbinomial(ξ, i)
        # println("v2=$v2, ξ=$ξ, i=$i")
        if isinf(v2) continue end
        v3 = logbinomial(ν-ξ, ς-i)
        # println("v3=$v3")
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

"""lt code upper bound on failure probability from Schotsch 2013 (1)"""
function ltfailure_upper(γ; k, Ω, q)
    rv = 0.0
    for w in 1:K
        v1 = logbinomial(k, w)
        # println("v1=$v1")
        if isinf(v1) continue end

        v2 = (w-1) * log(q-1)

        # inner term (sum over the distribution)
        v3 = 0.0
        for d in 1:k
            v = pdf(Ω, d)
            v *= krawtchouk(w; ν=k, ς=d, q=q)
            v /= krawtchouk(0; ν=k, ς=d, q=q)
            # println("inner v=$v")
            v3 += v
        end
        v3 *= (q-1)/q
        v3 += 1/q
        v3 = k*γ * log(v3)
        # println("v3=$v3")

        v = v1+v2+v3
        # println("v=$v")
        rv += exp(v)
    end
    return max(min(rv, 1.0), 0.0)
end

"""plot the robust soliton"""
function plot_distribution(K=300, M=100, δ=1e-3)
    Ω = Soliton(K, M, δ)
    xs = 1:1:K
    plt.plot(xs, pdf.(Ω, xs))
    return sum(pdf.(Ω, xs))
end

"""plot the upper bound"""
function plot_upper(K=300, M=50, δ=1e-1/2)
    Ω = Soliton(K, M, δ)
    q = 2
    γs = range(1.0, 1.1, length=10)
    fs = ltfailure_upper.(γs, k=K, Ω=Ω, q=q)
    plt.semilogy(γs, fs)
    plt.grid()
    plt.xlim(1.0, 1.1)
    plt.ylim(1e-6, 1)
    return
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
