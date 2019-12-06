"""return log(n!/k!)"""
function logfactorial(n::Integer, k::Integer=1)
    n >= 0 || throw(DomainError(n, "n must be non-negative"))
    k >= 0 || throw(DomainError(k, "k must be non-negative"))
    k <= n || throw(DomainError(k, "k must be <= n"))
    if iszero(n) return 0.0 end
    rv = zero(n)
    for i in k+1:n
        rv += log(i)
    end
    return rv
end

"""return ≈log(n!/k!), computed using the method from Batir2010"""
function logfactorial_approx(n::Int, k::Int=1)
    n >= 0 || throw(DomainError(n, "n must be non-negative"))
    k >= 0 || throw(DomainError(k, "k must be non-negative"))
    k <= n || throw(DomainError(k, "k must be <= n"))
    if iszero(n) return 0.0 end

    r = 1/2 * log(2pi)
    r += n * log(n) - n

    # compute last term one step at a time while watching for overflow
    # r += 1/2 * log(n+1/6+1/(72n) - 31/(6480n^2) - 139/(155520n^3) + 9871 / (6531840n^4))
    r2 = n + 1/6
    tmp = 1/(72n)
    if !isnan(tmp) && !isinf(tmp) && tmp > 0
        r2 += tmp
    end
    tmp = 31/(6480n^2)
    if !isnan(tmp) && !isinf(tmp) && tmp > 0
        r2 -= tmp
    end
    tmp = 139/(155520n^3)
    if !isnan(tmp) && !isinf(tmp) && tmp > 0
        r2 -= tmp
    end
    tmp = 9871 / (6531840n^4)
    if !isnan(tmp) && !isinf(tmp) && tmp > 0
        r2 += tmp
    end
    r += 1/2 * log(r2)

    if k != one(k)
        r -= logfactorial_approx(k)
    end
    return r
end

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

"""plot the exact and approximate expressions for log(n!/k!)"""
function logfactorial_plot()
    ns = 2:10:10000
    ks = round.(Int, ns./2)
    exact = logfactorial.(ns, ks)
    approx = logfactorial_approx.(ns, ks)
    plt.plot(ns, exact, label="exact")
    plt.plot(ns, approx, "--", label="approx")
    plt.grid()
    plt.legend()
end
