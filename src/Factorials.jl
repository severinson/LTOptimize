"""return log(n!/k!)"""
function logfactorial(n::Integer, k::Integer=min(n, 1); exact=(n-k)<100)
    n >= 0 || throw(DomainError((n, k), "n must be non-negative"))
    k >= 0 || throw(DomainError((n, k), "k must be non-negative"))
    k <= n || throw(DomainError((n, k), "k must be <= n"))
    if iszero(n) return 0.0 end
    if !exact return logfactorial_approx(n, k) end
    rv = zero(n)
    for i in k+1:n
        rv += log(i)
    end
    return rv
end

"""return ≈log(n!/k!), computed using the method from Batir2010"""
function logfactorial_approx(n::Integer, k::Integer=min(n, 1))
    n >= 0 || throw(DomainError((n, k), "n must be non-negative"))
    k >= 0 || throw(DomainError((n, k), "k must be non-negative"))
    k <= n || throw(DomainError((n, k), "k must be <= n"))
    if iszero(n) return 0.0 end
    rv = 1/2 * log(2pi)
    rv += n * log(n) - n
    # the final term is
    # 1/2 * log(n+1/6+1/(72n) - 31/(6480n^2) - 139/(155520n^3) + 9871 / (6531840n^4))
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
    rv += 1/2 * log(r2)

    if k > 1
        rv -= logfactorial(k)
    end
    return rv
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
function plot_logfactorial()
    ns = 2:10:1000000
    ks = round.(Int, ns./2)
    # exact = logfactorial.(ns, ks, exact=true)
    approx = logfactorial.(ns, ks, exact=false)
    lib = SpecialFunctions.logfactorial.(ns) .- SpecialFunctions.logfactorial.(ks)
    # plt.plot(ns, exact, label="exact")
    plt.plot(ns, approx, "-", label="approx")
    plt.plot(ns, lib, "--", label="lib")
    plt.grid()
    plt.legend()
end

function plot_logbinomial()
    ns = 2:10:1000
    ks = round.(Int, ns.*0.9)
    mine = logbinomial.(ns, ks)
    println(mine)
    lib = logabsbinomial.(ns, ks)
    println([x[1] for x in lib])
    # return
    plt.plot(ns, mine, label="mine")
    plt.plot(ns, [x[1] for x in lib], ".", label="lib")
    plt.grid()
    plt.legend()
end
