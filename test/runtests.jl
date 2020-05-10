using Test, LTOptimize

for n in 1:20
    for k in 0:n
        try
            @test logfactorial(n, k, exact=true) ≈ log(factorial(n)/factorial(k))
        catch e
            println("Failure for (n, k) = $((n, k))")
            rethrow()
        end
    end
end

for n in 1:20
    for k in 0:n
        try
            @test isapprox(logfactorial(n, k, exact=true), logfactorial(n, k, exact=false), atol=1e-3)
        catch e
            println("Failure for (n, k) = $((n, k))")
            rethrow()
        end
    end
end

for n in 1:20
    for k in 0:n
        try
            @test isapprox(logbinomial(n, k), log(binomial(n, k)))
        catch e
            println("Failure for (n, k) = $((n, k))")
            rethrow()
        end
    end
end
