using Test, LTOptimize

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
