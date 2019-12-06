using Test, LTOptimize

for n in 1:20
    for k in 0:n
        @test isapprox(log(factorial(n)/factorial(k)), logfactorial(n, k))
    end
end

for n in 1:20
    for k in 0:n
        @test isapprox(log(binomial(n, k)), logbinomial(n, k))
    end
end
