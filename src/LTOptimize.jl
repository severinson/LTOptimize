module LTOptimize

using PyPlot, Distributions

export logfactorial, logfactorial_approx, logbinomial

include("Soliton.jl")
include("Factorials.jl")
include("Bounds.jl")

end # module
