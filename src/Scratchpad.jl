using PyPlot
using SpecialFunctions

include("Soliton.jl")

"""plot the robust Soliton distribution"""
function plot_distribution(k=600, M=40, δ=1e-6)
    Ω = Soliton(k, M, δ)
    xs = 1:1:k
    plt.plot(xs, pdf.(Ω, xs))
    return sum(pdf.(Ω, xs))
end

function plot_bounds(k=600, M=40, δ=1e-6)
    Ω = Soliton(k, M, δ)
    ds = [d for d in 1:k if pdf(Ω, d)>1e-3]
    ps = pdf.(Ω, ds)
    ps ./= sum(ps)

    γs = range(1, 1.04, length=20)
    for q in [2, 3, 4, 5, 6, 7, 8]
        @views fs = upperbound_ltfailure.(γs, k=k, q=q, ds=ds, ps=ps)
        plt.semilogy(k.*(γs.-1), fs, ".-", label="GF($q)")
    end
    fs = lowerbound_ltfailure.(γs, k=k, ds=ds, ps=ps)
    plt.semilogy(k.*(γs.-1), fs, "k-", label="LB")

    plt.grid()
    plt.xlim(0, 15)
    plt.ylim(1e-4, 1)
    plt.xlabel("Reception overhead")
    plt.ylabel("Decoding failure probability")
    plt.legend()
    plt.tight_layout()
    return
end

function plot_bounds_uniform(k=50)

    ds = 1:k
    ps = [1/k for _ in 1:k]

    γs = range(1, 1.5, length=20)
    for q in [2, 4, 6, 8]
        println("q=$q")
        @views fs = upperbound_ltfailure.(γs, k=k, q=q, ds=ds, ps=ps)
        plt.semilogy(k.*(γs.-1), fs, ".-", label="GF($q)")
    end
    fs = lowerbound_ltfailure.(γs, k=k, ds=ds, ps=ps)
    plt.semilogy(k.*(γs.-1), fs, "k-", label="LB")

    plt.grid()
    # plt.xlim(0, 12)
    # plt.ylim(1e-4, 1)
    plt.xlabel("Reception overhead")
    plt.ylabel("Decoding failure probability")
    plt.legend()
    plt.tight_layout()
    return
end
