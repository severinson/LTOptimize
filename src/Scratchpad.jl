using PyPlot

"""plot the robust Soliton distribution"""
function plot_distribution(k=600, M=40, δ=1e-6)
    Ω = Soliton(k, M, δ)
    xs = 1:1:k
    plt.plot(xs, pdf.(Ω, xs))
    return sum(pdf.(Ω, xs))
end

"""plot the upper bound"""
function plot_upper(k=600, M=40, δ=1e-6)
    Ω = Soliton(k, M, δ)
    ps = pdf.(Ω, 1:k)
    ds = [d for d in 1:k if ps[d]>1e-3]
    @views ps[ds] ./= sum(ps[ds])
    γs = range(1, 1.1, length=20)
    for q in [2, 4, 8]
        @views fs = upperbound_ltfailure.(γs, k=k, q=q, ds=ds, ps=ps[ds])
        plt.semilogy(k.*(γs.-1), fs, ".-", label="GF($q)")
    end
    plt.grid()
    # plt.xlim(γs[1], γs[end])
    plt.ylim(1e-6, 1)
    plt.xlabel("Reception overhead")
    plt.ylabel("Decoding failure probability")
    plt.legend()
    plt.tight_layout()
    return
end
