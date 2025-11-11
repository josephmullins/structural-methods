using Optim
using Distributions

function solve_consumption(r,α,W,A,σ,ψ)
    Q = 1/ (1 - 1/(1+r))
    f(c) = (Q * c - Q * W^(1 + ψ) * α^ψ * c^(-σ*ψ) - A)^2
    r = Optim.optimize(f,0.,A+W)
    return r.minimizer
end
function simulate_data(σ,ψ,r,γ,N)
    ch = [0.3 0. 0.; 0.5 0.5 0.; 0.4 0.8 1.8]
    Σ = ch * ch'
    X = rand(MvNormal(Σ),N)
    Z = rand(N) .< 0.5
    α = exp.(X[1,:])
    W = exp.(X[2,:])
    W_net = exp.(Z .* γ) .* W
    A = exp.(X[3,:])
    C = [solve_consumption(r,α[i],W_net[i],A[i],σ,ψ) for i in eachindex(A)]
    @views H = exp.( X[1,:] .+ ψ .* log.(W_net) .- σ * ψ .* log.(C) )
    return (;α,W,A,C,H,W_net,Z)
end

# assume risk-aversion of 2 and frisch of 0.5
σ = 2.
ψ = 0.5
r = 0.05
γ = 0.2

N = 10_000

ψ_est = zeros(500)
ψ_ols = zeros(500)
for b in 1:500
    dat = simulate_data(σ,ψ,r,γ,N)
    #M = dat.A .< quantile(dat.A,0.75)
    M = dat.A .< median(dat.A)
    # construct instruments:
    Z = [ones(N) M dat.Z dat.Z .* M]
    X = [ones(N) M log.(dat.W_net) log.(dat.C)]
    # first stage:
    δ = inv(Z' * Z) * Z'*X
    Xh = Z * δ
    # second stage:
    β = inv(Xh' * Xh) * Xh' * log.(dat.H)
    ψ_est[b] = β[3]
    X = [ones(N) log.(dat.W_net) log.(dat.C)]
    β_ols = inv(X' * X) * X' * log.(dat.H)
    ψ_ols[b] = β_ols[2]
end
histogram(ψ_est,alpha=0.4)
histogram!(ψ_ols,alpha=0.4)
xlims!((-2,2))
# test:

# try this as well:
for b in 1:500
    dat = simulate_data(σ,ψ,r,γ,N)
    M1 = dat.A .< quantile(dat.A,0.25)
    M2 = dat.A .< quantile(dat.A,0.5)
    M3 = dat.A .< quantile(dat.A,0.75)
    # construct instruments:
    Z = [ones(N) M1 M2 M3 dat.Z dat.Z .* M1 dat.Z .* M2 dat.Z .* M3]
    X = [ones(N) M1 M2 M3 log.(dat.W_net) log.(dat.C)]
    # first stage:
    δ = inv(Z' * Z) * Z'*X
    Xh = Z * δ
    # second stage:
    β = inv(Xh' * Xh) * Xh' * log.(dat.H)
    ψ_est[b] = β[5]
end
histogram!(ψ_est,alpha=0.4)


