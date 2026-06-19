#!/usr/bin/env julia
#
# Run numerical monodromy on a given L-type JSON. Reports critical points,
# classification, orbits.
#
# Args: <ltype.json> <timeout_seconds> [seed]
using LinearAlgebra
using Printf
using JSON3
using Random

include(joinpath(@__DIR__, "custom_monodromy.jl"))

length(ARGS) >= 2 || error("Usage: run_numerical.jl <json> <timeout_seconds> [seed]")
path = ARGS[1]
timeout_s = parse(Float64, ARGS[2])
seed = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 0

println("Loading L-type from $path ...")
L = load_numerical_ltype(path)
n = L.n; d = L.d
println("n=$n, |star|=$(length(L.delaunay)), |T(Q)|=$(length(L.triangulation))")

# Initialize epsilons from JSON test_point_Q
j = JSON3.read(read(path, String))
test_Q = _mat(j.test_point_Q, Float64)
G_Delta_value(L, test_Q; set_epsilon=true)
println("Epsilons set from test_Q (det(test_Q)=$(det(test_Q)))")

# Start point: principal V^n form Q_0 = nI - J
qvec_start = Vector{Float64}(undef, d)
let kk = 1
    for i in 1:n, jj in i:n
        qvec_start[kk] = (i == jj) ? Float64(n) : -1.0
        kk += 1
    end
end
Q_start = Matrix{Float64}(undef, n, n); Q_from_qvec!(Q_start, qvec_start, n)
println("Q_start (the V^n principal form):")
display(Q_start)

# Find μ from ∂G/∂q_1 = μ · ∂det/∂q_1
g = zeros(Float64, d); grad_G_Delta!(g, L, qvec_start, Q_start)
detQ_start, gd = det_gradient_and_hessian(Q_start, n)
μ_start = g[1] / gd[1]
x_start = vcat(qvec_start, μ_start)
# Corresponding shift parameter such that F = 0 at start
p_start = zeros(Float64, d + 1)
@inbounds for k in 1:d
    p_start[k] = g[k] - μ_start * gd[k]
end
p_start[d+1] = detQ_start - 1
println("μ_start = $μ_start, max|p_start| = $(maximum(abs.(p_start)))")

# Sanity: F at (x_start, p_start) ≈ 0
sys = QuantizerLagrangianSystem(L)
u0 = zeros(ComplexF64, d + 1)
HC.evaluate!(u0, sys, ComplexF64.(x_start), ComplexF64.(p_start))
println("‖F(x_start, p_start)‖ = $(maximum(abs.(u0)))")

# Run custom monodromy
println("\nStarting custom monodromy (timeout=$(timeout_s)s, seed=$seed) ...")
sols = custom_monodromy(sys, ComplexF64.(x_start), ComplexF64.(p_start);
                        timeout=timeout_s,
                        seed=seed,
                        verbose=true)

# Filter to real (after walking back to p=p_start, F=0 ⟹ x is a critical point)
real_sols = Vector{Vector{Float64}}()
for cs in sols
    if all(abs(imag(z)) < 1e-5 for z in cs)
        push!(real_sols, Float64[real(z) for z in cs])
    end
end
println("Real solutions: $(length(real_sols)) of $(length(sols)) total")

# Classify
expected_G_denom = (2n + 1) / n
println("\nReal critical points:")
for (i, sol) in enumerate(real_sols)
    qvals = sol[1:d]; μ_v = sol[d+1]
    Q = Matrix{Float64}(undef, n, n); Q_from_qvec!(Q, qvals, n)
    Gd = G_Delta_value(L, Q)
    detQ = det(Q)
    G = Gd / detQ^expected_G_denom
    is_pd = isposdef(Symmetric(Q))
    @printf("  %3d: G=%+9.6f det=%+8.5f μ=%+9.6f pd=%s qvec=%s\n",
            i, G, detQ, μ_v, is_pd ? "yes" : "no ",
            join([@sprintf("%+7.4f", v) for v in qvals], ","))
end
