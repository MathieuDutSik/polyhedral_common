#!/usr/bin/env julia
# Run monodromy on the NUMERICAL Lagrangian system.
#
# Args: <ltype.json> <timeout_seconds>
import HomotopyContinuation as HC
using LinearAlgebra
using JSON3
using Printf

include(joinpath(@__DIR__, "quantizer_numerical.jl"))

length(ARGS) >= 2 || error("Usage: test_numerical_monodromy.jl <json> <timeout_seconds>")
path = ARGS[1]
timeout_s = parse(Float64, ARGS[2])

println("Loading L-type from $path ...")
L = load_numerical_ltype(path)
println("n=$(L.n), |star|=$(length(L.delaunay)), |T(Q)|=$(length(L.triangulation))")

# Initialize epsilons by evaluating at the JSON's test_point_Q.
j = JSON3.read(read(path, String))
test_Q = _mat(j.test_point_Q, Float64)
G_Delta_value(L, test_Q; set_epsilon=true)
println("Epsilons set from test_Q")

# Build start solution: principal V^n form Q_0 = nI - J (after row sign
# convention). qvec ordering is (q11, q12, q13, ..., q_{nn}).
n = L.n
d = L.d
qvec_start = Vector{Float64}(undef, d)
let kk = 1
    for i in 1:n, j in i:n
        qvec_start[kk] = (i == j) ? Float64(n) : -1.0
        kk += 1
    end
end

# Numeric Q from qvec
Q_start = Matrix{Float64}(undef, n, n)
Q_from_qvec!(Q_start, qvec_start, n)
println("Q_start = "); display(Q_start)

# Compute μ_start such that ∂G_Δ/∂q_1 = μ · ∂det/∂q_1 at qvec_start
g = zeros(Float64, d)
grad_G_Delta!(g, L, qvec_start, Q_start)
detQ_start, gd = det_gradient_and_hessian(Q_start, n)
μ_start = g[1] / gd[1]
println("μ_start = ", μ_start)

x_start = vcat(qvec_start, μ_start)

# Parameter values such that F(x_start, p) = 0:
#   F_k = g_k - μ·gd_k - p_k → p_k = g_k - μ·gd_k  (k=1..d)
#   F_{d+1} = det - 1 - p_{d+1}      → p_{d+1} = det - 1
p_start = zeros(Float64, d + 1)
@inbounds for k in 1:d
    p_start[k] = g[k] - μ_start * gd[k]
end
p_start[d+1] = detQ_start - 1
println("p_start residuals (should be ~0 with the right μ): ", maximum(abs.(p_start)))

# Wrap as HC.AbstractSystem
sys = QuantizerLagrangianSystem(L)
println("System sizes: nvars=", HC.nvariables(sys), " neq=", size(sys, 1),
        " npar=", HC.nparameters(sys))

# Verify HC.jl can evaluate
u = zeros(ComplexF64, d + 1)
U = zeros(ComplexF64, d + 1, d + 1)
x_c = ComplexF64.(x_start)
p_c = ComplexF64.(p_start)
HC.evaluate_and_jacobian!(u, U, sys, x_c, p_c)
println("‖u‖ at start (should be ≪ 1):  ", maximum(abs.(u)))
println("‖J‖ at start: ", maximum(abs.(U)))

# Try plain parameter homotopy first: HC.solve with start and target params.
# This may avoid the symbolic conversion that monodromy_solve forces.
p_target = ComplexF64.(zeros(d + 1))   # at p=0, we recover the target system
println("Trying parameter homotopy: start_params -> p=0")
result = HC.solve(sys, [ComplexF64.(x_start)];
                  start_parameters=ComplexF64.(p_start),
                  target_parameters=p_target,
                  show_progress=true,
                  threading=false)
println("\nMonodromy done. Found ", length(HC.solutions(result)), " complex solutions")
raw_sols = HC.solutions(result)
real_sols = [Float64[real(z) for z in cs] for cs in raw_sols
             if all(abs(imag(z)) < 1e-5 for z in cs)]
println("Real solutions: ", length(real_sols))
for (i, s) in enumerate(real_sols)
    qvals = s[1:d]
    μ_v = s[d+1]
    Q = Matrix{Float64}(undef, n, n); Q_from_qvec!(Q, qvals, n)
    G_Δ = G_Delta_value(L, Q)
    detQ = det(Q)
    G = G_Δ / detQ^((2n + 1) / n)
    posdef = isposdef(Q)
    @printf("  %3d: G=%+9.6f det=%+8.5f μ=%+9.6f pd=%s qvec=%s\n",
            i, G, detQ, μ_v, posdef ? "yes" : "no ",
            join([@sprintf("%+7.4f", v) for v in qvals], ","))
end
