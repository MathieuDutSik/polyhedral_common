#!/usr/bin/env julia
#
# Multistart Newton-Raphson critical-point finder for the Lagrangian
# system, using the numerical evaluator. Better than monodromy when the
# polynomial is too large to symbolically materialize: just refine from
# many known/random start points.

using LinearAlgebra
using Printf
using JSON3
using Random

include(joinpath(@__DIR__, "custom_monodromy.jl"))

# Build a starting (qvec, μ) at the principal V^n form Q_0 = nI - J.
function start_principal(L::NumericalLType)
    n = L.n; d = L.d
    qvec = Vector{Float64}(undef, d)
    let kk = 1
        for i in 1:n, j in i:n
            qvec[kk] = (i == j) ? Float64(n) : -1.0
            kk += 1
        end
    end
    Q = Matrix{Float64}(undef, n, n); Q_from_qvec!(Q, qvec, n)
    g = zeros(Float64, d); grad_G_Delta!(g, L, qvec, Q)
    _, gd = det_gradient_and_hessian(Q, n)
    μ = g[1] / gd[1]
    return vcat(qvec, μ)
end

# Build a starting (qvec, μ) for Q = I (the cubic Z^n point).
function start_cubic(L::NumericalLType; scale::Float64=1.0)
    n = L.n; d = L.d
    qvec = zeros(Float64, d)
    let kk = 1
        for i in 1:n, j in i:n
            qvec[kk] = (i == j) ? scale : 0.0
            kk += 1
        end
    end
    Q = Matrix{Float64}(undef, n, n); Q_from_qvec!(Q, qvec, n)
    g = zeros(Float64, d); grad_G_Delta!(g, L, qvec, Q)
    _, gd = det_gradient_and_hessian(Q, n)
    μ = abs(gd[1]) > 1e-9 ? g[1] / gd[1] : 0.0
    return vcat(qvec, μ)
end

# Random PD starting point with diagonal dominance.
function start_random(L::NumericalLType, rng::AbstractRNG;
                      diag_scale::Float64=2.0)
    n = L.n; d = L.d
    qvec = zeros(Float64, d)
    let kk = 1
        for i in 1:n, j in i:n
            if i == j
                qvec[kk] = diag_scale + rand(rng) * diag_scale
            else
                qvec[kk] = -rand(rng) * 0.5
            end
            kk += 1
        end
    end
    Q = Matrix{Float64}(undef, n, n); Q_from_qvec!(Q, qvec, n)
    g = zeros(Float64, d); grad_G_Delta!(g, L, qvec, Q)
    _, gd = det_gradient_and_hessian(Q, n)
    μ = abs(gd[1]) > 1e-9 ? g[1] / gd[1] : 0.0
    return vcat(qvec, μ)
end

# Refine x to a critical point by Newton on the system F(x, p=0) = 0.
# Returns (x, converged).
function refine_to_critical(sys::QuantizerLagrangianSystem,
                           x_init::AbstractVector;
                           tol::Float64=1e-9,
                           max_iter::Int=80)
    d = sys.L.d
    p_zero = zeros(ComplexF64, d + 1)
    x = ComplexF64.(x_init)
    return newton_correct(sys, x, p_zero; tol=tol, max_iter=max_iter)
end

# Classify a critical point via the constrained Hessian on the det=1 slice.
# Returns (kind_string, eigenvalues).
function classify_point(L::NumericalLType, qvec::AbstractVector{Float64}, μ::Float64)
    n = L.n; d = L.d
    Q = Matrix{Float64}(undef, n, n); Q_from_qvec!(Q, qvec, n)
    H_G = zeros(Float64, d, d); hess_G_Delta!(H_G, L, qvec, Q)
    # Hessian of det via FD on grad_det
    H_d = zeros(Float64, d, d)
    Qbuf = similar(Q)
    h = FD_H
    for k in 1:d
        qp = copy(qvec); qp[k] += h
        qm = copy(qvec); qm[k] -= h
        Q_from_qvec!(Qbuf, qp, n); _, gdp = det_gradient_and_hessian(Qbuf, n)
        Q_from_qvec!(Qbuf, qm, n); _, gdm = det_gradient_and_hessian(Qbuf, n)
        for j in 1:d
            H_d[j, k] = (gdp[j] - gdm[j]) / (2 * h)
        end
    end
    H_L = H_G .- μ .* H_d
    _, gd = det_gradient_and_hessian(Q, n)
    nrm = norm(gd)
    nrm == 0 && return ("degenerate", Float64[])
    nvec = gd ./ nrm
    M = hcat(reshape(nvec, d, 1), Matrix{Float64}(LinearAlgebra.I, d, d))
    F = qr(M)
    Qbasis = Matrix(F.Q)
    T = Qbasis[:, 2:end]
    H_proj = T' * H_L * T
    H_proj = (H_proj .+ H_proj') ./ 2
    evs = sort(real.(eigvals(H_proj)))
    n_pos = count(>(1e-8), evs)
    n_neg = count(<(-1e-8), evs)
    n_zero = length(evs) - n_pos - n_neg
    kind = if n_neg == 0 && n_zero == 0
        "local MIN"
    elseif n_pos == 0 && n_zero == 0
        "local MAX"
    elseif n_zero > 0
        "degenerate (zero:$n_zero, +:$n_pos, -:$n_neg)"
    else
        "SADDLE (+:$n_pos, -:$n_neg)"
    end
    return (kind, evs)
end

# Main multistart routine.
function multistart_search(L::NumericalLType;
                           n_random::Int=50,
                           seed::Int=0,
                           timeout::Float64=600.0,
                           dedup_atol::Float64=1e-4)
    sys = QuantizerLagrangianSystem(L)
    n = L.n; d = L.d
    t0 = time()
    rng = MersenneTwister(seed == 0 ? abs(rand(Int)) : seed)
    starts = Vector{Vector{Float64}}()
    push!(starts, start_principal(L))
    push!(starts, start_cubic(L; scale=1.0))
    for _ in 1:n_random
        push!(starts, start_random(L, rng; diag_scale=2.0))
    end
    found = Vector{Vector{Float64}}()
    n_converged = 0
    for (i, x0) in enumerate(starts)
        time() - t0 > timeout && (println("timeout reached, stopping"); break)
        x_complex = ComplexF64.(x0)
        x_ref, ok = refine_to_critical(sys, x_complex)
        ok || continue
        # Take real if imaginary part is small
        max_imag = maximum(abs.(imag.(x_ref)))
        max_imag > 1e-5 && continue
        x_real = Float64[real(z) for z in x_ref]
        # Verify by re-evaluation
        u = zeros(ComplexF64, d + 1)
        HC.evaluate!(u, sys, ComplexF64.(x_real), zeros(ComplexF64, d + 1))
        norm(u) > 1e-5 && continue
        n_converged += 1
        # Dedup vs found
        is_new = true
        for x in found
            if norm(x_real .- x) < dedup_atol
                is_new = false; break
            end
        end
        is_new && push!(found, x_real)
    end
    elapsed = time() - t0
    println("multistart_search: $(length(starts)) starts, $(n_converged) Newton-converged, $(length(found)) distinct critical points, $(round(elapsed; digits=1))s")
    return found
end

# ----- driver ------------------------------------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 2 || error("Usage: multistart_critical.jl <json> <timeout_s> [n_random] [seed]")
    path = ARGS[1]
    timeout_s = parse(Float64, ARGS[2])
    n_random = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 100
    seed = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 0

    println("Loading L-type: $path")
    L = load_numerical_ltype(path)
    n = L.n; d = L.d
    println("n=$n, |star|=$(length(L.delaunay)), |T(Q)|=$(length(L.triangulation))")

    j = JSON3.read(read(path, String))
    G_Delta_value(L, _mat(j.test_point_Q, Float64); set_epsilon=true)

    found = multistart_search(L; n_random=n_random, seed=seed, timeout=timeout_s)

    println("\n=== Critical points ===")
    for (i, x) in enumerate(found)
        qvec = x[1:d]; μ = x[d+1]
        Q = Matrix{Float64}(undef, n, n); Q_from_qvec!(Q, qvec, n)
        Gd = G_Delta_value(L, Q)
        detQ = det(Q)
        G = Gd / detQ^((2n + 1) / n)
        is_pd = isposdef(Symmetric(Q))
        kind, evs = is_pd ? classify_point(L, qvec, μ) : ("non-PD", Float64[])
        @printf("  %3d: G=%+9.6f det=%+8.5f μ=%+9.6f pd=%s %s\n",
                i, G, detQ, μ, is_pd ? "yes" : "no ", kind)
        @printf("       qvec=%s\n",
                join([@sprintf("%+7.4f", v) for v in qvec], ","))
        if length(evs) <= 8
            @printf("       eigenvalues: %s\n",
                    join([@sprintf("%+8.4f", e) for e in evs], ", "))
        end
    end
end
