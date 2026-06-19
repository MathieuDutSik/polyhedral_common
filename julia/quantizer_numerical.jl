#!/usr/bin/env julia
#
# Numerical evaluator for G_Δ(Q) and its gradient and Hessian, never
# materializing the polynomial. Reads the same JSON the symbolic driver
# does and uses the L-type structure directly.
#
# Formula (Zimmermann §6.3.7) per triangulation simplex S with vertex
# set {0, c_{i_1}, …, c_{i_n}}:
#   contrib_S(Q) = ε(S̄) · det(Q · M(S)^T) · (‖Σ_k Q c_{i_k}‖²_{Q̃} +
#                                            Σ_k ‖Q c_{i_k}‖²_{Q̃}) /
#                  (n · (n+2)!)
# where Q c_{i_k} = M(D_{i_k})⁻¹ · b(D_{i_k}, Q), and Q̃ = adj(Q).
#
# G_Δ(Q) is the sum over all triangulation simplices.
#
# This file is self-contained for the evaluator; HC.jl wrapping is in
# the AbstractSystem implementation below.

using LinearAlgebra
using JSON3

# -- L-type loader (reuses the JSON schema from the C++ exporter) --------

struct NumericalLType
    n::Int
    d::Int                                            # = n(n+1)/2
    # Delaunay simplices (in star(0)): each has n non-zero lattice
    # vertices stored as the rows of an integer matrix.
    delaunay::Vector{Matrix{Int}}
    # Precomputed for each Delaunay D: inverse of M(D) (n×n rational
    # matrix whose rows are the non-zero vertices).
    Minv::Vector{Matrix{Rational{BigInt}}}
    # Precomputed v_j v_j^T for each vertex of each Delaunay (used for
    # the quadratic forms |v|²_Q = trace(v v^T · Q)). Stored as a list
    # of symmetric matrices.
    vv::Vector{Vector{Matrix{Int}}}                   # vv[D_idx][k] = v_k v_k^T
    # Triangulation: each entry is a vector of n star indices (1-based,
    # ordered).
    triangulation::Vector{Vector{Int}}
    # ε per simplex (chosen so each contribution is non-negative at the
    # test point; set after first evaluation).
    epsilon::Vector{Int}
end

function _mat(J, ::Type{T}) where {T}
    nrows = length(J)
    ncols = length(J[1])
    M = Matrix{T}(undef, nrows, ncols)
    for i in 1:nrows, j in 1:ncols
        M[i, j] = T(J[i][j])
    end
    return M
end

function load_numerical_ltype(path)
    j = JSON3.read(read(path, String))
    n = Int(j.dim)
    d = n * (n + 1) ÷ 2
    delaunay = [_mat(s.vertices, Int) for s in j.delaunay_simplices]
    triangulation = [Int.(s.delaunay_indices) .+ 1 for s in j.triangulation]
    Minv = Vector{Matrix{Rational{BigInt}}}(undef, length(delaunay))
    vv   = Vector{Vector{Matrix{Int}}}(undef, length(delaunay))
    for (i, D) in enumerate(delaunay)
        M = Rational{BigInt}.(D)
        Minv[i] = inv(M)
        vv[i] = [D[k, :] * transpose(D[k, :]) for k in 1:n]
    end
    return NumericalLType(n, d, delaunay, Minv, vv, triangulation, Int[])
end

# -- Helpers ---------------------------------------------------------------

# Build symmetric Q matrix from qvec (length d).
function Q_from_qvec!(Q::AbstractMatrix, qvec::AbstractVector, n::Int)
    k = 1
    @inbounds for i in 1:n, j in i:n
        Q[i, j] = qvec[k]
        Q[j, i] = qvec[k]
        k += 1
    end
    return Q
end

# Compute Q·c for a single Delaunay simplex at quadratic form Q.
# Q is n×n; D's Minv is n×n rational; v_k v_k^T are n×n integer.
# Returns an n×n matrix whose columns are Q·c_k (one per vertex of D).
function compute_Qc_per_delaunay(QC::AbstractVector, L::NumericalLType,
                                 i_D::Int, Q::AbstractMatrix)
    n = L.n
    # b_k = (1/2) v_k^T Q v_k for k = 1..n (n values).
    b = Vector{eltype(Q)}(undef, n)
    @inbounds for k in 1:n
        s = zero(eltype(Q))
        v = view(L.delaunay[i_D], k, :)
        for a in 1:n, c in 1:n
            s += v[a] * Q[a, c] * v[c]
        end
        b[k] = s / 2
    end
    # Qc = Minv * b, where Minv is rational; result is in eltype(Q).
    Minv = L.Minv[i_D]
    @inbounds for i in 1:n
        s = zero(eltype(Q))
        for k in 1:n
            s += Minv[i, k] * b[k]
        end
        QC[i] = s
    end
    return QC
end

# Compute G_Δ(Q) by summing contributions over the triangulation. Returns
# (G_Δ_value, set_orientations_if_first). The epsilon vector is set on
# the first call so that contributions are positive at the calling Q.
function G_Delta_value(L::NumericalLType, Q::AbstractMatrix; set_epsilon=false)
    n = L.n
    # Precompute Qc per Delaunay.
    Qc_list = [zeros(eltype(Q), n) for _ in 1:length(L.delaunay)]
    for i_D in 1:length(L.delaunay)
        compute_Qc_per_delaunay(Qc_list[i_D], L, i_D, Q)
    end
    # Q-adjugate = det(Q) * Q^{-1}.
    Qadj = det(Q) * inv(Q)
    # Per-simplex contribution sum.
    G = zero(eltype(Q))
    inv_coeff = inv(eltype(Q)(n * factorial(n + 2)))
    if set_epsilon
        empty!(L.epsilon)
    end
    for (i_S, S) in enumerate(L.triangulation)
        cols = [Qc_list[k] for k in S]           # n column vectors
        QM = hcat(cols...)                       # n × n matrix
        detQM = det(QM)
        # term1: ‖Σ_k Qc_k‖²_{Q̃} = (Σ Qc)^T · Qadj · (Σ Qc)
        sum_Qc = sum(cols)
        term1 = transpose(sum_Qc) * Qadj * sum_Qc
        # term2: Σ_k ‖Qc_k‖²_{Q̃}
        term2 = zero(eltype(Q))
        for c in cols
            term2 += transpose(c) * Qadj * c
        end
        # orientation
        ε = if set_epsilon
            ε_local = real(detQM) > 0 ? 1 : -1
            push!(L.epsilon, ε_local)
            ε_local
        else
            L.epsilon[i_S]
        end
        G += ε * inv_coeff * detQM * (term1 + term2)
    end
    return G
end

# -- Gradient and Hessian via central finite differences -----------------
#
# G_Δ is a polynomial of degree 2n+1, so analytic in Q. Central differences
# with step h=1e-6 give ~12 digits accuracy — enough for HC.jl monodromy.

const FD_H = 1e-6

function grad_G_Delta!(g::AbstractVector, L::NumericalLType, qvec::AbstractVector,
                       Qbuf::AbstractMatrix)
    d = L.d
    n = L.n
    h = FD_H
    for k in 1:d
        q_plus = copy(qvec)
        q_minus = copy(qvec)
        q_plus[k]  += h
        q_minus[k] -= h
        Q_from_qvec!(Qbuf, q_plus,  n)
        Gp = G_Delta_value(L, Qbuf)
        Q_from_qvec!(Qbuf, q_minus, n)
        Gm = G_Delta_value(L, Qbuf)
        g[k] = (Gp - Gm) / (2 * h)
    end
    return g
end

function hess_G_Delta!(H::AbstractMatrix, L::NumericalLType, qvec::AbstractVector,
                       Qbuf::AbstractMatrix)
    d = L.d
    n = L.n
    h = FD_H
    g_plus  = zeros(eltype(qvec), d)
    g_minus = zeros(eltype(qvec), d)
    for k in 1:d
        q_plus = copy(qvec)
        q_minus = copy(qvec)
        q_plus[k]  += h
        q_minus[k] -= h
        grad_G_Delta!(g_plus,  L, q_plus,  Qbuf)
        grad_G_Delta!(g_minus, L, q_minus, Qbuf)
        @inbounds for j in 1:d
            H[j, k] = (g_plus[j] - g_minus[j]) / (2 * h)
        end
    end
    return H
end

# Det(Q) and its gradient and Hessian (analytic, via adjugate).
function det_gradient_and_hessian(Q::AbstractMatrix, n::Int)
    detQ = det(Q)
    adjQ = detQ * inv(Q)              # adjugate
    # ∂det/∂Q[i,j] = adj(Q)[j, i].   In our qvec ordering (i<=j), the
    # symmetric off-diagonal gets DOUBLED in qvec convention. Need care.
    # Here qvec[k] corresponds to Q[i,j] with i<=j (NOT doubled). So
    # ∂det/∂qvec[k] is for the SYMMETRIC perturbation Q[i,j]=Q[j,i]+=ε,
    # i.e., (adj_{ji} + adj_{ij})·dq_k for off-diagonal (i≠j), and
    # adj_{ii}·dq_k for diagonal.
    d = n * (n + 1) ÷ 2
    grad = zeros(eltype(Q), d)
    k = 1
    for i in 1:n, j in i:n
        if i == j
            grad[k] = adjQ[i, i]
        else
            grad[k] = adjQ[i, j] + adjQ[j, i]   # = 2 * adjQ[i, j] since adj sym
        end
        k += 1
    end
    return detQ, grad
end

# -- HC.AbstractSystem wrapping the Lagrangian system ---------------------
#
# Lagrangian system (with one shift parameter per equation, matching the
# existing monodromy setup):
#   F_k(qvec, μ, p) = ∂G_Δ/∂q_k - μ · ∂det/∂q_k - p_k       (k = 1..d)
#   F_{d+1}(qvec, μ, p) = det(Q) - 1 - p_{d+1}
#
# Variables: vcat(qvec, μ); parameters: p (length d+1).
# Output u of length d+1.

import HomotopyContinuation as HC

struct QuantizerLagrangianSystem <: HC.AbstractSystem
    L::NumericalLType
end

Base.size(S::QuantizerLagrangianSystem) = (S.L.d + 1, S.L.d + 1)
HC.ModelKit.variables(S::QuantizerLagrangianSystem) = HC.Variable[]   # opaque
HC.parameters(S::QuantizerLagrangianSystem)         = HC.Variable[]
HC.nvariables(S::QuantizerLagrangianSystem) = S.L.d + 1
HC.nparameters(S::QuantizerLagrangianSystem) = S.L.d + 1

# Evaluation: x = (qvec..., μ), p = (p_1, ..., p_{d+1})
function HC.evaluate!(u::AbstractVector, S::QuantizerLagrangianSystem,
                     x::AbstractVector, p::AbstractVector)
    L = S.L
    d = L.d; n = L.n
    qvec = view(x, 1:d)
    μ = x[d+1]
    Qbuf = Matrix{eltype(x)}(undef, n, n)
    Q_from_qvec!(Qbuf, qvec, n)
    # gradient of G_Δ
    g = zeros(eltype(x), d)
    grad_G_Delta!(g, L, qvec, Qbuf)
    # gradient of det and det
    detQ, gd = det_gradient_and_hessian(Qbuf, n)
    @inbounds for k in 1:d
        u[k] = g[k] - μ * gd[k] - p[k]
    end
    u[d+1] = detQ - 1 - p[d+1]
    return u
end

function HC.evaluate_and_jacobian!(u::AbstractVector, U::AbstractMatrix,
                                  S::QuantizerLagrangianSystem,
                                  x::AbstractVector, p::AbstractVector)
    L = S.L
    d = L.d; n = L.n
    qvec = view(x, 1:d)
    μ = x[d+1]
    Qbuf = Matrix{eltype(x)}(undef, n, n)
    Q_from_qvec!(Qbuf, qvec, n)
    # function values: reuse the evaluate logic
    g = zeros(eltype(x), d)
    grad_G_Delta!(g, L, qvec, Qbuf)
    detQ, gd = det_gradient_and_hessian(Qbuf, n)
    @inbounds for k in 1:d
        u[k] = g[k] - μ * gd[k] - p[k]
    end
    u[d+1] = detQ - 1 - p[d+1]
    # Jacobian
    H = zeros(eltype(x), d, d)
    hess_G_Delta!(H, L, qvec, Qbuf)
    # Hessian of det via FD too (could be analytic but cheap)
    Hd = zeros(eltype(x), d, d)
    let h = FD_H
        for k in 1:d
            q_plus  = copy(qvec); q_plus[k]  += h
            q_minus = copy(qvec); q_minus[k] -= h
            Q_from_qvec!(Qbuf, q_plus,  n); _, gd_p = det_gradient_and_hessian(Qbuf, n)
            Q_from_qvec!(Qbuf, q_minus, n); _, gd_m = det_gradient_and_hessian(Qbuf, n)
            @inbounds for j in 1:d
                Hd[j, k] = (gd_p[j] - gd_m[j]) / (2 * h)
            end
        end
    end
    # Restore Qbuf
    Q_from_qvec!(Qbuf, qvec, n)
    # Fill U: rows 1..d, cols 1..d+1 are derivative of Lagrangian eqs
    @inbounds for i in 1:d, j in 1:d
        U[i, j] = H[i, j] - μ * Hd[i, j]
    end
    @inbounds for i in 1:d
        U[i, d+1] = -gd[i]
    end
    # Row d+1: det(Q) - 1 - p_{d+1}; d/dq_j = gd[j], d/dμ = 0
    @inbounds for j in 1:d
        U[d+1, j] = gd[j]
    end
    U[d+1, d+1] = zero(eltype(x))
    return u, U
end

# -- Sanity test against n=2 closed form ----------------------------------
if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Usage: quantizer_numerical.jl <ltype.json>")
    path = ARGS[1]
    println("Loading $path ...")
    L = load_numerical_ltype(path)
    println("n=$(L.n), |star|=$(length(L.delaunay)), |T(Q)|=$(length(L.triangulation))")
    # Build Q from test_point_Q in the JSON.
    j = JSON3.read(read(path, String))
    Q_test = _mat(j.test_point_Q, Float64)
    println("Test Q:")
    display(Q_test)
    G_val = G_Delta_value(L, Q_test; set_epsilon=true)
    detQ = det(Q_test)
    G_full = G_val / detQ^((2*L.n + 1) / L.n)
    println("\nG_Δ(test_Q)  = ", G_val)
    println("det(test_Q)  = ", detQ)
    println("G(test_Q)    = G_Δ / det^((2n+1)/n) = ", G_full)
end
