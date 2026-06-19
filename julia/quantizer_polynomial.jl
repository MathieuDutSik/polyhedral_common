#!/usr/bin/env julia
#
# Reads an L-type description in JSON (see example_n2_V2.json) and builds the
# Zimmermann quantizer polynomial G_Δ(Q) symbolically.
#
# Formula (Zimmermann 2017 thesis, Theorem 6.3.7):
#
#   G_Δ(Q) = (1 / (n * (n+2)!)) *
#            Σ_{S ∈ T(Q), dim S = n}
#              ε(S̄) * det(Q · M(S)^T) *
#              ( ‖Σ_k Q·c_{i_k}‖²_{Q̃}  +  Σ_k ‖Q·c_{i_k}‖²_{Q̃} )
#
# where:
#   - n              = lattice dimension
#   - T(Q)           = pulling triangulation of DV(0), pyramidal with apex 0
#   - S              = {0, c_{i_1}, …, c_{i_n}}, c_{i_k} = circumcenter of a
#                       Delaunay n-simplex D_{i_k} ∈ star(0, 𝓓)
#   - M(S)           = matrix whose rows are the non-zero vertices of S in DV(0)
#   - Q̃ = adj(Q)    = det(Q) · Q⁻¹
#   - Q · c_{i_k}    = M(D_{i_k})⁻¹ · b(D_{i_k})   (Lemma 6.3.5;
#                       M(D) is integer, b is linear in Q-entries)
#
# The full quantizer constant is G(Q) = G_Δ(Q) / det(Q)^((2n+1)/n).
# By homogeneity, critical points may be sought on the slice det(Q) = 1.

using DynamicPolynomials
using MultivariatePolynomials
import HomotopyContinuation as HC
import LinearAlgebra
using JSON3
using LinearAlgebra
using Printf

# --- Polynomial helpers --------------------------------------------------------

"""
Return symmetric matrix `Q` whose distinct entries are independent polynomial
variables, along with the flat vector of those variables. Entries q[i,j] with
i ≤ j are independent; q[j,i] = q[i,j].
"""
function symmetric_Q(n::Int)
    DynamicPolynomials.@polyvar qvec[1:(n*(n+1) ÷ 2)]
    # Promote entries to Polynomial so cofactors/products land in a uniform type.
    P = typeof(qvec[1] * qvec[1] + qvec[1])
    Q = Matrix{P}(undef, n, n)
    k = 0
    for i in 1:n, j in i:n
        k += 1
        Q[i, j] = convert(P, qvec[k])
        Q[j, i] = Q[i, j]
    end
    return Q, collect(qvec)
end

"""
adj(Q): adjugate (classical adjoint) of an n×n matrix, computed by cofactor
expansion. Each entry is a polynomial of degree n-1 in the entries of Q.
"""
function adjugate(Q::AbstractMatrix)
    n = size(Q, 1)
    @assert size(Q, 2) == n
    function minor_det(i, j)
        rows = setdiff(1:n, i)
        cols = setdiff(1:n, j)
        return det_via_laplace(Q[rows, cols])
    end
    A = Matrix{eltype(Q)}(undef, n, n)
    for i in 1:n, j in 1:n
        A[i, j] = ((-1)^(i + j)) * minor_det(j, i)
    end
    return A
end

"""
Convert a polynomial that has no remaining variables (e.g. fully substituted)
into a plain number.
"""
function _to_number(p)
    # MultivariatePolynomials returns a "polynomial" with one constant term;
    # extract by evaluating against an empty assignment.
    pterm = polynomial(p)
    isapprox(pterm, zero(pterm)) && return 0 // 1
    # Sum the coefficients (only constant monomials survive after subs).
    return sum(MultivariatePolynomials.coefficient(t) for t in terms(pterm))
end

"""
det() for small symbolic matrices via Laplace expansion (DynamicPolynomials'
generic det relies on LU and complains about polynomial entries).
"""
function det_via_laplace(M::AbstractMatrix)
    n = size(M, 1)
    @assert size(M, 2) == n
    n == 0 && return one(eltype(M))
    n == 1 && return M[1, 1]
    if n == 2
        return M[1, 1] * M[2, 2] - M[1, 2] * M[2, 1]
    end
    s = zero(eltype(M))
    for j in 1:n
        cof = ((-1)^(1 + j)) * det_via_laplace(M[2:end, [1:j-1; j+1:n]])
        s += M[1, j] * cof
    end
    return s
end

# --- L-type description --------------------------------------------------------

struct LType{T}
    n::Int
    Q::Matrix{T}            # symbolic symmetric matrix
    qvec::Vector            # flat vector of distinct entries (length n(n+1)/2)
    delaunay::Vector{Matrix{Int}}        # each: n×n matrix, rows = non-zero vertices
    triangulation::Vector{Vector{Int}}   # each: list of n delaunay-simplex indices
    cone_facets_tspace::Vector{Vector{Rational{Int}}}
    tspace_basis::Vector{Matrix{Rational{Int}}}
    test_Q::Matrix{Rational{Int}}
end

"""
Decode a JSON 2-D array of numbers into a Matrix{T} (rows-of-rows).
"""
function _mat(J, ::Type{T}) where {T}
    nrows = length(J)
    ncols = length(J[1])
    M = Matrix{T}(undef, nrows, ncols)
    for i in 1:nrows, j in 1:ncols
        M[i, j] = T(J[i][j])
    end
    return M
end

function load_ltype(path::AbstractString)
    j = JSON3.read(read(path, String))
    n = Int(j.dim)
    Q, qvec = symmetric_Q(n)
    delaunay = [_mat(s.vertices, Int) for s in j.delaunay_simplices]
    triangulation = [Int.(s.delaunay_indices) for s in j.triangulation]
    cone_facets = [Rational{Int}.(collect(v)) for v in j.cone_facets_tspace]
    tspace_basis = [_mat(M, Rational{Int}) for M in j.tspace.basis]
    test_Q = _mat(j.test_point_Q, Rational{Int})
    return LType(n, Q, qvec, delaunay, triangulation, cone_facets, tspace_basis, test_Q)
end

# --- G_Δ construction ----------------------------------------------------------

"""
Compute the vector Q·c for the given Delaunay simplex D, in the symbolic ring.

For D with non-zero vertices s_1, ..., s_n (rows of `D_rows`), let
M(D) be the n×n integer matrix with rows s_k^T and b(D) the column vector
b_k = (1/2)·s_k^T Q s_k. Then Q·c = M(D)⁻¹ · b. Because M(D) has integer
entries, M(D)⁻¹ is rational of constant determinant; Qc has entries that are
polynomials of degree 1 in the entries of Q divided by det(M(D)).

We return the pair (Qc_numerator_polynomial_vector, det_M_int) so callers can
keep arithmetic exact.
"""
function compute_Qc(D_rows::Matrix{Int}, Q::AbstractMatrix, n::Int)
    M = Rational{BigInt}.(D_rows)                        # n × n int matrix
    detM = det(M)
    @assert detM != 0 "Delaunay simplex is degenerate"
    Minv = inv(M)                                        # rational entries
    half = 1 // 2
    b = [half * sum(D_rows[k, i] * D_rows[k, j] * Q[i, j] for i in 1:n, j in 1:n)
         for k in 1:n]
    Qc = Minv * b
    return Qc, detM
end

"""
Build the symbolic G_Δ polynomial. Returns (G_Δ, ε_array) where ε_array[s] ∈
{+1, -1} is the orientation sign chosen so that det(Q·M(S)^T) > 0 at the test
point for every triangulation simplex S.
"""
function build_G_Delta(L::LType)
    n = L.n
    Q = L.Q
    Qtilde = adjugate(Q)

    # Precompute Q·c for every Delaunay simplex (vector of polynomials)
    QcList = Vector{Tuple{Vector,Rational{BigInt}}}(undef, length(L.delaunay))
    for (i, D) in enumerate(L.delaunay)
        QcList[i] = compute_Qc(D, Q, n)
    end

    # Evaluator at the test point, used only to pick ε per simplex
    test_assignments = Dict{Any,Rational{Int}}()
    qidx = 1
    for i in 1:n, j in i:n
        test_assignments[L.qvec[qidx]] = L.test_Q[i, j]
        qidx += 1
    end
    eval_at_test(p) = subs(p, [k => v for (k, v) in test_assignments]...)

    Sn = 1
    for k in 2:(n + 2)
        Sn *= k
    end
    inv_coeff = 1 // (n * Sn)   # 1 / (n · (n+2)!)

    G = zero(eltype(Q))
    ε_list = Int[]
    for S in L.triangulation
        @assert length(S) == n "Triangulation simplex must have $n indices, got $(length(S))"
        cols = [QcList[k+1][1] for k in S]   # 0-indexed JSON → 1-indexed Julia
        # Q · M(S)^T has columns Q·c_{i_k}. Since `Minv` already has rational
        # entries, the polynomial det(QM) is exactly det(Q)·det(M(S)).
        QM = reduce(hcat, cols)
        detQM = det_via_laplace(QM)          # polynomial of degree n in Q
        # Orientation: choose ε so that ε · detQM > 0 at test point.
        val_at_test = Float64(_to_number(eval_at_test(detQM)))
        ε = val_at_test > 0 ? +1 : -1
        push!(ε_list, ε)

        sum_Qc = sum(cols[k] for k in 1:n)
        term1 = sum(sum_Qc[i] * Qtilde[i, j] * sum_Qc[j] for i in 1:n, j in 1:n)
        term2 = zero(eltype(Q))
        for k in 1:n
            v = cols[k]
            term2 += sum(v[i] * Qtilde[i, j] * v[j] for i in 1:n, j in 1:n)
        end
        G += (ε * inv_coeff) * detQM * (term1 + term2)
    end
    return G, ε_list
end

# --- main ---------------------------------------------------------------------

"""
Evaluate a symbolic polynomial at a numeric assignment of its variables.
"""
function eval_poly_at(p, L::LType)
    assignments = Dict{Any, Rational{BigInt}}()
    qidx = 1
    for i in 1:L.n, j in i:L.n
        assignments[L.qvec[qidx]] = Rational{BigInt}(L.test_Q[i, j])
        qidx += 1
    end
    return _to_number(subs(p, [k => v for (k, v) in assignments]...))
end

# --- Critical-point search via HomotopyContinuation -------------------------

"""
Convert a DynamicPolynomials polynomial into a HomotopyContinuation expression,
returning a tuple (hc_polynomial, hc_variables_in_qvec_order).
"""
function to_hc(p, L::LType)
    nQ = length(L.qvec)
    hc_vars = HC.Variable[]
    name_map = Dict{Any,Any}()
    for (k, q) in enumerate(L.qvec)
        v = HC.Variable(Symbol("Q$k"))
        push!(hc_vars, v)
        name_map[q] = v
    end
    # Walk terms of `p` and rebuild in HC expression space.
    hc_p = HC.Expression(0)
    for t in terms(polynomial(p))
        c = MultivariatePolynomials.coefficient(t)
        mon = MultivariatePolynomials.monomial(t)
        expvec = MultivariatePolynomials.exponents(mon)
        varsvec = MultivariatePolynomials.variables(mon)
        # Build the corresponding HC monomial.
        hc_mon = HC.Expression(Float64(c))
        for (var, e) in zip(varsvec, expvec)
            e == 0 && continue
            hc_v = name_map[var]
            hc_mon *= hc_v^Int(e)
        end
        hc_p += hc_mon
    end
    return hc_p, hc_vars
end

"""
Symbolic det(Q) for the symbolic matrix.
"""
function det_Q_poly(Q)
    n = size(Q, 1)
    return det_via_laplace(Q)
end

"""
Solve the Lagrangian system
    ∂G_Δ/∂q_k - μ · ∂det(Q)/∂q_k = 0,   det(Q) = 1
and return the list of real solutions (each a NamedTuple of q-values).
"""
function find_critical_points(L::LType, G; method=:polyhedral)
    n = L.n
    detQ = det_via_laplace(L.Q)

    # Convert G_Δ and det(Q) to HC variables.
    G_hc, hc_qvars = to_hc(G, L)
    detQ_hc, _    = to_hc(detQ, L)

    HC.@var mu
    eqs = HC.Expression[]
    for k in 1:length(hc_qvars)
        gk  = HC.differentiate(G_hc,    hc_qvars[k])
        dk  = HC.differentiate(detQ_hc, hc_qvars[k])
        push!(eqs, gk - mu * dk)
    end
    push!(eqs, detQ_hc - 1)

    if method === :polyhedral
        F = HC.System(eqs; variables = vcat(hc_qvars, [mu]))
        println("Solving Lagrangian system: ",
                length(eqs), " polynomial equations in ", length(eqs), " unknowns")
        result = HC.solve(F; show_progress = false)
        println("Path-tracking: total=", length(HC.solutions(result)),
                "  real=", HC.nreal(result))
        return result, hc_qvars
    elseif method === :monodromy
        # Richer parameter homotopy: introduce one shift parameter per
        # Lagrangian equation. F_i(x) - p_i = 0 (with the det equation
        # using det(Q) - p_det = 1 - p_det shift). At p = 0 we recover
        # the target system. The multi-dimensional parameter space gives
        # much richer monodromy than a single det parameter.
        d = length(hc_qvars)
        n_eqs = length(eqs)
        HC.@var p_param[1:n_eqs]
        eqs_param = [eqs[i] - p_param[i] for i in 1:n_eqs]
        F_param = HC.System(eqs_param;
                            variables = vcat(hc_qvars, [mu]),
                            parameters = collect(p_param))
        println("Monodromy: System has ", length(eqs_param),
                " equations in ", d+1, " variables, ", n_eqs, " parameters")
        # Starting solution: principal V^n form Q_0 = (n+1)I − J.
        x_start_vec = zeros(Float64, d)
        kk = 1
        for i in 1:n, j in i:n
            x_start_vec[kk] = (i == j) ? Float64(n) : -1.0
            kk += 1
        end
        G_grad_first = HC.differentiate(G_hc, hc_qvars[1])
        det_grad_first = HC.differentiate(detQ_hc, hc_qvars[1])
        assg = Dict(hc_qvars[k] => x_start_vec[k] for k in 1:d)
        gnum = Float64(HC.evaluate(G_grad_first, assg))
        dnum = Float64(HC.evaluate(det_grad_first, assg))
        mu_start = abs(dnum) > 1e-12 ? gnum / dnum :
                   error("∂det/∂q_1 vanishes at start point; pick different start")
        x_start = vcat(x_start_vec, mu_start)
        # Parameters at start: p_i = eqs[i] evaluated at (x_start), so that
        # eqs_param[i] = eqs[i] - p_i = 0 at the start.
        full_assg = merge(assg, Dict(mu => mu_start))
        p_start_vec = [Float64(HC.evaluate(eqs[i], full_assg)) for i in 1:n_eqs]
        println("Start: μ = ", mu_start, ", max |p_i| = ", maximum(abs.(p_start_vec)))
        # Wall-clock timeout (seconds) + bound on solution count, to avoid
        # explosive ghost-solution growth at larger n.
        budget_seconds = parse(Float64, get(ENV, "QUANT_MONODROMY_TIMEOUT", "300"))
        max_sols       = parse(Int,     get(ENV, "QUANT_MONODROMY_MAX_SOLS", "2000"))
        opts = HC.MonodromyOptions(; max_loops_no_progress = 6,
                                     timeout = budget_seconds,
                                     target_solutions_count = max_sols)
        seed_env = parse(UInt32, get(ENV, "QUANT_MONODROMY_SEED", "0"))
        kw_seed = seed_env == 0 ? (;) : (; seed = seed_env)
        result = HC.monodromy_solve(F_param,
                                    [ComplexF64.(x_start)],
                                    ComplexF64.(p_start_vec);
                                    options = opts,
                                    show_progress = true,
                                    threading = false,
                                    kw_seed...)
        println("Monodromy found: ", length(HC.solutions(result)), " solutions")
        return result, hc_qvars
    end
    error("Unknown method: $method")
end

if abspath(PROGRAM_FILE) == @__FILE__
    isempty(ARGS) && error("Usage: julia quantizer_polynomial.jl <ltype.json>")
    path = ARGS[1]
    println("Loading L-type from: $path")
    L = load_ltype(path)
    println("Dimension n = $(L.n)")
    println("# Delaunay simplices in star(0): $(length(L.delaunay))")
    println("# triangulation simplices T(Q):  $(length(L.triangulation))")
    println("Building G_Δ symbolically...")
    G, ε = build_G_Delta(L)
    println("Orientations chosen: $ε")
    println()
    println("G_Δ = $G")
    println()

    # Compare against the explicit closed-form polynomial from eq (6.3) when n=2.
    if L.n == 2
        # eq (6.3) gives G(Q) = num/(24·det^{3/2}) and G_Δ = det · num / 24.
        # In Julia variable naming: qvec[1]=q11, qvec[2]=q12, qvec[3]=q22.
        q11, q12, q22 = L.qvec[1], L.qvec[2], L.qvec[3]
        num_zim = q11^2*q22 + q11*q22^2 - 2*q11*q12^2 - 2*q22*q12^2 - 2*q12^3
        det_Q = q11*q22 - q12^2
        G_zim = det_Q * num_zim * (1 // 24)
        diff = G - G_zim
        # Use the evaluation-at-test-point heuristic to check if `diff` is zero
        # (sufficient for a sanity-check; an exact symbolic equality would need
        # polynomial normalization).
        diff_val_test = Float64(eval_poly_at(diff, L))
        println("G_Δ(test) computed = $(Float64(eval_poly_at(G, L)))")
        println("G_Δ(test) expected = $(Float64(eval_poly_at(G_zim, L)))")
        println("|G - G_zim|(test) = $(abs(diff_val_test))")
        # Try a second random-ish point in the cone to verify equality more robustly.
        let q11v=3//1, q22v=5//1, q12v=-1//1
            assg = Dict(L.qvec[1]=>q11v, L.qvec[2]=>q12v, L.qvec[3]=>q22v)
            gv  = Float64(_to_number(subs(G,     [k=>v for (k,v) in assg]...)))
            gzv = Float64(_to_number(subs(G_zim, [k=>v for (k,v) in assg]...)))
            println("G_Δ at (q11=3,q22=5,q12=-1): computed=$gv expected=$gzv")
        end
    end

    println()
    println("=== Critical-point search via HomotopyContinuation ===")
    method = length(ARGS) >= 2 && ARGS[2] == "monodromy" ? :monodromy : :polyhedral
    println("Method: $method")
    result, hc_qvars = find_critical_points(L, G; method=method)

    # Cone-interior filter and quantizer-constant evaluation.
    facets = L.cone_facets_tspace
    n = L.n
    expG_denom_exp = (2*n + 1) // n   # G = G_Δ / det(Q)^{(2n+1)/n}

    # Helper: rebuild symmetric Q from qvec.
    function Q_from_qvec(qvals)
        Qm = zeros(Float64, n, n)
        kk = 1
        for i in 1:n, j in i:n
            Qm[i, j] = qvals[kk]; Qm[j, i] = Qm[i, j]
            kk += 1
        end
        return Qm
    end
    function qvec_from_Q(Qm)
        v = zeros(Float64, n*(n+1)÷2)
        kk = 1
        for i in 1:n, j in i:n
            v[kk] = Qm[i, j]; kk += 1
        end
        return v
    end

    # Eval G_Δ at a numeric qvec
    function eval_G_Δ(qvals)
        Qassign = Dict{Any,Float64}()
        kk = 1
        for i in 1:n, j in i:n
            Qassign[L.qvec[kk]] = qvals[kk]
            kk += 1
        end
        return Float64(_to_number(subs(G, [k=>v for (k,v) in Qassign]...)))
    end

    # Classify every real critical point.
    println("\nAll real critical points:")
    all_pts = NamedTuple[]
    # Monodromy results don't support real_solutions() directly; extract from raw solutions
    raw_sols = HC.solutions(result)
    real_sols = Vector{Vector{Float64}}()
    for cs in raw_sols
        if all(abs(imag(z)) < 1e-6 for z in cs)
            push!(real_sols, Float64[real(z) for z in cs])
        end
    end
    println("(monodromy/polyhedral returned ", length(raw_sols), " total, ",
            length(real_sols), " real)")
    for sol in real_sols
        qvals = collect(sol[1:end-1])
        mu_v  = sol[end]
        Qnum = Q_from_qvec(qvals)
        is_pd = all(det(Qnum[1:k, 1:k]) > 1e-10 for k in 1:n)
        # Facet values, with tolerance for incidence detection.
        fac_vals = [Float64(sum(facets[i][k] * qvals[k] for k in 1:length(qvals)))
                    for i in 1:length(facets)]
        n_incident = count(x -> abs(x) < 1e-8, fac_vals)
        any_neg = any(x -> x < -1e-8, fac_vals)
        in_closure = !any_neg                # in closed cone (no facet violated)
        in_interior = in_closure && n_incident == 0
        detv = det(Qnum)
        G_Δ_val = eval_G_Δ(qvals)
        G_val   = G_Δ_val / detv^Float64(expG_denom_exp)
        push!(all_pts, (qvals=qvals, Qnum=Qnum, mu_v=mu_v,
                        is_pd=is_pd, in_closure=in_closure,
                        in_interior=in_interior, n_incident=n_incident,
                        G_Δ=G_Δ_val, G=G_val))
        loc = !is_pd ? "non-PD" :
              in_interior ? "INTERIOR" :
              in_closure ? "boundary (incident facets=$(n_incident))" :
              "outside-cone"
        @printf("  G=%+10.7f  Q=[%s]  %s\n",
                G_val, join([@sprintf("%+7.4f", v) for v in qvals], ", "), loc)
    end

    # Keep only PD critical points (in closure of the cone).
    pd_pts = [p for p in all_pts if p.is_pd && p.in_closure]
    println("\nPD critical points (in closed cone): ", length(pd_pts),
            "  (interior=", count(p -> p.in_interior, pd_pts), ")")

    # --- Aut(L) orbit grouping ------------------------------------------
    println("\n--- Aut(L) orbit grouping ---")
    # Lift each Aut(L) generator to act on Q: Q → g · Q · g^T (matching the
    # C++ matrix_in_t_space convention).
    j_raw = JSON3.read(read(ARGS[1], String))
    autL_gens = Vector{Matrix{Float64}}()
    for g in j_raw.autL_generators
        push!(autL_gens, _mat(g, Float64))
    end
    println("Loaded ", length(autL_gens), " Aut(L) generator(s) from JSON.")

    # Generate the full group (BFS from generators, dedup as matrices).
    function mat_key(M::Matrix{Float64})
        return tuple(round.(M; digits=8)...)
    end
    function close_group(gens::Vector{Matrix{Float64}}, max_order::Int=1000)
        Id = Matrix{Float64}(LinearAlgebra.I, size(gens[1])...)
        seen = Set([mat_key(Id)])
        elems = [Id]
        q = [Id]
        while !isempty(q)
            X = pop!(q)
            for g in gens
                M = g * X
                k = mat_key(M)
                if k ∉ seen
                    push!(seen, k)
                    push!(elems, M)
                    push!(q, M)
                    if length(elems) > max_order
                        return elems
                    end
                end
            end
        end
        return elems
    end
    group_elems = close_group(autL_gens, 1000)
    println("Closed Aut(L) has ", length(group_elems), " elements.")

    # Normalize a Q matrix to det=1 for comparison.
    function norm_Q(Q)
        d = det(Q)
        d > 0 || return nothing
        return Q ./ d^(1/n)
    end

    function are_equiv(Q1, Q2, gens; tol=1e-5)
        Q1n = norm_Q(Q1); Q2n = norm_Q(Q2)
        (Q1n === nothing || Q2n === nothing) && return false
        for g in gens
            if norm(g * Q1n * transpose(g) - Q2n) < tol
                return true
            end
        end
        return false
    end

    orbits = Vector{Vector{Int}}()   # indices into pd_pts
    assigned = falses(length(pd_pts))
    for i in 1:length(pd_pts)
        assigned[i] && continue
        orb = [i]; assigned[i] = true
        for j in (i+1):length(pd_pts)
            assigned[j] && continue
            if are_equiv(pd_pts[i].Qnum, pd_pts[j].Qnum, group_elems)
                push!(orb, j); assigned[j] = true
            end
        end
        push!(orbits, orb)
    end

    function is_near_rational(qvals; maxden=24, tol=1e-4)
        for v in qvals
            ok = false
            for d in 1:maxden
                if abs(v - round(v*d)/d) < tol
                    ok = true; break
                end
            end
            ok || return false
        end
        return true
    end

    println("\nDistinct Aut(L)-orbits among PD critical points: ", length(orbits))

    # --- Hessian classification on the slice det(Q) = 1 ----------------
    # At a critical point we have ∇G_Δ = μ ∇det. The constrained Hessian is
    # H_L = H_{G_Δ} − μ·H_{det}; its signature on the tangent space
    # {v : v · ∇det = 0} classifies min / max / saddle on the det=1 slice.
    d = length(L.qvec)
    detQ_sym = det_via_laplace(L.Q)
    grad_det_sym = [differentiate(detQ_sym, L.qvec[i]) for i in 1:d]
    H_G_sym = [differentiate(differentiate(G, L.qvec[i]), L.qvec[j])
               for i in 1:d, j in 1:d]
    H_det_sym = [differentiate(grad_det_sym[i], L.qvec[j]) for i in 1:d, j in 1:d]

    function eval_at(poly, qvals)
        ass = Dict{Any,Float64}()
        for (k, q) in enumerate(L.qvec)
            ass[q] = qvals[k]
        end
        return Float64(_to_number(subs(poly, [k=>v for (k,v) in ass]...)))
    end

    function classify_critical_point(p)
        H_G  = [eval_at(H_G_sym[i, j], p.qvals)  for i in 1:d, j in 1:d]
        H_det = [eval_at(H_det_sym[i, j], p.qvals) for i in 1:d, j in 1:d]
        grad_det = [eval_at(grad_det_sym[i], p.qvals) for i in 1:d]
        H_L = H_G .- Float64(p.mu_v) .* H_det
        # Tangent space basis: orthonormal complement of grad_det
        nrm = norm(grad_det)
        nrm == 0 && return ("degenerate", Float64[])
        nvec = grad_det ./ nrm
        # Build orthonormal basis of R^d, then drop the direction along nvec
        # Use QR of a matrix with first column = nvec
        M = hcat(reshape(nvec, d, 1), Matrix{Float64}(LinearAlgebra.I, d, d))
        F = qr(M)
        Qbasis = Matrix(F.Q)
        T = Qbasis[:, 2:end]  # d × (d-1), columns are tangent space basis
        H_L_tan = T' * H_L * T
        H_L_tan = (H_L_tan + H_L_tan') ./ 2   # symmetrize numerically
        evs = sort(LinearAlgebra.eigvals(H_L_tan))
        n_pos = count(e -> e >  1e-8, evs)
        n_neg = count(e -> e < -1e-8, evs)
        n_zero = length(evs) - n_pos - n_neg
        kind = if n_neg == 0 && n_zero == 0
            "local MIN"
        elseif n_pos == 0 && n_zero == 0
            "local MAX"
        elseif n_zero > 0
            "degenerate (zeros: $n_zero, pos: $n_pos, neg: $n_neg)"
        else
            "SADDLE (pos: $n_pos, neg: $n_neg)"
        end
        return (kind, evs)
    end

    println("\nClassification (Hessian signature on tangent space of det=1):")
    for (idx, orb) in enumerate(orbits)
        rep = pd_pts[orb[1]]
        rat = is_near_rational(rep.qvals) ? "rational" : "irrational"
        loc = rep.in_interior ? "INTERIOR" :
              "boundary (incident facets=$(rep.n_incident))"
        kind, evs = classify_critical_point(rep)
        @printf("  Orbit %2d: |orb|=%2d  G=%+10.7f  %s  %s\n",
                idx, length(orb), rep.G, loc, rat)
        @printf("           Q=[%s]\n",
                join([@sprintf("%+7.4f", v) for v in rep.qvals], ", "))
        @printf("           %s\n", kind)
        @printf("           eigenvalues of constrained Hessian: %s\n",
                join([@sprintf("%+8.5f", e) for e in evs], ", "))
    end
end

