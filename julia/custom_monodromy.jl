#!/usr/bin/env julia
#
# Custom Newton-based monodromy solver using our numerical Lagrangian
# evaluator (in quantizer_numerical.jl). Bypasses HC.jl's symbolic-only
# AbstractSystem conversion.
#
# Strategy: Given (x_start, p_start) such that F(x_start, p_start) = 0,
# we walk loops p_start → p_mid → p_start in C^{d+1} (random complex
# midpoints). Each loop tracks every known solution; new endpoints (mod
# rounding) are added to the solution set. Repeat until no new solutions
# appear after several loops.

using LinearAlgebra
using Random

include(joinpath(@__DIR__, "quantizer_numerical.jl"))

# ---------------------------------------------------------------------------
# Newton corrector: given (x, p) refine x so that F(x, p) ≈ 0.

function newton_correct(sys::QuantizerLagrangianSystem,
                        x::AbstractVector{ComplexF64},
                        p::AbstractVector{ComplexF64};
                        tol::Float64 = 1e-10,
                        max_iter::Int = 30)
    n_eq = sys.L.d + 1
    u = zeros(ComplexF64, n_eq)
    J = zeros(ComplexF64, n_eq, n_eq)
    x_cur = copy(x)
    for iter in 1:max_iter
        HC.evaluate_and_jacobian!(u, J, sys, x_cur, p)
        nrm = norm(u)
        if nrm < tol
            return x_cur, true
        end
        local dx
        try
            dx = J \ u
        catch
            return x_cur, false
        end
        x_cur -= dx
        if norm(dx) < 1e-13
            HC.evaluate!(u, sys, x_cur, p)
            return x_cur, norm(u) < 1e3 * tol
        end
    end
    HC.evaluate!(u, sys, x_cur, p)
    return x_cur, norm(u) < 1e3 * tol
end

# Predict-correct path tracker from (x_start, p_start) to p_end via a
# straight line in parameter space.
function track_path(sys::QuantizerLagrangianSystem,
                    x_start::AbstractVector{ComplexF64},
                    p_start::AbstractVector{ComplexF64},
                    p_end::AbstractVector{ComplexF64};
                    n_steps::Int = 60,
                    verbose::Bool = false)
    n_eq = sys.L.d + 1
    x = ComplexF64.(x_start)
    u = zeros(ComplexF64, n_eq)
    J = zeros(ComplexF64, n_eq, n_eq)
    for i in 1:n_steps
        t = i / n_steps
        p = (1 - t) * p_start + t * p_end
        # Predictor: linear extrapolation (= last x; basic predictor)
        # Corrector: Newton
        x_new, ok = newton_correct(sys, x, p; tol=1e-10, max_iter=30)
        if !ok
            verbose && println("  track_path: Newton failed at t=$t")
            return nothing
        end
        x = x_new
    end
    return x
end

# Tolerance-based dedup: a new solution is recognized as different from
# all existing ones if no existing solution is within `atol` (Euclidean
# distance in C^{d+1}).
function is_new_solution(x_new::AbstractVector{ComplexF64},
                         existing::Vector{Vector{ComplexF64}};
                         atol::Float64 = 1e-3)
    for x in existing
        if norm(x_new .- x) < atol
            return false
        end
    end
    return true
end

# Custom monodromy: from a single start (x_start, p_start) (where F=0),
# walk random loops in p-space and grow the solution set.
function custom_monodromy(sys::QuantizerLagrangianSystem,
                         x_start::AbstractVector{ComplexF64},
                         p_start::AbstractVector{ComplexF64};
                         timeout::Float64 = 300.0,
                         max_loops_no_progress::Int = 8,
                         loop_radius::Float64 = 2.0,
                         n_steps_per_segment::Int = 100,
                         dedup_atol::Float64 = 1e-3,
                         seed::Int = 0,
                         verbose::Bool = true)
    seed != 0 && Random.seed!(seed)
    p_dim = length(p_start)
    sols = Vector{Vector{ComplexF64}}()
    push!(sols, ComplexF64.(x_start))
    t_start = time()
    loops_done = 0
    loops_no_progress = 0
    while loops_no_progress < max_loops_no_progress
        time() - t_start > timeout && (verbose && println("  monodromy timeout"); break)
        loops_done += 1
        # Generate random complex midpoint (large enough to cross
        # discriminants in p-space).
        Δ = loop_radius * randn(ComplexF64, p_dim)
        p_mid = p_start + Δ
        n_added = 0
        # Snapshot current sols before this loop iteration to avoid an
        # infinite cascade if new sols are produced mid-loop.
        sols_snapshot = copy(sols)
        for x in sols_snapshot
            time() - t_start > timeout && break
            x_mid = track_path(sys, x, p_start, p_mid;
                               n_steps=n_steps_per_segment)
            x_mid === nothing && continue
            x_back = track_path(sys, x_mid, p_mid, p_start;
                                n_steps=n_steps_per_segment)
            x_back === nothing && continue
            # Verify that x_back is actually a solution (F=0 at p_start).
            u_check = zeros(ComplexF64, sys.L.d + 1)
            HC.evaluate!(u_check, sys, x_back, p_start)
            if norm(u_check) > 1e-4
                continue   # path drifted; reject
            end
            if is_new_solution(x_back, sols; atol=dedup_atol)
                push!(sols, x_back)
                n_added += 1
            end
        end
        if n_added > 0
            loops_no_progress = 0
            verbose && println("  loop $loops_done: +$n_added new (total=$(length(sols)))")
        else
            loops_no_progress += 1
            verbose && println("  loop $loops_done: no new (no_progress=$loops_no_progress, total=$(length(sols)))")
        end
    end
    elapsed = time() - t_start
    verbose && println("custom_monodromy done: $(length(sols)) solutions in $loops_done loops, $(round(elapsed; digits=1))s")
    return sols
end
