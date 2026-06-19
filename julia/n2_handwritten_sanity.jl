#!/usr/bin/env julia
#
# n=2 sanity check for the quantizer critical-point pipeline.
#
# Uses Zimmermann's hand-derived polynomial (eq 6.3 of his 2017 PhD thesis):
#
#   G(Q) = (q11^2 q22 + q11 q22^2 - 2 q11 q12^2 - 2 q22 q12^2 - 2 q12^3)
#          / ( 24 * (q11 q22 - q12^2)^(3/2) )
#
# G is homogeneous of degree 0, so all critical rays can be normalized to
# det(Q) = 1. At det(Q) = 1, the constrained Lagrangian system
#
#   ∇num(Q) = mu * ∇det(Q),     det(Q) = 1
#
# is a polynomial system in (q11, q22, q12, mu). Solve it with
# HomotopyContinuation.jl, filter to real positive-definite solutions, check
# that the A_2 critical point (Q ~ [[2,-1],[-1,2]] normalized to det=1, i.e.
# (2/sqrt(3), 2/sqrt(3), -1/sqrt(3))) is among them with G = 5*sqrt(3)/108.

using HomotopyContinuation
using LinearAlgebra
using Printf

@var q11 q22 q12 mu

num = q11^2*q22 + q11*q22^2 - 2*q11*q12^2 - 2*q22*q12^2 - 2*q12^3
det_Q = q11*q22 - q12^2

vars = [q11, q22, q12]
eqs = Expression[differentiate(num, v) - mu * differentiate(det_Q, v) for v in vars]
push!(eqs, det_Q - 1)

F = System(eqs; variables = [q11, q22, q12, mu])

println("Solving Lagrangian system...")
result = solve(F; show_progress = false)

println("Path-tracking summary:")
println("  total solutions found  : ", length(solutions(result)))
println("  real solutions         : ", nreal(result))

println("\nReal critical points:")
expected_G = 5*sqrt(3)/108
let hits_A2 = 0
    for sol in real_solutions(result)
        q11_v, q22_v, q12_v, mu_v = sol
        detv = q11_v*q22_v - q12_v^2
        numv = q11_v^2*q22_v + q11_v*q22_v^2 - 2*q11_v*q12_v^2 - 2*q22_v*q12_v^2 - 2*q12_v^3
        G_v = numv / (24 * abs(detv)^(3/2))
        posdef = q11_v > 0 && detv > 0
        matches_A2 = posdef && abs(G_v - expected_G) < 1e-8
        if matches_A2
            hits_A2 += 1
        end
        @printf("  Q=(%+9.6f, %+9.6f, %+9.6f)  det=%+8.6f  mu=%+9.6f  G=%+9.6f  pd=%s  A2=%s\n",
                q11_v, q22_v, q12_v, detv, mu_v, G_v,
                posdef ? "yes" : "no ",
                matches_A2 ? "yes" : "no ")
    end
    global _hits_A2 = hits_A2
end
hits_A2 = _hits_A2

@printf("\nExpected A_2 quantizer value: G = 5*sqrt(3)/108 = %.10f\n", expected_G)
@printf("A_2 hits among real critical points: %d\n", hits_A2)

if hits_A2 == 0
    println("FAIL: did not find the A_2 critical point.")
    exit(1)
else
    println("OK: A_2 critical point identified.")
end
