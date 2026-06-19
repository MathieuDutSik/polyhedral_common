#!/usr/bin/env julia
#
# Run ONE monodromy pass on a given L-type JSON and append/write critical
# points to an output text file. Usage:
#
#   julia --project monodromy_one_run.jl <json> <seed> <timeout_seconds> <out_path>

using DynamicPolynomials
using MultivariatePolynomials
import HomotopyContinuation as HC
using JSON3
using LinearAlgebra
using Printf
using Random

include(joinpath(@__DIR__, "quantizer_lib.jl"))

if length(ARGS) < 4
    println(stderr, "Usage: monodromy_one_run.jl <json> <seed> <timeout_s> <out_path>")
    exit(1)
end

json_path = ARGS[1]
seed = parse(Int, ARGS[2])
timeout_s = parse(Float64, ARGS[3])
out_path = ARGS[4]

Random.seed!(seed)
println("monodromy_one_run: json=$json_path seed=$seed timeout=$(timeout_s)s out=$out_path")

L = load_ltype(json_path)
println("Loaded L-type: dim=$(L.n), |star|=$(length(L.delaunay)), |T(Q)|=$(length(L.triangulation))")
println("Building G_Δ...")
G, _ε = build_G_Delta(L)
println("G_Δ built, terms=", length(terms(polynomial(G))))

result, hc_qvars = find_critical_points_monodromy(L, G; seed=UInt32(seed),
                                                  timeout=timeout_s,
                                                  max_sols=10000)

raw_sols = HC.solutions(result)
real_sols = Vector{Vector{Float64}}()
for cs in raw_sols
    if all(abs(imag(z)) < 1e-6 for z in cs)
        push!(real_sols, Float64[real(z) for z in cs])
    end
end
println("Found ", length(raw_sols), " complex, ", length(real_sols), " real")

# Write minimal CSV: one line per real solution = qvec values, μ value
open(out_path, "w") do io
    println(io, "# json=", json_path, " seed=", seed, " budget=", timeout_s, "s")
    println(io, "# complex=", length(raw_sols), " real=", length(real_sols))
    println(io, "# dim=", L.n)
    for sol in real_sols
        println(io, join(string.(sol), ","))
    end
end
println("Wrote ", length(real_sols), " real solutions to ", out_path)
