#!/usr/bin/env julia
# Run one monodromy seed pass on a given JSON, save raw real solutions.
# Args: <json> <seed> <budget_seconds> <out_path>
import HomotopyContinuation as HC
using DynamicPolynomials
using MultivariatePolynomials
using LinearAlgebra
using JSON3
using Printf
using Random

include(joinpath(@__DIR__, "quantizer_polynomial.jl"))
