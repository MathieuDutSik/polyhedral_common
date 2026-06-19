#!/usr/bin/env julia
#
# Batch monodromy runner. Reads an L-type JSON, runs monodromy with the
# given seed and time budget, writes a per-run summary file as JSON.
#
# Args: JSON-path, seed (int), budget-seconds, output-summary-path

push!(LOAD_PATH, @__DIR__)
include(joinpath(@__DIR__, "quantizer_polynomial.jl"))

# The script above runs at top-level when ARGS[1] is the input json. We
# subvert that here: re-load with adjusted ARGS so monodromy is invoked.
# (Simpler: just import functions explicitly and call them.)
# Since quantizer_polynomial.jl is structured as a script, we instead
# invoke it as a subprocess in run_batch.sh. This file is unused.
