#!/bin/bash
# The genera of pairs (L, sigma) by characteristic polynomial: for each line
# "LAT sig mass" of the table, GENUS_GlueIsometry factors on the genus of
# charpoly/LAT.txt with the factors of charpoly/LAT_sig.spec must give the
# total mass of the oracle (CharPolyOracle.g: the sum over the classes
# Lambda of the genus and the sigma in O(Lambda) of characteristic
# polynomial sig of 1 / |O(Lambda)|), every genus of pairs closed. For an
# odd lattice (I4, I8, I9, A2I2) the pairs of its even sublattice are
# enumerated, then their odd overlattices (IsometryGlueFactors.h).
#
# The data: MakeCharPolyTests.py (spec files, the jobs of the piece files
# for HermitianPieces.jl, the full table charpoly_tests.txt); the default
# table charpoly_tests_ci.txt is a part of it, the cases of a few seconds.
#
# Usage: RunCharPolyTests.sh [table] [timeout in seconds, default 600]
# (a case beyond the timeout is a failure)
set -e
BIN=${BIN:-../../src_genera/GENUS_GlueIsometry}
TABLE=${1:-charpoly_tests_ci.txt}
TMAX=${2:-600}
n_case=0
n_fail=0
while read -r lat sig expected; do
    case "$lat" in
        ''|'#'*) continue ;;
    esac
    n_case=$((n_case + 1))
    start=$(date +%s)
    status=0
    out=$(timeout "$TMAX" $BIN factors charpoly/${lat}.txt \
              charpoly/${lat}_${sig}.spec Summary stdout 2>/dev/null) || status=$?
    end=$(date +%s)
    total=$(echo "$out" | sed -n 's/^total mass = //p')
    complete=$(echo "$out" | sed -n 's/^complete = //p')
    echo "--- $lat $sig: total mass $total (expecting $expected), $((end - start))s ---"
    if [ "$status" != 0 ] || [ "$total" != "$expected" ] || [ "$complete" != "true" ]; then
        echo "ERROR: $lat $sig gave $total, complete=$complete, status=$status" >&2
        n_fail=$((n_fail + 1))
    fi
done < "$TABLE"
echo "$n_case cases, $n_fail failures"
if [ "$n_fail" != 0 ]; then
    exit 1
fi
echo "All characteristic polynomial tests passed"
touch CI_CONCLUSION
