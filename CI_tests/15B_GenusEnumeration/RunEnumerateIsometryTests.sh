#!/bin/bash
# Enumeration of the classes of pairs (L, g) in the genus of a pair
# (GENUS_EnumerateIsometry), stopped by the mass of GENUS_IsometryMass.
#
# - E8+E8 with -1: the pairs are the lattices, so this is the genus of the
#   even unimodular lattices of rank 16 (E8+E8 and D16+);
# - the genus of rank 8 and determinant 25 with -1: its 3 classes;
# - E8+E8 with a root reflection: 2 classes, the centralizers being
#   2 |W(E7)| |W(E8)| and 2^15 16! / 240;
# - E8+E8 with an element of order 3 of an A2 subsystem: 2 classes, a
#   Hermitian component over Z[zeta_3] of mass 1/6, found at p = 2;
# - the genus of determinant 25 with a root reflection: 4 classes of pairs,
#   glued at 2.
# In every case the mass must close: complete = true.
set -e
BIN=../../src_genera/GENUS_EnumerateIsometry

run_case() {
    lat=$1
    isom=$2
    expected=$3
    echo "--- ($lat, $isom), expecting $expected classes of pairs ---"
    out=$($BIN eiso_${lat}.txt eiso_${isom}.txt eiso_mass_${isom}.txt \
              Summary stdout 2>/dev/null)
    echo "$out"
    nb=$(echo "$out" | sed -n 's/^classes of pairs = //p')
    complete=$(echo "$out" | sed -n 's/^complete = //p')
    if [ "$nb" != "$expected" ] || [ "$complete" != "true" ]; then
        echo "ERROR: ($lat, $isom) gave $nb classes, complete=$complete" >&2
        exit 1
    fi
}

run_case E8E8 mI16 2
run_case d25 mI8 3
run_case E8E8 refl16 2
run_case E8E8 ord3 2
run_case d25 refl8 4

echo "All isometry enumeration tests passed"
touch CI_CONCLUSION
