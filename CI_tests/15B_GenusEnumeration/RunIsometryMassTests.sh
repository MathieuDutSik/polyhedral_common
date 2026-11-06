#!/bin/bash
# The mass of the genus of a lattice with an isometry (GENUS_IsometryMass).
#
# First, masses known by hand: A1+A1 with the swap (1/4), E8 with -1
# (1/|W(E8)|) and with a root reflection (1/(2 |W(E7)|)), and E8 with
# products of two and three orthogonal reflections, whose glue factors are
# 2 and 6 (the centralisers W(A1^2) x W(D6) extended by the swap, and so
# on).
#
# Then the check mode, independent of the glue computation: on a lattice
# whose genus has one class, the elements of O(L) of the characteristic
# polynomial of g are grouped by the mass of their pair, and for each group
# the count divided by |O(L)| times the mass must be an integer, the number
# of genera of pairs in it. The involutions of A2+A2, D4 and D5 give glue
# factors 2, 3 and 6; the elements of order 3, 4 and 6 of A2, D4 and E6
# have Hermitian components of rank one over Z[zeta_3] or Z[i], of mass
# 1/6 or 1/4 (isom_herm_rank1.txt), glued to orthogonal ones.
set -e
BIN=../../src_genera/GENUS_IsometryMass

mass_case() {
    lat=$1
    isom=$2
    expected=$3
    echo "--- mass of ($lat, $isom), expecting $expected ---"
    out=$($BIN mass isom_${lat}.txt isom_${isom}.txt none stdout 2>/dev/null)
    mass=$(echo "$out" | sed -n 's/^  mass:=\(.*\));$/\1/p')
    echo "mass = $mass"
    if [ "$mass" != "$expected" ]; then
        echo "ERROR: ($lat, $isom) gave the mass $mass, expected $expected" >&2
        exit 1
    fi
}

check_case() {
    lat=$1
    isom=$2
    echo "--- check of ($lat, $isom) ---"
    out=$($BIN check isom_${lat}.txt isom_${isom}.txt isom_herm_rank1.txt \
              stdout 2>/dev/null)
    echo "$out"
    if ! echo "$out" | grep -q "^check: every ratio is an integer$"; then
        echo "ERROR: the check of ($lat, $isom) failed" >&2
        exit 1
    fi
    if echo "$out" | grep "^count=" | grep -v -q " ratio=1 "; then
        echo "ERROR: the check of ($lat, $isom) has a ratio other than 1" >&2
        exit 1
    fi
}

mass_case A1A1 swap 1/4
mass_case E8 minusI 1/696729600
mass_case E8 refl 1/5806080
mass_case E8 E8_inv2 1/184320
mass_case E8 E8_inv3 1/18432

check_case A2A2 A2A2_inv2
check_case D4 D4_inv1
check_case D4 D4_inv2
check_case D5 D5_inv2
check_case A2 A2_w01
check_case D4 D4_w01
check_case D4 D4_w012
check_case D4 D4_w0123
check_case E6 E6_w01

echo "All isometry mass tests passed"
touch CI_CONCLUSION
