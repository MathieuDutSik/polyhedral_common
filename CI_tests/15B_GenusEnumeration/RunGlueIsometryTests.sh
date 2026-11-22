#!/bin/bash
# The seeds of the enumeration of pairs by gluing (GENUS_GlueIsometry): the
# sigma-stable primitive extensions of (L_1, g_1) (+) (L_2, g_2), each
# genus of pairs met enumerated by its mass.
#
# - A1 with -1 and E7 with 1, determinant 1: E8 with a root reflection,
#   one class, centralizer 2 |W(E7)|;
# - D4 with 1 and D4 with -1, any determinant: the glue orders 1, 2 and 4,
#   three genera of one class each, the last one E8 with an involution of
#   type D4 (centralizer 6 |W(D4)|^2);
# - A1 with -1 and E7+E8 with 1, determinant 1: the genus of II_16 with a
#   root reflection, whose second class (on D16+) the enumeration finds;
# - A2 with an element of order 3 and E6+E8 with 1, determinant 1: the
#   genus of II_16 with an element of order 3 of an A2, two classes;
# - D4+D4 with 1 and with -1, determinant 1: glue (Z/2)^4, one class.
# In every case all the genera must close their mass: complete = true.
set -e
BIN=../../src_genera/GENUS_GlueIsometry

glue_case() {
    l1=$1
    i1=$2
    l2=$3
    i2=$4
    det=$5
    expected_genera=$6
    expected_classes=$7
    echo "--- ($l1, $i1) + ($l2, $i2), det $det, expecting $expected_genera genera, $expected_classes classes ---"
    out=$($BIN enumerate glue_${l1}.txt glue_${i1}.txt glue_${l2}.txt \
              glue_${i2}.txt $det Summary stdout 2>/dev/null)
    echo "$out"
    nb_genera=$(echo "$out" | sed -n 's/^genera of pairs = //p')
    nb_classes=$(echo "$out" | sed -n 's/^classes of pairs = //p')
    complete=$(echo "$out" | sed -n 's/^complete = //p')
    if [ "$nb_genera" != "$expected_genera" ] ||
           [ "$nb_classes" != "$expected_classes" ] ||
           [ "$complete" != "true" ]; then
        echo "ERROR: ($l1, $i1) + ($l2, $i2) gave $nb_genera genera, $nb_classes classes, complete=$complete" >&2
        exit 1
    fi
}

glue_case A1 A1_mI E7 E7_id 1 1 1
glue_case D4 D4_id D4 D4_mI any 3 3
glue_case A1 A1_mI E7E8 E7E8_id 1 1 2
glue_case A2 A2_w E6E8 E6E8_id 1 1 2
glue_case D4D4 D4D4_id D4D4 D4D4_mI 1 1 1

# Every genus of pairs (L, sigma), sigma an involution with a +1-eigenspace
# of rank r and L in a given genus, from the genera of the pieces (even, of
# ranks r and n - r, d_+ d_- = d 4^a, the odd local symbols adding up to
# those of L): the total mass must be the sum over the classes Lambda of
# the genus of #{sigma in O(Lambda) of that type} / |O(Lambda)|, computed
# in GAP (W(E8) for E8; the two classes A2^4 and its neighbour of
# determinant 81).
invol_case() {
    lat=$1
    r=$2
    expected=$3
    out=$($BIN involutions invol_${lat}.txt $r Summary stdout 2>/dev/null)
    total=$(echo "$out" | sed -n 's/^total mass = //p')
    complete=$(echo "$out" | sed -n 's/^complete = //p')
    echo "--- involutions of rank $r on the genus of $lat: total mass $total (expecting $expected) ---"
    if [ "$total" != "$expected" ] || [ "$complete" != "true" ]; then
        echo "ERROR: involutions $lat r=$r gave $total, complete=$complete" >&2
        exit 1
    fi
}
invol_case E8 0 1/696729600
invol_case E8 1 1/5806080
invol_case E8 2 1/184320
invol_case E8 3 1/18432
invol_case E8 4 37/221184
invol_case E8 5 1/18432
invol_case E8 8 1/696729600
invol_case A2x4 1 13/82944
invol_case A2x4 2 805/497664
invol_case A2x4 3 229/27648
invol_case A2x4 4 6925/331776
invol_case A2x4 7 13/82944

# The order-d isometries sigma = 1 (+) sigma_2 with char poly
# (x - 1)^r Phi_d^b on a genus: the +1-pieces from the genus symbols, the
# pieces (L_2, sigma_2) one per hermitian genus over Z[zeta_d], made by
# HermitianPieces.jl (Hecke.jl): d = 5, 8, 10, 12 over a real centre
# K_0 != Q, the trace form integral for the scales in the inverse of the
# absolute different. The expected masses are the sums of
# 1 / |O(Lambda)| over the classes Lambda of the genus and the elements of
# O(Lambda) of that char poly (GAP).
herm_case() {
    lat=$1
    d=$2
    b=$3
    r=$4
    expected=$5
    out=$($BIN pieces invol_${lat}.txt $r herm_${lat}_d${d}_b${b}.txt Summary stdout 2>/dev/null)
    total=$(echo "$out" | sed -n 's/^total mass = //p')
    complete=$(echo "$out" | sed -n 's/^complete = //p')
    echo "--- Phi_$d^$b on the genus of $lat: total mass $total (expecting $expected) ---"
    if [ "$total" != "$expected" ] || [ "$complete" != "true" ]; then
        echo "ERROR: pieces $lat d=$d b=$b gave $total, complete=$complete" >&2
        exit 1
    fi
}
herm_case E8 3 1 6 1/311040
herm_case E8 3 2 4 1/2592
herm_case E8 3 3 2 1/7776
herm_case E8 3 4 0 1/155520
herm_case E8 4 1 6 0
herm_case E8 4 2 4 1/18432
herm_case E8 4 3 2 0
herm_case E8 4 4 0 1/46080
herm_case E8 6 4 0 1/155520
herm_case E8 5 1 4 1/1200
herm_case E8 5 2 0 1/600
herm_case E8 8 2 0 1/192
herm_case E8 10 2 0 1/600
herm_case E8 12 2 0 1/288
herm_case A2x4 3 1 6 19/62208
herm_case A2x4 3 2 4 193/20736
herm_case A2x4 3 3 2 109/15552
herm_case A2x4 3 4 0 19/31104

# The genera of even lattices of rank n and determinant d against
# SageMath, for n <= 8 and d <= 64
TSE=../../src_genera/Test_GenusSymbolEnumeration
ours=$(for n in 1 2 3 4 5 6 7 8; do $TSE count $n 1 64; done)
if [ "$ours" != "$(grep -v '^#' genus_symbol_counts_sage.txt)" ]; then
    echo "ERROR: the numbers of genera differ from those of SageMath" >&2
    exit 1
fi
echo "--- the numbers of genera of 512 (rank, determinant) agree with SageMath ---"
ours=$(for n in 1 2 3 4 5 6 7 8; do $TSE countodd $n 1 64; done)
if [ "$ours" != "$(grep -v '^#' genus_symbol_counts_odd_sage.txt)" ]; then
    echo "ERROR: the numbers of odd genera differ from those of SageMath" >&2
    exit 1
fi
echo "--- the numbers of odd genera of 512 (rank, determinant) agree with SageMath ---"

# The genera of hermitian lattices over Q(zeta_d) / Q(zeta_d)^+ against
# hermitian_genera of Hecke.jl (HermitianGenusCounts.jl), d not a power of 2
THG=../../src_genera/Test_HermitianGenus
ours=$(grep -v '^#' hermitian_genus_counts_hecke.txt | while read -r d b N c; do $THG count $d $b $N; done)
theirs=$(grep -v '^#' hermitian_genus_counts_hecke.txt | awk '{print $1, $2, $3, $4}')
if [ "$ours" != "$theirs" ]; then
    echo "ERROR: the numbers of hermitian genera differ from those of Hecke" >&2
    exit 1
fi
echo "--- the numbers of hermitian genera of $(echo "$theirs" | wc -l | tr -d ' ') cases agree with Hecke ---"
# The primes of Q(zeta_d)^+ and their type in Q(zeta_d), from the Dirichlet
# characters (AbelianField.h), against the explicit formulas
if ! r=$($THG decomposition 2>&1 | tail -1) || [ "${r##*errors=}" != "0" ]; then
    echo "ERROR: the prime decompositions differ: $r" >&2
    exit 1
fi
echo "--- the prime decompositions agree with the formulas ($r) ---"

# The pieces of GENUS_HermitianPieces (HermitianLattice.h, without Hecke)
# against those of Hecke.jl: every genus of each file must be found
THL=../../src_genera/Test_HermitianLattice
n_cmp=0
while read -r d b N f; do
    case "$d" in ''|'#'*) continue ;; esac
    r=$($THL compare "$f" "$d" "$b" "$N" 2>/dev/null | tail -1)
    if [ "${r##*missing=}" != "0" ]; then
        echo "ERROR: the hermitian pieces of $f: $r" >&2
        exit 1
    fi
    n_cmp=$((n_cmp + 1))
done < hermitian_pieces_compare_ci.txt
echo "--- the hermitian pieces of $n_cmp files of Hecke are all found ---"

echo "All gluing tests passed"
touch CI_CONCLUSION
