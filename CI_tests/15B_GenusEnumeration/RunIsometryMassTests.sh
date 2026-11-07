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
# E8+E8 with the exchange of the factors: the glue is the graph of an
# anti-isometry between two copies of the discriminant form of E8(2), the
# plus type quadratic space of dimension 8 over F_2, so the glue factor is
# |O+_8(2)| = 348364800, beyond the listing and computed by the graph
# automorphism groups; the mass is 1/(2 |W(E8)|), that of the single pair.
mass_case E8E8_16 swap16 1/1393459200

check_case A2A2 A2A2_inv2
check_case D4 D4_inv1
check_case D4 D4_inv2
check_case D5 D5_inv2
check_case A2 A2_w01
check_case D4 D4_w01
check_case D4 D4_w012
check_case D4 D4_w0123
check_case E6 E6_w01

# The same check with all the masses computed (no mass file), on glued
# Hermitian components of rank two, where the image of the local unitary
# group is taken to be all the equivariant isometries of the discriminant
# form (Question 6 of the paper): an element of order 4 of D5 (Z[i] in
# rank 2, glued to the fixed line) and one of order 6 of E6 (Z[zeta_3] in
# rank 2 glued to Z[zeta_6] in rank 1).
check_case_computed() {
    lat=$1
    isom=$2
    echo "--- check of ($lat, $isom), masses computed ---"
    out=$($BIN check isom_${lat}.txt isom_${isom}.txt none stdout 2>/dev/null)
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
check_case_computed D5 D5_g5
check_case_computed E6 E6_g24

# The masses of the Hermitian lattices of rank one over Q(zeta_d), h^- / w,
# and standard masses, checked against Hecke.jl: a prime (relative class
# number 695 for 47), a power of 2 (the dyadic prime), composite d (no
# ramified prime) and rank two.
herm_case() {
    d=$1
    m=$2
    key=$3
    expected=$4
    out=$($BIN hermitian $d $m 2>/dev/null)
    val=$(echo "$out" | sed -n "s/.*${key}:=\([0-9/]*\).*/\1/p")
    echo "--- hermitian d=$d m=$m $key=$val (expecting $expected) ---"
    if [ "$val" != "$expected" ]; then
        echo "ERROR: hermitian d=$d m=$m gave $key=$val" >&2
        exit 1
    fi
}
# Hermitian lattices of higher rank, as trace lattices with the action of
# zeta_d (from Hecke.jl's trace_lattice_with_isometry), with the masses of
# Hecke.jl: Q(zeta_12) in rank 3, Q(zeta_15) in rank 2 with a prime near
# 9.3 10^10 in the determinant, and the ramified prime of Q(zeta_7) and
# Q(zeta_9) in rank 2.
herm_lattice_case() {
    name=$1
    d=$2
    expected=$3
    out=$($BIN mass herm_${name}_gram.txt herm_${name}_isom.txt none stdout \
              2>/dev/null)
    val=$(echo "$out" | grep "rec(d:=$d," | sed -n 's/.*mass:=\([0-9/]*\).*/\1/p')
    echo "--- Hermitian lattice $name over Q(zeta_$d): mass=$val (expecting $expected) ---"
    if [ "$val" != "$expected" ]; then
        echo "ERROR: the Hermitian lattice $name gave the mass $val" >&2
        exit 1
    fi
}
herm_lattice_case cases7 12 52486783453/384
herm_lattice_case cases12 15 1244881/900
herm_lattice_case ram11 7 691590/49
herm_lattice_case ram15 9 1200500/81
# the dyadic ramified prime: of Q(zeta_8) for a lattice maximal there, of
# Q(i) in rank 3 and of Q(zeta_16) for lattices that are not (the chain
# method with Kirschmer's local isometry test)
herm_lattice_case dy21 8 1741/96
herm_lattice_case dy10 4 1333/16
herm_lattice_case dy37 16 28824005/32
herm_case 47 1 mass 695/94
herm_case 19 1 mass 1/38
herm_case 16 1 mass 1/16
herm_case 12 1 mass 1/12
herm_case 39 1 mass 1/39
herm_case 3 2 standard_mass 1/72
herm_case 16 2 standard_mass 5/768

echo "All isometry mass tests passed"
touch CI_CONCLUSION
