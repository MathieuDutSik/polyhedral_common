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

# Question 6 decided: the image of the local centraliser on the
# discriminant modules of the glued Hermitian components is computed, by
# lifting the solutions of g A g^T = A in the commutant (or by -1 and sigma
# when they already generate all the equivariant isometries), and found
# to be all of them. Here for D5_g5 and E6_g24, and for the isometry of
# order 19 of the determinant-30 genus (rank 19, glue Z/19).
image_case() {
    lat=$1
    isom=$2
    expected=$3
    echo "--- image of the local group for ($lat, $isom), expecting mass $expected ---"
    out=$($BIN mass isom_${lat}.txt isom_${isom}.txt none stdout 2>/dev/null)
    mass=$(echo "$out" | sed -n 's/^  mass:=\(.*\));$/\1/p')
    assumed=$(echo "$out" | sed -n 's/^  image_assumed:=\(.*\),$/\1/p')
    echo "mass = $mass, image_assumed = $assumed"
    if [ "$mass" != "$expected" ] || [ "$assumed" != "false" ]; then
        echo "ERROR: ($lat, $isom) gave mass $mass, image_assumed $assumed" >&2
        exit 1
    fi
}
image_case D5 D5_g5 1/64
image_case E6 E6_g24 1/144
image_case det30_L557 det30_g19 1/38

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
# Hermitian lattices over abelian CM fields that are not cyclotomic,
# given by the trace Gram matrix and the actions of a Z-basis of O_E, the
# field as the subfield of Q(zeta_n) fixed by H, with the masses of
# Hecke.jl: Q(sqrt -7), Q(sqrt -15) (odd ramified primes), Q(sqrt -2) and
# Q(sqrt -6) (the dyadic ramified prime, with relative different exponent
# 3), Q(i, sqrt 5) (dyadic, residue field F_4), and Q(sqrt 5, sqrt -7),
# where 2 is a common index divisor.
hab_case() {
    tag=$1
    n=$2
    H=$3
    expected=$4
    out=$($BIN hermitianablattice hab_${tag}_gram.txt hab_${tag}_omega.txt $n $H \
              2>/dev/null)
    val=$(echo "$out" | sed -n 's/return rec(mass:=\(.*\));/\1/p')
    echo "--- Hermitian lattice $tag over the subfield of Q(zeta_$n): mass=$val (expecting $expected) ---"
    if [ "$val" != "$expected" ]; then
        echo "ERROR: the Hermitian lattice $tag gave the mass $val" >&2
        exit 1
    fi
}
hab_case q7r3 7 1,2,4 4802/3
hab_case q15r3 15 1,2,4,8 426720
hab_case s5m7r2 35 1,4,9,11,16,29 25
hab_case s5m7r3 35 1,4,9,11,16,29 257985
hab_case q2r3 8 1,3 564375/32
hab_case q6r2 24 1,5,7,11 16
hab_case i5r2 20 1,9 300
# ramified primes of K that are not alone over q: 11, 3 and 19 split in
# K and ramify in E, and 2 splits in Q(sqrt 17) and ramifies in E
hab_case s5m11 55 4 5/12
hab_case s13m3 39 4 5/9
hab_case s17m1 68 9 3/4
hab_case s17m2 136 9,19 1950
hab_case s5m19 95 4 500
# The local masses (the products of the local factors) of quadratic
# lattices over totally real abelian fields, as given by Hecke.jl
# (local_mass): odd primes ramified or split, the dyadic primes unramified
# (Cho) and ramified (maximal, modular, chain) over 2.
quad_case() {
    tag=$1
    n=$2
    H=$3
    expected=$4
    out=$($BIN quadablocal quad_${tag}_gram.txt quad_${tag}_omega.txt $n $H \
              2>/dev/null)
    val=$(echo "$out" | sed -n 's/.*local_mass:=\([^,]*\),.*/\1/p')
    echo "--- quadratic lattice $tag over a real subfield of Q(zeta_$n): local mass=$val (expecting $expected) ---"
    if [ "$val" != "$expected" ]; then
        echo "ERROR: the quadratic lattice $tag gave the local mass $val" >&2
        exit 1
    fi
}
quad_case s5m2 5 1 9/2
quad_case z7pm4 7 1 1183
quad_case z15pm4 15 1 19113046875/2
quad_case z9pm3 9 1 567
quad_case s2m4 8 7 8760960
quad_case s3m4 12 11 6417840
quad_case s2s5m4 40 9,31 74880000000
quad_case z16pm4 16 1 5164032
quad_case z24pm4 24 1 503884800
quad_case z32pm3 32 1 9375024
quad_case s6m2 24 5,19 1/4
quad_case s7m3 28 3,9 155370096
quad_case s6m2b 24 5,19 1/8
# 2 splits completely in Q(sqrt 17, sqrt 41) and Q(sqrt 17, sqrt 33): a
# common index divisor, the idempotents refined from several elements
quad_case s17s41m2 697 2,21 1
quad_case s17s41m3 697 2,21 27720225/16
quad_case s17s33m2 561 2,25 1/16
quad_case s17s33m3 561 2,25 3645/16
# The masses of quadratic lattices over totally real abelian fields, as
# given by Hecke.jl (mass), with E = K(sqrt D) for the even ranks.
quadmass_case() {
    tag=$1
    n=$2
    H=$3
    expected=$4
    out=$($BIN quadablocal quadm_${tag}_gram.txt quadm_${tag}_omega.txt $n $H \
              2>/dev/null)
    val=$(echo "$out" | sed -n 's/return rec(mass:=\([^,]*\),.*/\1/p')
    echo "--- quadratic lattice $tag over a real subfield of Q(zeta_$n): mass=$val (expecting $expected) ---"
    if [ "$val" != "$expected" ]; then
        echo "ERROR: the quadratic lattice $tag gave the mass $val" >&2
        exit 1
    fi
}
quadmass_case s5m2 5 1 1
quadmass_case s2m2 8 7 45/4
quadmass_case s3m2 12 11 36
quadmass_case z7pm2 7 1 7/4
quadmass_case z9pm2 9 1 21/4
quadmass_case z15pm2 15 1 75/4
quadmass_case z16pm2 16 1 1/8
quadmass_case s5m3 5 1 15925/4
quadmass_case z7pm3 7 1 245245/12
quadmass_case z16pm3 16 1 5197365/16
# E = K(sqrt d) abelian but not K(sqrt D): Q(zeta_5) and Q(zeta_20)^+ over
# Q(sqrt 5), Q(zeta_16)^+ and the cyclic quartic CM field of Q(zeta_16)
# over Q(sqrt 2), found from the automorphisms of K and the split primes
quadmass_case s5abm2 5 1 1/4
quadmass_case s5abm4 5 1 205/32
quadmass_case s2abm2 8 7 1/4
quadmass_case s2abm4 8 7 205/192
# E not abelian: L(1 - r, chi_{E/K}) numerically (QuadraticHeckeL.h), the
# mass made exact by the Minkowski multiple, as Hecke.jl does
quadmass_case s13nabm2 13 4 99/4
quadmass_case z15pnabm2 15 1 1/2
# A common index divisor with primes of residue degree 2: 2 in
# Q(sqrt 5, sqrt 17) and Q(sqrt 13, sqrt 17), each prime of K over 2 with
# its own generator
quadmass_case s5s17m3 85 1,4,9,16,19,21,26,36,49,59,64,66,69,76,81,84 3375/2
quadmass_case s13s17m2 221 1,4,9,16,25,30,35,36,38,42,43,49,53,55,64,66,69,77,81,87,94,100,101,103,118,120,121,127,134,140,144,152,155,157,166,168,172,178,179,183,185,186,191,196,205,212,217,220 5
# A common index divisor ramified in K: 2 in Q(sqrt 3) K_31, K_31 the cubic
# field of conductor 31, three primes over 2 with e = 2 and f = 1
quadmass_case s3c31m2 372 1,23,35,47,61,85,95,97,109,157,215,263,275,277,287,311,325,337,349,371 3
# The same with the other methods for the L-values: lvalue=double (every
# term in double precision) and lvalue=arb (the integer certified by Arb),
# the latter skipped when the programs are built without flint
# (make ENABLE_FLINT_SUPPORT=1)
quadmass_method_case() {
    method=$1
    tag=$2
    n=$3
    H=$4
    expected=$5
    # without flint, lvalue=arb stops with an error, hence the || true
    out=$($BIN quadablocal quadm_${tag}_gram.txt quadm_${tag}_omega.txt $n $H \
              lvalue=$method 2>&1 || true)
    if echo "$out" | grep -q "built without flint"; then
        echo "--- quadratic lattice $tag with $method: skipped, no flint ---"
        return
    fi
    val=$(echo "$out" | sed -n 's/return rec(mass:=\([^,]*\),.*/\1/p')
    echo "--- quadratic lattice $tag with $method: mass=$val (expecting $expected) ---"
    if [ "$val" != "$expected" ]; then
        echo "ERROR: the quadratic lattice $tag gave with $method the mass $val" >&2
        exit 1
    fi
}
quadmass_method_case double s13nabm2 13 4 99/4
quadmass_method_case double z15pnabm2 15 1 1/2
quadmass_method_case arb s13nabm2 13 4 99/4
quadmass_method_case arb z15pnabm2 15 1 1/2
# The Hermitian and quadratic masses over abelian fields are generic in the
# rational type: the same masses with mpq_class, boost mpq_rational and,
# when built with flint, fmpq_class.
types_case() {
    mode=$1
    tag=$2
    n=$3
    H=$4
    out=$(../../src_genera/Test_MassTypes $mode ${tag}_gram.txt ${tag}_omega.txt $n $H 2>/dev/null)
    echo "--- the mass of $tag over the rational types: $(echo $out | tr '\n' ' ') ---"
    if ! echo "$out" | grep -q "^agree$"; then
        echo "ERROR: the rational types disagree on $tag" >&2
        exit 1
    fi
}
types_case herm hab_q7r3 7 1,2,4
types_case herm hab_s17m2 136 9,19
types_case quad quadm_z7pm2 7 1
types_case quad quad_s2m4 8 7
types_case quad quad_s17s41m3 697 2,21
herm_case 47 1 mass 695/94
herm_case 19 1 mass 1/38
herm_case 16 1 mass 1/16
herm_case 12 1 mass 1/12
herm_case 39 1 mass 1/39
herm_case 3 2 standard_mass 1/72
herm_case 16 2 standard_mass 5/768

# Finite abelian groups (GENUS_IsometryMass group), the components being
# the Galois orbits of characters: on E8+E8, <-1, r> has the mass of the
# root reflection r alone; two orthogonal reflections (in one summand or in
# both, the same genus of pairs), an element of order 3 times a commuting
# reflection, and two orthogonal reflections of the genus of determinant
# 25. The masses are those closed by GENUS_EnumerateIsometry from one seed.
group_mass_case() {
    lat=$1
    gens=$2
    expected=$3
    echo "--- mass of ($lat, <$gens>), expecting $expected ---"
    out=$($BIN group isom_${lat}.txt isom_${lat}_${gens}.txt stdout 2>/dev/null)
    mass=$(echo "$out" | sed -n 's/^  mass:=\(.*\));$/\1/p')
    echo "mass = $mass"
    if [ "$mass" != "$expected" ]; then
        echo "ERROR: ($lat, <$gens>) gave the mass $mass, expected $expected" >&2
        exit 1
    fi
}
group_mass_case E8E8_16 mIrefl 691/1156946589646848000
group_mass_case E8E8_16 klein1 42151/385648863215616000
group_mass_case E8E8_16 klein2 42151/385648863215616000
group_mass_case E8E8_16 c6 691/397303087104000
group_mass_case d25 klein 155/9216
# Non-abelian groups: W(A2) = S3 and W(A3) = S4 generated by root
# reflections of E8+E8. The reflection representations have End_G = Q and
# multiplicity one, so their masses are 1/2 and exact; both masses are
# closed by GENUS_EnumerateIsometry from one seed (2 and 3 classes).
group_mass_case E8E8_16 S3 691/20659760529408000
group_mass_case E8E8_16 S4 691/765176315904000
# In higher multiplicity the local factors [O(X_p) : Stab(L_p)] at the p
# dividing |G| enter: W(A2) acting diagonally on the A2^3 of E6 (the
# multiplicity lattice gives 1/24, the factor is 2) and on two orthogonal
# A2 of E8+E8 (glue factor 8); both masses are closed by the enumeration.
group_mass_case E6diag A2diag 1/12
# The Frobenius group F21 = C7 : C3 on seven coordinates: its faithful
# rational representation has End_G = Q(sqrt -7), a CM field that is not
# cyclotomic, so the multiplicity lattice is Hermitian over Q(sqrt -7) (of
# rank one on A6 and A8, of rank two on A6 + A6, local factor 1 at 7).
group_mass_case A6 F21 1/2
group_mass_case A8 F21 1/2
group_mass_case A6A6 F21 1/8
group_mass_case E8E8_16 S3diag2 143/1203948748800
# The dihedral group D5 = <(12345), (25)(34)> in W(A4): End_G has the
# totally real centre Q(sqrt 5), so the multiplicity lattice is quadratic
# over Q(sqrt 5). Diagonally on A4 + A4 (rank two, mass 1/8), and on two
# A4, one in each E8 of E8+E8 (glued at 5 to the fixed lattice, mass
# 13/921600, closed by GENUS_EnumerateIsometry with 2 classes).
group_mass_case A4A4 D5 1/8
group_mass_case E8E8D5 D5 13/921600
# Rank three over O_K, where the local factor at the ramified prime 5 is
# found by the pi-adic search over O_K (the solutions modulo 5 of the
# trace form are about 10^4 times too many): D5 diagonally on A4^3 (factor
# 1, the image of order 30000), and on the three A4 of E8 + A4 (factor 2,
# so the mass 1/8 is that of the seed, with a centraliser of order 8)
group_mass_case A4A4A4 D5 1/48
group_mass_case E8A4 D5 1/8

# The brute-force check for the groups <g, tau>, tau in C_{O(L)}(g); the
# ratios are integers but not 1, as (g, tau) and (g, g tau) generate the
# same group.
group_check_case() {
    lat=$1
    isom=$2
    echo "--- group check of ($lat, $isom) ---"
    out=$($BIN groupcheck isom_${lat}.txt isom_${isom}.txt stdout 2>/dev/null)
    echo "$out"
    if ! echo "$out" | grep -q "^check: every ratio is an integer$"; then
        echo "ERROR: the group check of ($lat, $isom) failed" >&2
        exit 1
    fi
}
group_check_case D4 D4_inv1
group_check_case D4 D4_w01
group_check_case D5 D5_inv2

echo "All isometry mass tests passed"
touch CI_CONCLUSION
