#!/bin/bash
# The methods with their vector families (method:FAMILY), and the conjugacy
# of isometries in O(L) (LATT_TestConjugation).
# - the order of the automorphism group of A2, D4, E6, A2+A2 and of the
#   pair (A2+A2, swap of the summands) by every method that applies;
# - the reflections in two roots of A2 (g, h, rows: g A g^T = A) are
#   conjugate in O(A2), a reflection and its negative are not.
set -e
BIN=../../src_latt
TMP=$(mktemp -d)
trap 'rm -rf $TMP' EXIT

printf '1\n2 2\n2 -1\n-1 2\n' > $TMP/A2
printf '1\n4 4\n2 -1 0 0\n-1 2 -1 -1\n0 -1 2 0\n0 -1 0 2\n' > $TMP/D4
printf '1\n6 6\n2 -1 0 0 0 0\n-1 2 -1 0 0 0\n0 -1 2 -1 0 -1\n0 0 -1 2 -1 0\n0 0 0 -1 2 0\n0 0 -1 0 0 2\n' > $TMP/E6
printf '1\n4 4\n2 -1 0 0\n-1 2 0 0\n0 0 2 -1\n0 0 -1 2\n' > $TMP/A2A2
# the swap s of the two summands, as the second form A s^T
printf '2\n4 4\n2 -1 0 0\n-1 2 0 0\n0 0 2 -1\n0 0 -1 2\n4 4\n0 0 2 -1\n0 0 -1 2\n2 -1 0 0\n-1 2 0 0\n' > $TMP/A2A2swap

METHODS_ONE="graph graph:characteristic graph:iterated_shortest graph:root_iterated_shortest graph:shells_zbasis graph:zspan_closure iterated_shortest plesken_souvignier plesken_souvignier:shells_zbasis plesken_souvignier:characteristic heuristic vinberg weyl_chamber"
METHODS_SEVERAL="graph graph:iterated_shortest graph:shells_zbasis iterated_shortest plesken_souvignier plesken_souvignier:shells_zbasis plesken_souvignier:zspan_closure heuristic"

check_order() {
    lat=$1
    expected=$2
    methods=$3
    for m in $methods; do
        out=$($BIN/LATT_Automorphism gmp $TMP/$lat GAP_order stdout $m 2>/dev/null)
        echo "  $lat $m: $out"
        if [ "$out" != "return $expected;" ]; then
            echo "ERROR: |Aut($lat)| by $m is $out, expected $expected" >&2
            exit 1
        fi
    done
}
echo "--- automorphism group orders ---"
check_order A2 12 "$METHODS_ONE"
check_order D4 1152 "$METHODS_ONE"
check_order E6 103680 "$METHODS_ONE"
check_order A2A2 288 "$METHODS_ONE"
# the centralizer of the swap in O(A2+A2) = O(A2) wr S2: the diagonal O(A2)
# times the swap, order 24
check_order A2A2swap 24 "$METHODS_SEVERAL"

echo "--- conjugacy in O(A2) ---"
printf '2 2\n2 -1\n-1 2\n' > $TMP/A2m
printf '2 2\n-1 0\n1 1\n' > $TMP/g
printf '2 2\n1 1\n0 -1\n' > $TMP/h
printf '2 2\n1 0\n-1 -1\n' > $TMP/mg
for m in $METHODS_SEVERAL; do
    yes=$($BIN/LATT_TestConjugation gmp $TMP/A2m $TMP/g $TMP/h GAP stdout $m 2>/dev/null)
    no=$($BIN/LATT_TestConjugation gmp $TMP/A2m $TMP/g $TMP/mg GAP stdout $m 2>/dev/null)
    echo "  $m: $yes / $no"
    if [ "$yes" = "return false;" ] || [ -z "$yes" ] || [ "$no" != "return false;" ]; then
        echo "ERROR: LATT_TestConjugation with $m" >&2
        exit 1
    fi
done
rc=0
$BIN/LATT_TestConjugation gmp $TMP/A2m $TMP/g $TMP/h GAP stdout weyl_chamber \
    > /dev/null 2>&1 || rc=$?
if [ "$rc" != "1" ]; then
    echo "ERROR: weyl_chamber was not refused by LATT_TestConjugation" >&2
    exit 1
fi
echo "All method and family tests passed"
touch CI_CONCLUSION_FAMILIES
