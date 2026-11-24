#!/bin/bash
# The programs built with the SANITY_CHECK macros of src_genera (see the
# 15B workflow), on cases that reach the checks: the residue fields F_q
# (ScalarFiniteField.h), the maximal orders (Round 2 against the search), the
# Howell forms and the local images of the multiplicity factors (against
# brute force and a listing), the glue arithmetic and the decomposition of
# the primes in Q(zeta_d) (against the explicit formulas). A check that
# fails stops the program with an error, so these are the regular tests,
# run again with the checks.
set -e
set -o pipefail
rm -f CI_CONCLUSION
HP=../../src_genera/GENUS_HermitianPieces
MASS=../../src_genera/GENUS_IsometryMass
if [ "$($HP primes 150 100 2>/dev/null)" != "$(grep -v '^#' hermitian_primes_150_100.txt)" ]; then
    echo "ERROR: the prime decompositions differ under the checks" >&2
    exit 1
fi
echo "--- the prime decompositions, with the checks ---"
grep -v '^#' maximal_order_tests.txt | while read -r f disc; do
    [ -z "$f" ] && continue
    out=$($MASS maximalorder $f 2>/dev/null)
    if [ "${out##* }" != "$disc" ]; then
        echo "ERROR: the maximal order of Z[x] / ($f) gave '$out'" >&2
        exit 1
    fi
done
echo "--- the maximal orders, with the checks ---"
bash RunIsometryMassTests.sh > /dev/null
echo "--- the isometry masses, with the checks ---"
bash RunCharPolyTests.sh sanity_tests.txt | tail -1
rm -f CI_CONCLUSION
echo "All sanity checks passed"
touch CI_CONCLUSION
