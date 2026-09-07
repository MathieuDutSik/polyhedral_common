#!/bin/bash
# Enumerate a few genera and check the class number and the mass certificate.
#
# The three cases are small on purpose so that the test stays quick, but they
# exercise the parts that matter: the p = 2 neighbour construction on an even
# lattice of odd determinant, the deduplication (the second and third cases
# have several classes), and the mass as the stopping criterion.
#
# Expected values come from Hecke: 1, 2 and 3 classes respectively.
set -e
BIN=../../src_genera/GENUS_Enumerate

run_case() {
    name=$1
    expected=$2
    echo "--- $name (expecting $expected classes) ---"
    out=$($BIN gmp ${name}_genus.txt ${name}_lattice.txt ${name}_mass.txt \
              Summary stdout 2>/dev/null)
    echo "$out"
    nb=$(echo "$out" | sed -n 's/^class number = //p')
    complete=$(echo "$out" | sed -n 's/^complete = //p')
    if [ "$nb" != "$expected" ]; then
        echo "ERROR: $name gave $nb classes, expected $expected" >&2
        exit 1
    fi
    if [ "$complete" != "true" ]; then
        echo "ERROR: $name did not reach the mass, so the enumeration is not" \
             "certified complete" >&2
        exit 1
    fi
}

run_case E8 1
run_case rank8_det9 2
run_case rank8_det25 3

# A negative case. With a target mass that cannot be reached -- which is what
# a missing spinor genus, or an unsuitable prime, looks like -- the program
# must terminate and report the enumeration as INCOMPLETE. Reporting
# completeness here, or looping, would mean the mass certificate is not doing
# its job, and every result of the program rests on that certificate.
echo "--- E8 with an unreachable mass (expecting complete = false) ---"
out=$($BIN gmp E8_genus.txt E8_lattice.txt E8_wrongmass.txt Summary stdout \
          2>/dev/null)
echo "$out"
complete=$(echo "$out" | sed -n 's/^complete = //p')
if [ "$complete" != "false" ]; then
    echo "ERROR: an unreachable mass was reported as complete=$complete" >&2
    exit 1
fi

echo "All genus enumeration tests passed"
touch CI_CONCLUSION
