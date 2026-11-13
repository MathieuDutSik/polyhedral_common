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

# The namelist of a case: its DATA and METHOD blocks (name.nml) and a SYSTEM
# block sending the Summary to stdout, written to a temporary directory. The
# file names in name.nml are relative to this directory.
TMPD=$(mktemp -d)
trap 'rm -rf "$TMPD"' EXIT
run_nml() {
    { cat $1.nml; printf '\n&SYSTEM\n OutFormat = "Summary"\n OutFile = "stdout"\n/\n'; } \
        > $TMPD/$1.nml
    $BIN $TMPD/$1.nml 2>/dev/null
}

run_case() {
    name=$1
    expected=$2
    echo "--- $name (expecting $expected classes) ---"
    out=$(run_nml $name)
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
# The same genus without FileMass: the mass is computed from the seed by the
# mass formula, and must close at the same 3 classes.
run_case rank8_det25_computedmass 3
# An odd genus, the odd unimodular lattices of rank 12 (Z^12, E8 + Z^4 and
# D12+), enumerated through the genus of their even sublattices.
run_case I12_odd 3

# A negative case. With a target mass that cannot be reached -- which is what
# a missing spinor genus, or an unsuitable prime, looks like -- the program
# must terminate and report the enumeration as INCOMPLETE. Reporting
# completeness here, or looping, would mean the mass certificate is not doing
# its job, and every result of the program rests on that certificate.
echo "--- E8 with an unreachable mass (expecting complete = false) ---"
out=$(run_nml E8_wrongmass)
echo "$out"
complete=$(echo "$out" | sed -n 's/^complete = //p')
if [ "$complete" != "false" ]; then
    echo "ERROR: an unreachable mass was reported as complete=$complete" >&2
    exit 1
fi

echo "All genus enumeration tests passed"
touch CI_CONCLUSION
