#!/bin/bash
# Enumerate a few genera and check the class number and the mass certificate.
#
# The cases are small on purpose so that the test stays quick, but they
# exercise the parts that matter: the p = 2 neighbor construction on an even
# lattice of odd determinant, the deduplication (most cases have several
# classes), the mass computed by the program as the stopping criterion
# (checked against the mass given by Hecke where name_mass.txt exists), the
# seed built from a genus symbol, an odd genus, and the choice of the prime
# for a genus with two spinor genera.
#
# The expected class numbers come from Hecke and DefLatDB.
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
    # The mass is computed by the program; where the mass of the genus is
    # known independently (from Hecke, name_mass.txt as "num den"), the two
    # must agree.
    if [ -f ${name}_mass.txt ]; then
        known=$(tr -s ' \n' ' ' < ${name}_mass.txt | sed 's/ *$//' | sed 's# #/#')
        target=$(echo "$out" | sed -n 's/^target mass = //p')
        if [ "$target" != "$known" ]; then
            echo "ERROR: $name has the computed mass $target, the known one is" \
                 "$known" >&2
            exit 1
        fi
    fi
}

run_case E8 1
run_case rank8_det9 2
run_case rank8_det25 3
# The same genus from its symbol alone: the seed is built from the symbol
# (GenusSymbol.h), and the mass is computed from that seed.
run_case rank8_det25_symbol 3
# An odd genus, the odd unimodular lattices of rank 12 (Z^12, E8 + Z^4 and
# D12+), enumerated through the genus of their even sublattices.
run_case I12_odd 3
# The same odd genus from its symbol alone.
run_case I12_symbol 3
# A genus with two spinor genera (DefLatDB c02371). The default prime of
# out_prime would be 5 anyway, but the program must find from the spinor
# kernels that 5 connects the two spinor genera, and close the mass with
# the single seed built from the symbol.
run_case spinor_two 3

# A negative case. The same genus at the prime 7, whose neighbour graph
# leaves out one of the two spinor genera: the program must terminate and
# report the enumeration as INCOMPLETE. Reporting completeness here, or
# looping, would mean the mass certificate is not doing its job, and every
# result of the program rests on that certificate.
echo "--- spinor_disconnected (expecting complete = false) ---"
out=$(run_nml spinor_disconnected)
echo "$out"
complete=$(echo "$out" | sed -n 's/^complete = //p')
if [ "$complete" != "false" ]; then
    echo "ERROR: a genus missing a spinor genus was reported as" \
         "complete=$complete" >&2
    exit 1
fi

echo "All genus enumeration tests passed"
touch CI_CONCLUSION
