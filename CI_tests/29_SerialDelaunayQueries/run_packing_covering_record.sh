#!/bin/bash
# Packing-covering constant of the dimension-5 record, through
# LATT_SerialPeriodicDelaunay.
#
# The record is a 2-periodic NON-LATTICE point set X = Z^5 + {0, c} whose
# packing-covering constant gamma = mu / rho is strictly below
#
#     gamma_5 = sqrt(3/2 + sqrt(13)/6) = 1.4494568681327953...,
#
# the best value over LATTICES, attained by Horvath's Ho_5. So the
# dimension-5 packing-covering problem is not solved by a lattice.
#
# What is checked here is the rational witness: an explicit point set with
# Gram matrix and coset of denominator 10^5. Being rational it is handled
# exactly, and being a legitimate point set in its own right the inequality
# it satisfies is a complete result -- no perturbation argument is needed to
# get back to the optimum. See PackingCoveringRecord/README.md for the
# algebraic optimum itself, which is a different matter.
#
# Usage:  ./run_packing_covering_record.sh [output directory]

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/PackingCoveringRecord"
WORK_DIR="${1:-$(mktemp -d)}"
mkdir -p "$WORK_DIR"

PROG="$SCRIPT_DIR/../../src_delaunay/LATT_SerialPeriodicDelaunay"
if [ ! -x "$PROG" ]; then
    PROG="$(command -v LATT_SerialPeriodicDelaunay || true)"
fi
if [ -z "$PROG" ] || [ ! -x "$PROG" ]; then
    echo "LATT_SerialPeriodicDelaunay not found: build it in src_delaunay first"
    exit 1
fi

NML="$WORK_DIR/record.nml"
OUT="$WORK_DIR/record.out"
ERR="$WORK_DIR/record.err"
cat > "$NML" <<NMLEOF
&SYSTEM
 OutFormat = "PackingCoveringGAP"
 OutFile = "$OUT"
 max_runtime_second = 0
/

&DATA
 arithmetic = "gmp"
 GRAMfile = "$DATA_DIR/Gram_witness"
 FileCosets = "$DATA_DIR/Cosets_witness"
/
NMLEOF

echo "running LATT_SerialPeriodicDelaunay on the rational witness ..."
"$PROG" "$NML" 2> "$ERR"
echo
cat "$OUT"
echo

# gamma^2 is an exact rational, so the comparison with gamma_5 is exact:
#   gamma < gamma_5  <=>  6 gamma^2 - 9 < sqrt(13)  <=>  (6 gamma^2 - 9)^2 < 13,
# both sides being positive.
python3 - "$OUT" <<'PYEOF'
import re, sys
from fractions import Fraction as F
txt = open(sys.argv[1]).read()
def get(key):
    m = re.search(key + r":=\s*(-?\d+(?:/\d+)?)", txt)
    if not m:
        print("could not read " + key + " from the output"); sys.exit(1)
    return F(m.group(1))
nb = int(re.search(r"nb:=(\d+)", txt).group(1))
g2 = get("SquarePackingCoveringConstant")
mu2, lam2 = get("SquareCoveringRadius"), get("SquareMinimalDistance")
lhs = 6*g2 - 9
print(f"orbits of Delaunay cells : {nb}")
print(f"mu^2                     : {float(mu2):.13f}")
print(f"lambda^2                 : {float(lam2):.13f}")
print(f"gamma^2 = 4 mu^2/lambda^2: {float(g2):.13f}")
print(f"gamma                    : {float(g2)**0.5:.13f}")
print(f"gamma_5                  : 1.4494568681328")
print()
print(f"6 gamma^2 - 9     = {lhs}")
print(f"                  = {float(lhs):.12f}  (> 0: {lhs > 0})")
print(f"(6 gamma^2 - 9)^2 = {float(lhs*lhs):.12f}  vs 13")
ok = (g2 == 4*mu2/lam2) and lhs > 0 and lhs*lhs < 13
print()
print("VERDICT: " + ("gamma < gamma_5, PROVED in exact rational arithmetic"
                     if ok else "NOT proved"))
sys.exit(0 if ok else 1)
PYEOF
