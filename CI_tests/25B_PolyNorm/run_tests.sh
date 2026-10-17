#!/bin/bash
# Packing and covering scalars of polytopes with respect to Z^n: the
# branch and bound covering radius is compared with the brute force
# enumeration of Cslovjecsek, Malikiosis, Naszodi and Schymura when that
# enumeration is small enough, and both scalars are compared with their
# known values. Columns: file, alpha, mu, largest brute force size.
set -e
PROG=../../src_polynorm/POLYNORM_TestCovering
if [ ! -x "$PROG" ]; then
  echo "The program $PROG is missing" >&2
  exit 1
fi
rm -f CI_CONCLUSION
while read -r file alpha mu max_bf; do
  case "$file" in
    ''|\#*) continue ;;
  esac
  echo "=== $file (alpha=$alpha mu=$mu)"
  if ! $PROG gmp "$file" "$alpha" "$mu" "$max_bf" > "$file.log" 2>&1; then
    echo "Failure on $file" >&2
    tail -20 "$file.log" >&2
    exit 1
  fi
  grep '^TEST' "$file.log"
done <<'LIST'
# Dimension 2: the brute force is run on every case.
square.ext 1 1 2000000
diamond.ext 1/2 1 2000000
triangle.ext 1 2 2000000
hexagon.ext 1/2 2/3 2000000
pentagon.ext 1/4 1/2 2000000
asym2d.ext 1/2 5/6 2000000
# Dimension 3: the brute force runs on the cube (4 million systems).
cube3.ext 1 1 5000000
octa3.ext 1/2 3/2 2000000
simplex3.ext 1 3 2000000
asym3d.ext 1/2 4/3 2000000
frac3d.ext 2/3 12/5 2000000
# Dimension 4.
cube4.ext 1 1 2000000
cross4.ext 1/2 2 2000000
simplex4.ext 1 4 2000000
cell24.ext 1/2 1 2000000
LIST

# Optimal lattices. The densest lattice packings of the tetrahedron (18/49,
# Hoylman) and of the octahedron (18/19, Minkowski), as given by Betke and
# Henk, and the thinnest lattice covering of the octahedron (9/8 with the
# body-centered cubic lattice, Lian and Xue) are evaluated directly, and
# must be certified strict local optima.
PROG_CERT=../../src_polynorm/POLYNORM_LatticeCertify
PROG_OPT=../../src_polynorm/POLYNORM_OptimalLattice
for prog in "$PROG_CERT" "$PROG_OPT"; do
  if [ ! -x "$prog" ]; then
    echo "The program $prog is missing" >&2
    exit 1
  fi
done
check_line() {
  # $1 log file, $2 the expected line (fixed string)
  if ! grep -qF -- "$2" "$1"; then
    echo "Failure: expected the line '$2' in $1" >&2
    cat "$1" >&2
    exit 1
  fi
}
echo "=== certification of known optimal lattices"
$PROG_CERT gmp tetra_betke_henk.ext tetra_betke_henk.lat stdout > cert_tetra.log 2>/dev/null
check_line cert_tetra.log "Packing: alpha=1 density=18/49"
check_line cert_tetra.log "certified_local_optimum=1 (strict local optimum on every branch)"
$PROG_CERT gmp octa3.ext octa3_betke_henk.lat stdout > cert_octa.log 2>/dev/null
check_line cert_octa.log "Packing: alpha=1 density=18/19"
check_line cert_octa.log "certified_local_optimum=1 (strict local optimum on every branch)"
$PROG_CERT gmp octa3.ext octa3_bcc.lat stdout > cert_bcc.log 2>/dev/null
check_line cert_bcc.log "Covering: mu=3/2 density=9/8"
check_line cert_bcc.log "certified_local_optimum=1 (strict local optimum)"
echo "TEST certification: OK"

# The search from a few random starts (the generator is portable, so the
# seed fixes the run) must find the known optimum exactly, polished and
# certified. Columns: file, problem, number of starts, density.
while read -r file problem nstart dens; do
  case "$file" in
    ''|\#*) continue ;;
  esac
  echo "=== $file $problem (density=$dens)"
  log="opt_${file%.ext}_$problem.log"
  if ! $PROG_OPT gmp "$file" "$problem" "$nstart" 1 none stdout > "$log" 2>/dev/null; then
    echo "Failure of POLYNORM_OptimalLattice on $file $problem" >&2
    exit 1
  fi
  best=$(sed -n '/^Best lattice/,$p' "$log")
  if ! echo "$best" | grep -qF "exact=$dens"; then
    echo "Failure on $file $problem: the best density is not $dens" >&2
    cat "$log" >&2
    exit 1
  fi
  if ! echo "$best" | grep -qF "certified_local_optimum=1"; then
    echo "Failure on $file $problem: the optimum is not certified" >&2
    cat "$log" >&2
    exit 1
  fi
  echo "TEST $file $problem: OK"
done <<'LIST'
# Dimension 2: Fary (triangle), Sriamorn (quadrilateral). The packing
# problem is also solved by the exact algorithm.
triangle.ext packing 2 2/3
asym2d.ext packing 2 3/4
triangle.ext covering 2 3/2
asym2d.ext covering 2 4/3
# Dimension 3: Betke and Henk.
simplex3.ext packing 3 18/49
octa3.ext packing 3 18/19
LIST

# The algorithm of Betke and Henk: densest lattice packings of 3-polytopes,
# over the rationals and the fields Q(sqrt 2), Q(sqrt 5), against the
# values of their paper. Columns: arith, file, density (12 digits).
PROG_BH=../../src_polynorm/POLYNORM_BetkeHenk
if [ ! -x "$PROG_BH" ]; then
  echo "The program $PROG_BH is missing" >&2
  exit 1
fi
while read -r arith file dens; do
  case "$arith" in
    ''|\#*) continue ;;
  esac
  echo "=== Betke-Henk $file ($arith, density $dens)"
  log="bh_${file%.ext}.log"
  if ! $PROG_BH "$arith" "$file" stdout 4 > "$log" 2>/dev/null; then
    echo "Failure of POLYNORM_BetkeHenk on $file" >&2
    exit 1
  fi
  check_line "$log" "Densest lattice packing density ~ $dens"
  check_line "$log" "Every candidate of smaller determinant is resolved"
  echo "TEST Betke-Henk $file: OK"
done <<'LIST'
# Octahedron 18/19, truncated cube 9/(5 + 3 sqrt 2), dodecahedron (2 + tau)/4.
gmp octa3.ext 0.947368421053
Qsqrt2 truncated_cube_qsqrt2.ext 0.973747687989
Qsqrt5 dodecahedron_qsqrt5.ext 0.904508497187
LIST
# The maximal kissing number: the number of touching translates of the
# densest-contact lattice packing. The square and the cube are
# parallelotopes and reach 3^n - 1; the triangle and the hexagon reach 6,
# the most a plane convex body other than a parallelogram can.
# Columns: file, kissing number.
PROG_KISS=../../src_polynorm/POLYNORM_MaxKissing
if [ ! -x "$PROG_KISS" ]; then
  echo "The program $PROG_KISS is missing" >&2
  exit 1
fi
while read -r file kiss; do
  case "$file" in
    ''|\#*) continue ;;
  esac
  echo "=== kissing $file (kiss=$kiss)"
  log="kiss_${file%.ext}.log"
  if ! $PROG_KISS gmp "$file" text "$log" 2>/dev/null; then
    echo "Failure of POLYNORM_MaxKissing on $file" >&2
    exit 1
  fi
  check_line "$log" "Maximal kissing number kiss=$kiss"
  if grep -q WARNING "$log"; then
    echo "Failure on $file: undecided configurations remain" >&2
    cat "$log" >&2
    exit 1
  fi
  echo "TEST kissing $file: OK"
done <<'LIST'
square.ext 8
triangle.ext 6
hexagon.ext 6
cube3.ext 26
LIST

echo "All tests passed" > CI_CONCLUSION
