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
  if ! $PROG "$file" "$alpha" "$mu" "$max_bf" > "$file.log" 2>&1; then
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
echo "All tests passed" > CI_CONCLUSION
