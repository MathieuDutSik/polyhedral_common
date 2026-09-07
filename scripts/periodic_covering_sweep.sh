#!/bin/bash
# Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#
# Sweep random periodic point sets of a given dimension, running the
# random-walk covering search on each, until one beats the record or the
# configurations run out.
#
# Each configuration is drawn by PERIODIC_RandomCosets and searched by
# PERIODIC_LookForRecordCovering; see doc/COVERING_OPTIMIZATION.md. Neither
# failing to find a record nor a configuration that admits no admissible
# coset set stops the sweep: both are normal and are logged.
#
# Usage:
#   scripts/periodic_covering_sweep.sh -n 5 -o /tmp/sweep5 [options]
#
#   -n DIM        the dimension (required)
#   -o OUTDIR     where the per-configuration files go (required)
#   -k COUNT      how many configurations to try (default 20)
#   -t SECONDS    search budget per configuration (default 600)
#   -N "LIST"     denominators to draw from (default "2 3 4")
#   -m "LIST"     coset counts to draw from (default "2 3 4")
#   -r RECORD     the density to beat, or "auto" (default auto)
#   -w STEPS      random jumps used to escape a local minimum (default 20)
#   -b DIR        the directory holding the binaries (default src_delaunay)
#
# The summary of every configuration is appended to OUTDIR/summary.txt and
# printed.
#
# Exit status: 0 when a record was found, 1 when the sweep completed without
# finding one, 2 on a usage or setup error. A status of 1 is a normal
# outcome, not a failure: in dimensions 3 to 5 the best covering is expected
# to be the lattice one, so no periodic configuration should beat it.

set -u

DIM=""
OUTDIR=""
COUNT=20
BUDGET=600
NLIST="2 3 4"
MLIST="2 3 4"
RECORD="auto"
WALK=20
BINDIR=""

while getopts "n:o:k:t:N:m:r:w:b:" opt; do
  case $opt in
    n) DIM=$OPTARG ;;
    o) OUTDIR=$OPTARG ;;
    k) COUNT=$OPTARG ;;
    t) BUDGET=$OPTARG ;;
    N) NLIST=$OPTARG ;;
    m) MLIST=$OPTARG ;;
    r) RECORD=$OPTARG ;;
    w) WALK=$OPTARG ;;
    b) BINDIR=$OPTARG ;;
    *) echo "see the header of $0 for the usage" >&2; exit 2 ;;
  esac
done

if [ -z "$DIM" ] || [ -z "$OUTDIR" ]; then
  echo "the dimension (-n) and the output directory (-o) are required" >&2
  exit 2
fi

if [ -z "$BINDIR" ]; then
  BINDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../src_delaunay" && pwd)"
fi
GEN="$BINDIR/PERIODIC_RandomCosets"
SEARCH="$BINDIR/PERIODIC_LookForRecordCovering"
for prog in "$GEN" "$SEARCH"; do
  if [ ! -x "$prog" ]; then
    echo "missing binary $prog; build it with 'make -C src_delaunay -f Makefile_periodic'" >&2
    exit 2
  fi
done

mkdir -p "$OUTDIR"
SUMMARY="$OUTDIR/summary.txt"
NARR=($NLIST)
MARR=($MLIST)

found=1
for ((i = 0; i < COUNT; i++)); do
  N=${NARR[$((RANDOM % ${#NARR[@]}))]}
  M=${MARR[$((RANDOM % ${#MARR[@]}))]}
  tag="cfg${i}_N${N}_m${M}"
  cosets="$OUTDIR/$tag.cosets"
  # A given (N, m) may admit no admissible configuration at all, which is
  # not a failure of the sweep.
  if ! "$GEN" "$DIM" "$N" "$M" 500 "$cosets" > "$OUTDIR/$tag.gen.log" 2>&1; then
    echo "$tag: no admissible coset configuration" | tee -a "$SUMMARY"
    continue
  fi
  nml="$OUTDIR/$tag.nml"
  out="$OUTDIR/$tag.result.g"
  cat > "$nml" <<EOF
&SYSTEM
 max_runtime_second = $BUDGET
 ApplyStdUnitbuf = T
 Prefix = "$OUTDIR/$tag.hits/"
 OutFile = "$out"
/
&DATA
 arithmetic = "gmp"
 FileDualDescription = "unset"
 FileCosets = "$cosets"
/
&SEARCH
 RecordToBeat = "$RECORD"
 n_walk_steps = $WALK
/
&TSPACE
 TypeTspace = "Classic"
 ClassicDim = $DIM
/
EOF
  if ! "$SEARCH" "$nml" > "$OUTDIR/$tag.search.log" 2>&1; then
    echo "$tag: the search failed, see $OUTDIR/$tag.search.log" | tee -a "$SUMMARY"
    continue
  fi
  line=$(grep -o 'found_record:=[a-z]*\|best_density:=[0-9.e+-]*' "$out" | tr '\n' ' ')
  echo "$tag: $line" | tee -a "$SUMMARY"
  if grep -q 'found_record:=true' "$out"; then
    echo "RECORD FOUND with $cosets, see $out" | tee -a "$SUMMARY"
    found=0
    break
  fi
done

if [ $found -ne 0 ]; then
  echo "no record found over $COUNT configurations; the best density of each" \
       "is in $SUMMARY" | tee -a "$SUMMARY"
fi
exit $found
