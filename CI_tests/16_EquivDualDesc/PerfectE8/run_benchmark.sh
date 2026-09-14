#!/bin/bash
# Dual description of the perfect domain of E8.
#
# PerfectE8_saving.nml has max_runtime = 1800, so the program stops every half
# hour and flushes the whole recursion tree (Saving_Polyhedral) and the bank
# (Saving_Bank) to disk. Restarting picks up where it left off, so the loop
# below is a checkpointed version of a single long run: the enumeration is over
# when the "orbits" output file appears.

rm -rf Saving_Bank Saving_Polyhedral orbits
mkdir Saving_Bank Saving_Polyhedral

while [ ! -f orbits ]; do
  date
  ../../../src_dualdesc/POLY_SerialDualDesc PerfectE8_saving.nml 2>&1 | tee -a err2_polydualdesc
  rc=${PIPESTATUS[0]}
  if [ $rc -ne 0 ]; then
    echo "POLY_SerialDualDesc failed with rc=$rc" >&2
    exit $rc
  fi
done
