G(3,1,3) P025 RealAlgebraic reproducibility test
=================================================

Purpose
-------
This package provides a small reproducibility case showing different behavior
between a previous working RealAlgebraic build of POLY_SerialDualDesc and a
current upstream-master build.

This is reported as a reproducible behavior difference only.  It does not
establish that any particular upstream commit is the cause.

Test case
---------
Group: G(3,1,3)
Path: P025
Parameter: t = 7/8

Exact input files
-----------------
input/G3_1_3_path_G313-P025_k7_d8_orbit.ext
input/G3_1_3_path_G313-P025_k7_d8.grp

The original namelist is preserved as:

input/input_G3_1_3_path_G313-P025_k7_d8.original.nml

A portable copy with package-local paths is:

input/input_G3_1_3_path_G313-P025_k7_d8.nml

Observed behavior
-----------------

previous_working_build/

  5/5 normal terminations.

  Incidence type:
  G313-N03

  Same incidence certificate in all five runs:
  d9768bf30e8d1cba7969d0bf1f760569894689fcbdb6414b41f022f0eaecacf4


current_upstream_master/

  0/5 normal terminations.

  Each run reports:

  Failed to find an approximant that allows to conclude,
  please produce better approximants

Reproduction
------------

From this top-level directory:

  /path/to/POLY_SerialDualDesc \
      input/input_G3_1_3_path_G313-P025_k7_d8.nml

Contents
--------

previous_working_build/
  build_info.txt
  run_01.log
  run_02.log
  run_03.log
  run_04.log
  run_05.log

current_upstream_master/
  build_info.txt
  run_01.log
  run_02.log
  run_03.log
  run_04.log
  run_05.log

SHA256SUMS.txt contains package checksums.
