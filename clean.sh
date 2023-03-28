#!/bin/bash
set -e
(cd src_poly && make clean)
(cd src_poly && make -f Makefile_double_cddlib clean)
(cd src_polydecomp && make clean)
(cd src_latt && make clean)
(cd src_short && make clean)
(cd src_dualdesc && make clean)
(cd src_sparse_solver && make clean)
(cd src_spectra && make clean)
(cd src_copos && make clean)
#(cd src_ctype_mpi && make clean)
(cd src_perfect && make clean)
(cd src_indefinite && make clean)
(cd src_lorentzian && make clean)
(cd src_perfect_mpi && make clean)
(cd src_poincare_polyhedron && make clean)
