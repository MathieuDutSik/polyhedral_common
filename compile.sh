#!/bin/bash
set -e
(cd src_group && make clean && make)
(cd src_poly && make clean && make)
(cd src_poly && make -f Makefile_double_cddlib clean && make -f Makefile_double_cddlib)
(cd src_polydecomp && make clean && make)
(cd src_latt && make clean && make)
(cd src_short && make clean && make)
(cd src_dualdesc && make clean && make)
(cd src_sparse_solver && make clean && make)
(cd src_spectra && make clean && make)
(cd src_copos && make clean && make)
(cd src_ctype && make clean && make)
(cd src_isotropy && make clean && make)
(cd src_lorentzian && make clean && make)
(cd src_spherical_code && make clean && make)
(cd src_poincare_polyhedron && make clean && make)
(cd src_indefinite && make clean && make)
echo "Normal termination of compile.sh"

