#!/bin/bash

rm -rf Saving_Bank_nproc*
rm -rf Saving_Polyhedral_nproc*

mkdir -p RUN0
mkdir -p RUN1
mkdir -p RUN2

mpirun --oversubscribe -np 10 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_direct.nml
mv log_* RUN0/

mpirun --oversubscribe -np 10 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_inittriv.nml
mv log_* RUN1/

mpirun --oversubscribe -np 10 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_guess.nml
mv log_* RUN2/

