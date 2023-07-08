#!/bin/bash

rm -f Result
rm -f log_*
rm -rf Saving_Bank_nproc*_rank*
rm -rf Saving_Polyhedral_nproc*_rank*

mpirun -np 2 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_canonic.nml
mpirun -np 2 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_canonic_initial_triv.nml

