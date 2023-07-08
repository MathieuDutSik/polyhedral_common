#!/bin/bash

rm -f Result
rm -f log_*
rm -rf Saving_Bank_nproc*_rank*
rm -rf Saving_Polyhedral_nproc*_rank*

mpirun -np 2 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_canonic.nml
mv log_2_0 log_blockA_0
mv log_2_1 log_blockA_1


mpirun -np 2 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_canonic_initial_triv.nml
mv log_2_0 log_blockB_0
mv log_2_1 log_blockB_1
