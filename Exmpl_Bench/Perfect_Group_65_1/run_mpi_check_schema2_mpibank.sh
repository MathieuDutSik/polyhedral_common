#!/bin/bash

rm -f Result
rm -f log_*
rm -rf Saving_Bank_nproc*_rank*
rm -rf Saving_Polyhedral_nproc*_rank*

mpirun -np 3 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_mpibank_canonic.nml
mv log_3_0 log_blockA_0
mv log_3_1 log_blockA_1
mv log_3_2 log_blockA_2


mpirun -np 3 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_mpibank_canonic_initial_triv.nml
mv log_3_0 log_blockB_0
mv log_3_1 log_blockB_1
mv log_3_2 log_blockB_2
