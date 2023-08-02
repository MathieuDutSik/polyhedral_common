#!/bin/bash

rm -f Result
rm -f log_*
rm -rf Saving_Bank_nproc*_rank*
rm -rf Saving_Polyhedral_nproc*_rank*

mpirun -np 10 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_mpibank_canonic.nml
mv log_10_0 log_blockA_0
mv log_10_1 log_blockA_1
mv log_10_2 log_blockA_2
mv log_10_3 log_blockA_3
mv log_10_4 log_blockA_4
mv log_10_5 log_blockA_5
mv log_10_6 log_blockA_6
mv log_10_7 log_blockA_7
mv log_10_8 log_blockA_8
mv log_10_9 log_blockA_9


mpirun -np 10 ../../src_dualdesc/POLY_MPI_DualDesc Perfect_mpi_mpibank_canonic_initial_triv.nml
mv log_10_0 log_blockA_0
mv log_10_1 log_blockA_1
mv log_10_2 log_blockA_2
mv log_10_3 log_blockA_3
mv log_10_4 log_blockA_4
mv log_10_5 log_blockA_5
mv log_10_6 log_blockA_6
mv log_10_7 log_blockA_7
mv log_10_8 log_blockA_8
mv log_10_9 log_blockA_9
