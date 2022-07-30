#!/bin/bash

rm -rf Saving_Bank*
rm -rf Saving_Polyhedral*


mpirun -np 1 ../../src_dualdesc/POLY_MPI_DualDesc Main_mpi.nml
