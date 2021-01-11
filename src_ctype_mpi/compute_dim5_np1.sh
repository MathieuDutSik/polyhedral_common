#!/bin/bash
echo "1 : CTYP_PrepareInitialFile"
./CTYP_PrepareInitialFile create_netcdf5_np1.nml

echo "        -----"
echo "2 : CTYP_MPI_Enumeration_c"
./CTYP_MPI_Enumeration_c ctype_enum5.nml

echo "        -----"
echo "3 : CTYP_ComputeInvariant"
./CTYP_ComputeInvariant WORK_5_ WORK_OUT_5_

echo "        -----"
echo "4 : CTYP_PrepareAdjacencyFile"
./CTYP_PrepareAdjacencyFile WORK_OUT_5_ WORK_ADJ_5_

echo "        -----"
echo "5 : CTYP_MPI_EnumerationAdjacencies"
./CTYP_MPI_EnumerationAdjacencies ctype_adj_enum5.nml

