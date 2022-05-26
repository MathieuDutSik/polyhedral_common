// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "CtypeMPI_types.h"
#include "Namelist.h"
#include "NumberTheory.h"
#include <boost/mpi.hpp>
#include <netcdf>
#include <numeric>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  int irank_i;
  MPI_Comm_rank(MPI_COMM_WORLD, &irank_i);
  size_t irank = irank_i;
  //
  if (argc != 3) {
    std::cerr << "Number of argument is = " << argc << "\n";
    std::cerr << "This program is used as\n";
    std::cerr << "CTYP_PrepareAdjacencyFile [WORK_] [WORK_ADJ_]\n";
    return -1;
  }
  std::string PrefixFile = argv[1];
  std::string PrefixAdjFile = argv[2];
  //
  // Creating the netcdf output files.
  //
  std::string eFile = PrefixFile + std::to_string(irank) + ".nc";
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read, netCDF::NcFile::nc4);
  netCDF::NcVar varNbAdj = dataFile.getVar("nb_adjacent");
  size_t n_ctype = varNbAdj.getDim(0).getSize();
  std::cerr << "n_ctype=" << n_ctype << "\n";
  std::vector<int> ListNbAdjacent(n_ctype);
  std::vector<size_t> start{0};
  std::vector<size_t> count{n_ctype};
  varNbAdj.getVar(start, count, ListNbAdjacent.data());
  //
  int TotalNbAdjacencies =
      std::accumulate(ListNbAdjacent.begin(), ListNbAdjacent.end(), 0);
  std::cerr << "TotalNbAdjacencies = " << TotalNbAdjacencies << "\n";
  //
  std::string eFileAdj = PrefixAdjFile + std::to_string(irank) + ".nc";
  netCDF::NcFile dataFileAdj(eFileAdj, netCDF::NcFile::replace,
                             netCDF::NcFile::nc4);
  std::cerr << "step 1\n";
  //
  netCDF::NcDim eDimNbCtype = dataFileAdj.addDim("number_ctype", n_ctype);
  std::cerr << "step 2\n";
  netCDF::NcDim eDimNbAdjacencies =
      dataFileAdj.addDim("number_adjacencies", TotalNbAdjacencies);
  std::cerr << "step 3\n";
  //
  std::vector<std::string> LDimA{"number_ctype"};
  std::vector<std::string> LDimB{"number_adjacencies"};
  //
  netCDF::NcVar varStatus = dataFileAdj.addVar("status", "byte", LDimA);
  varStatus.putAtt("long_name", "status of the C-types");
  std::vector<int8_t> ListStatusAdj(n_ctype, 0);
  varStatus.putVar(start, count, ListStatusAdj.data());
  std::cerr << "step 4\n";
  //
  netCDF::NcVar varIdxProc = dataFileAdj.addVar("idx_proc", "byte", LDimB);
  varIdxProc.putAtt("long_name", "processor of the adjacent C-type");
  std::cerr << "step 5\n";
  //
  netCDF::NcVar varIdxAdj = dataFileAdj.addVar("idx_adjacent", "int", LDimB);
  varIdxAdj.putAtt("long_name", "index of the adjacent C-type");
  std::cerr << "step 6\n";
  //
  MPI_Finalize();
}
