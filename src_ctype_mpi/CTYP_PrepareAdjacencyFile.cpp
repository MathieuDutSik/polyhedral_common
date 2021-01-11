#include "CtypeMPI_types.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include <netcdf>
#include <numeric>

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int irank_i;
  MPI_Comm_rank(MPI_COMM_WORLD,&irank_i);
  size_t irank=irank_i;
  //
  if (argc != 3) {
    std::cerr << "Number of argument is = " << argc << "\n";
    std::cerr << "This program is used as\n";
    std::cerr << "CTYP_PrepareAdjacencyFile [WORK_] [WORK_ADJ_]\n";
    return -1;
  }
  std::string PrefixFile=argv[1];
  std::string PrefixAdjFile=argv[2];
  //
  // Creating the netcdf output files.
  //
  std::string eFile=PrefixFile + std::to_string(irank) + ".nc";
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read, netCDF::NcFile::nc4);
  netCDF::NcVar varNbAdj=dataFile.getVar("nb_adjacent");
  size_t n_ctype = varNbAdj.getDim(0).getSize();
  std::vector<int> ListNbAdjacent(n_ctype);
  std::vector<size_t> start{0};
  std::vector<size_t> count{n_ctype};
  varNbAdj.getVar(start, count, ListNbAdjacent.data());
  //
  int TotalNbAdjacencies = std::accumulate(ListNbAdjacent.begin(), ListNbAdjacent.end(), 0);
  //
  std::string eFileAdj=PrefixAdjFile + std::to_string(irank) + ".nc";
  netCDF::NcFile dataFileAdj(eFileAdj, netCDF::NcFile::replace, netCDF::NcFile::nc4);
  //
  netCDF::NcDim eDimNbCtype=dataFileAdj.addDim("number_ctype", n_ctype);
  netCDF::NcDim eDimNbAdjacencies=dataFileAdj.addDim("number_adjacencies", TotalNbAdjacencies);
  //
  std::vector<std::string> LDimA{"number_ctype"};
  std::vector<std::string> LDimB{"number_adjacencies"};
  //
  netCDF::NcVar varStatus = dataFileAdj.addVar("status", "byte", LDimA);
  varNbAdj.putAtt("long_name", "status of the C-types");
  std::vector<int8_t> ListStatusAdj(n_ctype, 0);
  varStatus.putVar(start, count, ListStatusAdj.data());
  //
  netCDF::NcVar varIdxProc = dataFileAdj.addVar("idx_proc", "byte", LDimB);
  varNbAdj.putAtt("long_name", "status of the C-types");
  //
  netCDF::NcVar varIdxAdj = dataFileAdj.addVar("idx_adjacent", "int", LDimB);
  varNbAdj.putAtt("long_name", "status of the C-types");
  //
  MPI_Finalize();
}
