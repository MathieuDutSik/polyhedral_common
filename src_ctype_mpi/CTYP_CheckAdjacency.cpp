#include "CtypeMPI_types.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include <boost/mpi.hpp>

#include <netcdf>

//#define DEBUG_INFO



int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int irank_i, n_pes_i;
  MPI_Comm_size(MPI_COMM_WORLD, &n_pes_i);
  MPI_Comm_rank(MPI_COMM_WORLD,&irank_i);
  size_t irank=irank_i;
  size_t n_pes=n_pes_i;
  //
  std::string Prefix = argv[1];
  std::string PrefixAdj = argv[1];
  std::vector<int> ListAdj_Glob;
  std::vector<int> ListNbCtype(n_pes);
  std::vector<int> ListNbAdj_Glob;
  // Reading number of adjacencies
  for (size_t i_pes=0; i_pes<n_pes; i_pes++) {
    std::string eFile = Prefix + std::to_string(i_pes) + ".nc";
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
    netCDF::NcVar varNbAdj=dataFile.getVar("idx_proc");
    size_t n_ctype = varNbAdj.getDim(0).getSize();
    ListNbCtype[i_pes] = n_ctype;
    std::vector<int> ListNbAdj_Loc(n_ctype);
    std::vector<size_t> start{0};
    std::vector<size_t> count{n_ctype};
    varNbAdj.getVar(start, count, ListNbAdj_Loc.data());
    for (auto & eVal : ListNbAdj_Loc)
      ListNbAdj_Glob.push_back(eVal);
  }
  int n_ctype_glob = std::accumulate(ListNbCtype.begin(), ListNbCtype.end(), 0);
  std::vector<int> ListNbAdjShift_Glob(n_ctype_glob);
  ListNbAdjShift_Glob[0] = 0;
  for (int i=1; i<n_ctype_glob; i++)
    ListNbAdjShift_Glob[i] = ListNbAdjShift_Glob[i-1] + ListNbAdj_Glob[i-1];
  std::vector<int> ListShiftNbCtype(n_pes,0);
  for (size_t i_pes=1; i_pes<n_pes; i_pes++)
    ListShiftNbCtype[i_pes] = ListShiftNbCtype[i_pes-1] + ListNbCtype[i_pes-1];
  // Reading the adjacencies themselves
  for (size_t i_pes=0; i_pes<n_pes; i_pes++) {
    std::string eFileAdj = PrefixAdj + std::to_string(i_pes) + ".nc";
    netCDF::NcFile dataFileAdj(eFileAdj, netCDF::NcFile::read);
    netCDF::NcVar varIdxProc=dataFileAdj.getVar("idx_proc");
    netCDF::NcVar varIdxAdjacent=dataFileAdj.getVar("idx_adjacent");
    size_t n_adjacencies = varIdxProc.getDim(0).getSize();
    std::vector<int8_t> ListIdxProc_Loc(n_adjacencies);
    std::vector<int8_t> ListIdxAdjacent_Loc(n_adjacencies);
    for (size_t i=0; i<n_adjacencies; i++) {
      int8_t iProc = ListIdxProc_Loc[i];
      int iAdj_Loc = ListIdxAdjacent_Loc[i];
      int iAdj_Glob = ListShiftNbCtype[iProc] + iAdj_Loc;
      ListAdj_Glob.push_back(iAdj_Glob);
    }
  }
  // Now checking the adjacency is coherent
  auto IsVerticesAdjacent=[&](int const& pos1, int const& pos2) -> bool {
    int nbadj = ListNbAdj_Glob[pos1];
    int shift = ListNbAdjShift_Glob[pos1];
    for (int iadj=0; iadj<nbadj; iadj++) {
      int eAdj = ListAdj_Glob[shift + iadj];
      if (eAdj == pos2)
        return true;
    }
    return false;
  };
  int nb_error=0;
  for (size_t i_loc=0; i_loc<ListNbCtype[irank]; i_loc++) {
    int i_glob = ListShiftNbCtype[irank] + i_loc;
    int nbadj = ListNbAdj_Glob[i_glob];
    int shift = ListNbAdjShift_Glob[i_glob];
    for (int iadj=0; iadj<nbadj; iadj++) {
      int eAdj = ListAdj_Glob[shift + iadj];
      if (!IsVerticesAdjacent(eAdj, i_glob))
        nb_error++;
    }
  }
  std::cerr << "nb_error=" << nb_error << "\n";

  MPI_Finalize();
}
