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
  if (argc != 3) {
    std::cerr << "Usage:\n";
    std::cerr << "CTYP_CheckAdjacency [WORK_OUT_] [WORK_ADJ]\n";
    return 0;
  }
  std::string Prefix = argv[1];
  std::string PrefixAdj = argv[2];
  std::cerr << "Prefix=" << Prefix << " PrefixAdj=" << PrefixAdj << "\n";
  std::vector<int> ListAdj_Glob;
  std::vector<size_t> ListNbCtype(n_pes);
  std::vector<int> ListNbAdj_Glob;
  std::vector<int> ListNbFree_Glob;
  // Reading number of adjacencies
  for (size_t i_pes=0; i_pes<n_pes; i_pes++) {
    std::string eFile = Prefix + std::to_string(i_pes) + ".nc";
    std::cerr << "eFile=" << eFile << "\n";
    netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
    netCDF::NcVar varNbAdj=dataFile.getVar("nb_adjacent");
    size_t n_ctype = varNbAdj.getDim(0).getSize();
    std::cerr << "n_ctype=" << n_ctype << "\n";
    ListNbCtype[i_pes] = n_ctype;
    std::vector<int> ListV_Loc(n_ctype);
    std::vector<size_t> start{0};
    std::vector<size_t> count{n_ctype};
    varNbAdj.getVar(start, count, ListV_Loc.data());
    for (auto & eVal : ListV_Loc)
      ListNbAdj_Glob.push_back(eVal);
    //
    netCDF::NcVar varNbFree=dataFile.getVar("nb_free");
    varNbFree.getVar(start, count, ListV_Loc.data());
    for (auto & eVal : ListV_Loc)
      ListNbFree_Glob.push_back(eVal);
  }
  std::cerr << "We have ListNbAdj_Glob\n";
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
    std::cerr << "eFileAdj=" << eFileAdj << "\n";
    netCDF::NcFile dataFileAdj(eFileAdj, netCDF::NcFile::read);
    netCDF::NcVar varIdxProc=dataFileAdj.getVar("idx_proc");
    netCDF::NcVar varIdxAdjacent=dataFileAdj.getVar("idx_adjacent");
    size_t n_adjacencies = varIdxProc.getDim(0).getSize();
    std::vector<int8_t> ListIdxProc_Loc(n_adjacencies);
    std::vector<int8_t> ListIdxAdjacent_Loc(n_adjacencies);
    std::vector<size_t> start{0};
    std::vector<size_t> count{n_adjacencies};
    varIdxProc.getVar(start, count, ListIdxProc_Loc.data());
    varIdxAdjacent.getVar(start, count, ListIdxAdjacent_Loc.data());
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
  int nb_error = 0;
  int max_diff_nb_free = 0;
  for (size_t i_loc=0; i_loc<ListNbCtype[irank]; i_loc++) {
    int i_glob = ListShiftNbCtype[irank] + i_loc;
    int nbadj = ListNbAdj_Glob[i_glob];
    int shift = ListNbAdjShift_Glob[i_glob];
    std::cerr << "i_loc=" << i_loc << " i_glob=" << i_glob << " nbadj=" << nbadj << " shift=" << shift << "\n";
    int eFree = ListNbFree_Glob[i_glob];
    for (int iadj=0; iadj<nbadj; iadj++) {
      int eAdj = ListAdj_Glob[shift + iadj];
      int fFree = ListNbFree_Glob[eAdj];
      if (!IsVerticesAdjacent(eAdj, i_glob)) {
        std::cerr << "  eAdj=" << eAdj << " i_glob=" << i_glob << "\n";
        nb_error++;
      }
      int delta_nb_free = T_abs(eFree - fFree);
      max_diff_nb_free = T_max(max_diff_nb_free, delta_nb_free);
    }
  }
  std::cerr << "nb_error=" << nb_error << " max_diff_nb_free=" << max_diff_nb_free << "\n";

  MPI_Finalize();
}
