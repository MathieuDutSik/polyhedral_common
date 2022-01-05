#include "NumberTheory.h"
#include "CtypeMPI_types.h"
#include "Namelist.h"
#include <boost/mpi.hpp>

#include <netcdf>

//#define DEBUG_INFO



int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);
  int irank_i;
  MPI_Comm_rank(MPI_COMM_WORLD,&irank_i);
  size_t irank=irank_i;
  //
  //  int n_pes_i;
  //  MPI_Comm_size(MPI_COMM_WORLD, &n_pes_i);
  //  size_t n_pes=n_pes_i;
  //
  if (argc != 3) {
    std::cerr << "Usage:\n";
    std::cerr << "CTYP_SearchMartrix [WORK_OUT_] [FileI]\n";
    return 0;
  }
  std::string Prefix = argv[1];
  std::string FileI = argv[2];
  std::cerr << "Prefix=" << Prefix << " FileI=" << FileI << "\n";
  using T = int8_t;
  std::ifstream is(FileI);
  MyMatrix<T> M = ReadMatrix<T>(is);
  //
  std::string eFile = Prefix + std::to_string(irank) + ".nc";
  std::cerr << "eFile=" << eFile << "\n";
  netCDF::NcFile dataFile(eFile, netCDF::NcFile::read);
  netCDF::NcVar varNbAdj=dataFile.getVar("nb_adjacent");
  netCDF::NcVar varCtype=dataFile.getVar("Ctype");
  size_t n_ctype = varCtype.getDim(0).getSize();
  size_t n_vect  = varCtype.getDim(1).getSize();
  size_t n       = varCtype.getDim(2).getSize();
  //
  std::cerr << "n_ctype=" << n_ctype << " n_vect=" << n_vect << " n=" << n << "\n";
  std::vector<T> V(n_vect, n);
  auto IsMatching=[&]() -> bool {
    int idx=0;
    for (size_t i_vect=0; i_vect<n_vect; i_vect++)
      for (size_t i=0; i<n; i++) {
        if (M(i_vect, i) != V[idx])
          return false;
      }
    return true;
  };
  int64_t nb_found_loc = 0;
  for (size_t ic=0; ic<n_ctype; ic++) {
    std::vector<size_t> start{ic, 0, 0};
    std::vector<size_t> count{1, n_vect, n};
    varCtype.getVar(start, count, V.data());
    if (IsMatching())
      nb_found_loc++;
  }
  int64_t nb_found_tot = 0;
  int ierr1 = MPI_Allreduce(&nb_found_loc, &nb_found_tot, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
  if (ierr1 != MPI_SUCCESS) {
    std::cerr << "ierr1 wrongly set\n";
    throw TerminalException{1};
  }
  std::cerr << "irank=" << irank << " nb_found_loc=" << nb_found_loc << " nb_found_tot=" << nb_found_tot << "\n";
  MPI_Finalize();
}
