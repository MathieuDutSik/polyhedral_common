// clang-format off
#include "NumberTheory.h"
#include "CtypeMPI_types.h"
#include "Namelist.h"
#include "hash_functions.h"
#include "rational.h"
#include "sparse-map/include/tsl/sparse_map.h"
#include <boost/mpi.hpp>
#include <netcdf>
#include <unordered_map>
// clang-format on

template<typename T>
std::pair<TypeCtypeExch<T>,int> GetLocalFreenessMinimum(int const& dim) {
  std::cerr << "GetLocalFreenessMinimum, step 1\n";
  MyMatrix<T> M = GetPrincipalDomain<T>(dim);
  std::cerr << "M=\n";
  WriteMatrix(std::cerr, M);
  std::cerr << "GetLocalFreenessMinimum, step 2\n";
  TypeCtypeExch<T> TheCtypeArr{M};
  int curr_nb_free = dim * (dim + 1) / 2;
  int iter = 0;
  std::cerr << "GetLocalFreenessMinimum, step 3\n";
  using Tidx = int32_t;
  while(true) {
    std::cerr << "GetLocalFreenessMinimum, Before generation of adjacent domains\n";
    bool canonicalize = false;
    std::vector<TypeCtypeExch<T>> ListAdj = CTYP_Kernel_GetAdjacentCanonicCtypes<T,Tidx>(TheCtypeArr, canonicalize);
    std::vector<int> ListNbFree;
    for (auto & eAdj : ListAdj) {
      int nb_free = CTYP_GetNumberFreeVectors(eAdj);
      //      std::cerr << " nb_free=" << nb_free << "\n";
      ListNbFree.push_back(nb_free);
    }
    int the_min = VectorMin(ListNbFree);
    if (the_min >= curr_nb_free) {
      return {TheCtypeArr,curr_nb_free};
    }
    curr_nb_free = the_min;
    std::vector<TypeCtypeExch<T>> ListMin;
    for (size_t u=0; u<ListNbFree.size(); u++) {
      if (ListNbFree[u] == the_min) {
        ListMin.push_back(ListAdj[u]);
      }
    }
    int n_min = ListMin.size();
    int pos = random() % n_min;
    std::cerr << "iter=" << iter << " curr_nb_free=" << curr_nb_free << " n_min=" << n_min << " pos=" << pos << "\n";
    TheCtypeArr = ListAdj[pos];
    iter++;
  }


}




int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using Tint = int64_t;
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CTYP_LookForNoFreeVector dim n_try FileOut\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "dim     : the dimension\n";
      std::cerr << "ntry    : the umber of tries\n";
      std::cerr << "FileOut : the file in output\n";
      throw TerminalException{1};
    }
    int dim = ParseScalar<size_t>(argv[1]);
    size_t n_try = ParseScalar<size_t>(argv[2]);
    std::string FileOut = argv[3];
    std::unordered_set<MyMatrix<Tint>> set_mat;
    std::map<int, size_t> list_mult;
    for (size_t i_try=0; i_try<n_try; i_try++) {
      std::pair<TypeCtypeExch<Tint>,int> pair = GetLocalFreenessMinimum<Tint>(dim);
      MyMatrix<Tint> CanMat = LinPolytopeAntipodalIntegral_CanonicForm<Tint>(pair.first.eMat);
      set_mat.insert(CanMat);
      list_mult[pair.second]++;
    }
    std::cerr << "list_mult =";
    for (auto & kv : list_mult)
      std::cerr << " (" << kv.first << "/" << kv.second << ")";
    std::cerr << "\n";
    //
    std::ofstream os(FileOut);
    size_t n_sol = set_mat.size();
    os << n_sol << "\n";
    for (auto & eMat : set_mat) {
      WriteMatrix(os, eMat);
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_lrs\n";
    exit(e.eVal);
  }
  runtime(time1);
  
}

