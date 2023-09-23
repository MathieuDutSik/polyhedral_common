// clang-format off
#include "NumberTheory.h"
#include "CtypeMPI_types.h"
#include "Namelist.h"
#include "hash_functions.h"
#include "rational.h"
#include <unordered_map>
// clang-format on



template<typename T>
TypeCtypeExch<T> RandomWalk(TypeCtypeExch<T> const& eMat) {
  using Tidx = int32_t;
  TypeCtypeExch<T> WorkT = eMat;
  for (int iter=0; iter<50; iter++) {
    bool canonicalize = false;
    std::vector<TypeCtypeExch<T>> ListAdj = CTYP_Kernel_GetAdjacentCanonicCtypes<T,Tidx>(WorkT, canonicalize);
    int n_adj = ListAdj.size();
    int pos = random() % n_adj;
    WorkT = ListAdj[pos];
  }
  return WorkT;
}

template<typename T>
void GetLocalFreenessMinimum(int const& dim, int const& max_s, std::string const& Prefix) {
  std::cerr << "GetLocalFreenessMinimum, step 1\n";
  MyMatrix<T> M = GetPrincipalDomain<T>(dim);
  std::cerr << "M=\n";
  WriteMatrix(std::cerr, M);
  std::cerr << "GetLocalFreenessMinimum, step 2\n";
  TypeCtypeExch<T> TheCtypeArr{M};
  int curr_nb_free = dim * (dim + 1) / 2;
  int iter1 = 0;
  int iter2 = 0;
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
    if (the_min >= curr_nb_free && curr_nb_free > 0) {
      TheCtypeArr = RandomWalk(TheCtypeArr);
      curr_nb_free = CTYP_GetNumberFreeVectors(TheCtypeArr);
      std::cerr << "iter1=" << iter1 << " iter2=" << iter2 << " After RandomWalk curr_nb_free=" << curr_nb_free << "\n";
      iter1++;
      iter2 = 0;
    } else {
      curr_nb_free = the_min;
      std::vector<size_t> ListIdx;
      for (size_t u=0; u<ListNbFree.size(); u++) {
        if (ListNbFree[u] == the_min) {
          ListIdx.push_back(u);
        }
      }
      int n_min = ListIdx.size();
      int pos = random() % n_min;
      std::cerr << "iter1=" << iter1 << " iter2=" << iter2 << " curr_nb_free=" << curr_nb_free << " n_min=" << n_min << " pos=" << pos << "\n";
      TheCtypeArr = ListAdj[ListIdx[pos]];
      if (curr_nb_free <= max_s) {
        std::string FileOut = FindAvailableFileFromPrefix(Prefix);
        WriteMatrixFile(FileOut, TheCtypeArr.eMat);
      }
      if (curr_nb_free == 0) {
        return;
      }
      iter2++;
    }
  }


}




int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    using Tint = int64_t;
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "CTYP_LookForNoFreeVector  dim  n_try  max_s  FileOut\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "dim     : the dimension\n";
      std::cerr << "ntry    : the umber of tries\n";
      std::cerr << "max_s   : the maximum number of free vectors that still renders it interesting\n";
      std::cerr << "FileOut : the file in output\n";
      throw TerminalException{1};
    }
    int dim = ParseScalar<size_t>(argv[1]);
    size_t n_try = ParseScalar<size_t>(argv[2]);
    int max_s = ParseScalar<size_t>(argv[3]);
    std::string Prefix = argv[4];
    for (size_t i_try=0; i_try<n_try; i_try++) {
      GetLocalFreenessMinimum<Tint>(dim, max_s, Prefix);
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in CTYP_LookForNoFreeVector\n";
    exit(e.eVal);
  }
  runtime(time1);
}

