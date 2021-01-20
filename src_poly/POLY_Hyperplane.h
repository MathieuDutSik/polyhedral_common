#ifndef INCLUDE_POLY_HYPERPLANE
#define INCLUDE_POLY_HYPERPLANE

#include "POLY_cddlib.h"
#include <unordered_set>


// The list of region is described by a vector of + and -.
// If it is + then the encoding is by a 1. Otherwise it is by a 0.
template<typename T>
std::vector<Face> EnumerateHyperplaneRegions(MyMatrix<T> const& ListV)
{
  int n=ListV.cols();
  int n_vect=ListV.rows();
  auto GetSingleEntry=[&]() -> Face {
    auto try_vect=[&](MyVector<T> const& V) -> std::pair<bool,Face> {
      Face f(n_vect);
      for (int i_vect=0; i_vect<n_vect; i_vect++) {
        T eScal=0;
        for (int i=0; i<n; i++)
          eScal += V(i) * ListV(i_vect,i);
        if (eScal == 0)
          return {false,f};
        if (eScal > 0)
          f[i_vect] = 1;
      }
      return {true,f};
    };
    while(true) {
      MyVector<T> eV(n);
      for (int i=0; i<n; i++) {
        int eVal = rand() % 10;
        eV(i) = eVal;
      }
      std::pair<bool,Face> ePair = try_vect(eV);
      if (ePair.first)
        return ePair.second;
    }
  };
  std::unordered_set<Face> ListDone;
  std::unordered_set<Face> ListUndone;
  ListUndone.insert(GetSingleEntry());
  auto fInsert=[&](Face const& f) -> void {
    if (ListDone.count(f) == 0)
      return;
    ListUndone.insert(f);
  };
  auto ProcessAdjacent=[&](Face const& eF) -> void {
    MyMatrix<T> ListInequalities(n_vect,n);
    for (int i_vect=0; i_vect<n_vect; i_vect++) {
      int eSign=-1;
      if (eF[i_vect])
        eSign=1;
      for (int i=0; i<n; i++)
        ListInequalities(i_vect,i) = eSign * ListV(i_vect,i);
    }
    std::vector<int> ListIrred = cdd::RedundancyReductionClarkson(ListInequalities);
    for (auto & idx : ListIrred) {
      Face newF = eF;
      newF[idx] = 1 - eF[idx];
      fInsert(newF);
    }
  };
  while(true) {
    if (ListUndone.size() == 0)
      break;
    Face f = *(ListUndone.begin());
    ProcessAdjacent(f);
    ListUndone.erase(f);
    ListDone.insert(f);
  }
  std::vector<Face> ListFace;
  for (auto & eF : ListDone)
    ListFace.push_back(eF);
  return ListFace;
}


#endif
