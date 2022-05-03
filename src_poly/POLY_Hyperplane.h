#ifndef SRC_POLY_POLY_HYPERPLANE_H_
#define SRC_POLY_POLY_HYPERPLANE_H_

#include "POLY_c_cddlib.h"
#include "POLY_cddlib.h"
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

// The list of region is described by a vector of + and -.
// If it is + then the encoding is by a 1. Otherwise it is by a 0.
template <typename T>
vectface EnumerateHyperplaneRegions(MyMatrix<T> const &ListV) {
  int n = ListV.cols();
  int n_vect = ListV.rows();
#ifdef DEBUG_HYPERPLANE
  auto StringFace = [&](Face const &f) -> std::string {
    std::string strO = "[";
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      if (i_vect > 0)
        strO += ",";
      if (f[i_vect])
        strO += "+";
      else
        strO += "-";
    }
    strO += "]";
    return strO;
  };
#endif
  auto GetSingleEntry = [&]() -> Face {
    auto try_vect = [&](MyVector<T> const &V) -> std::pair<bool, Face> {
      Face f(n_vect);
      for (int i_vect = 0; i_vect < n_vect; i_vect++) {
        T eScal = 0;
        for (int i = 0; i < n; i++)
          eScal += V(i) * ListV(i_vect, i);
        if (eScal == 0)
          return {false, f};
        if (eScal > 0)
          f[i_vect] = 1;
      }
      return {true, f};
    };
    while (true) {
      MyVector<T> eV(n);
      for (int i = 0; i < n; i++) {
        int eVal = rand() % 10;
        eV(i) = eVal;
      }
      std::pair<bool, Face> ePair = try_vect(eV);
      if (ePair.first)
        return ePair.second;
    }
  };
  std::unordered_set<Face> ListDone;
  std::unordered_set<Face> ListUndone;
  ListUndone.insert(GetSingleEntry());
#ifdef USE_CDDLIB
  auto fInsert = [&](Face const &f) -> void {
    if (ListDone.count(f) != 0)
      return;
    ListUndone.insert(f);
  };
#endif
  auto ProcessAdjacent = [&](Face const &eF) -> void {
#ifdef DEBUG_HYPERPLANE
    std::cerr << "ProcessAdjacent eF=" << StringFace(eF) << "\n";
#endif
    MyMatrix<T> ListInequalities(n_vect, n + 1);
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      int eSign = -1;
      if (eF[i_vect])
        eSign = 1;
      ListInequalities(i_vect, 0) = 0;
      for (int i = 0; i < n; i++)
        ListInequalities(i_vect, i + 1) = eSign * ListV(i_vect, i);
    }
#ifdef DEBUG_HYPERPLANE
    std::cerr << "ListInequalities=\n";
    WriteMatrix(std::cerr, ListInequalities);
#endif
#ifdef USE_CDDLIB
    std::vector<int> ListIrred =
        cbased_cdd::RedundancyReductionClarkson(ListInequalities);
    //    std::vector<int> ListIrred =
    //    cdd::RedundancyReductionClarkson(ListInequalities);
#ifdef DEBUG_HYPERPLANE
    std::cerr << "len(ListIrred)=" << ListIrred.size() << "\n";
#endif
#ifdef CHECK_HYPERPLANE
    if (static_cast<int>(ListIrred.size()) < n) {
      std::cerr << "RankMat(...)=" << RankMat(ListInequalities) << "\n";
      std::cerr << "ListInequalities=\n";
      WriteMatrix(std::cerr, ListInequalities);
      std::cerr << "|ListIrred|=" << ListIrred.size() << "\n";
      std::cerr << "ListIrred is too small. It is a bug\n";
      throw TerminalException{1};
    }
#endif
    for (auto &idx : ListIrred) {
      Face newF = eF;
      newF[idx] = 1 - eF[idx];
#ifdef DEBUG_HYPERPLANE
      std::cerr << "idx=" << idx << " newF=" << StringFace(newF) << "\n";
#endif
      fInsert(newF);
    }
#else
    std::cerr << "We need to compile with USE_CDDLIB\n";
    throw TerminalException{1};
#endif
  };
  while (true) {
    if (ListUndone.size() == 0)
      break;
    Face f = *(ListUndone.begin());
    ProcessAdjacent(f);
    ListUndone.erase(f);
    ListDone.insert(f);
#ifdef DEBUG_HYPERPLANE
    std::cerr << " len(ListDone)=" << ListDone.size()
              << " len(ListUndone)=" << ListUndone.size() << "\n";
#endif
  }
  vectface ListFace;
  for (auto &eF : ListDone)
    ListFace.push_back(eF);
  return ListFace;
}

// clang-format off
#endif  // SRC_POLY_POLY_HYPERPLANE_H_
// clang-format on
