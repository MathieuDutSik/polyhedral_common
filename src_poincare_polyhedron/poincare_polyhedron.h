// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POINCARE_POLYHEDRON_TH_POINCARE_POLYHEDRON_H_
#define SRC_POINCARE_POLYHEDRON_TH_POINCARE_POLYHEDRON_H_

#include "Namelist.h"
#include "POLY_DirectDualDesc.h"

// The initial data for the Poincare Polyhedron Theorem
// ---a point x
// ---a list of group element which is of course assumed to generate the group
template <typename T> struct DataPoincare {
  MyVector<T> x;
  std::vector<MyMatrix<T>> ListGroupElt;
};

// The data structure for keeîng track of the group elements:
// ---Positive value (so element 0 correspond to 1, X to X+1, ...) are the
// elements themselves.
// ---Negative values correspond to their inverse (that is -1 correspond to
// inverse of generator 0)
// ---0 should never show up in the list.
// DI stands for "Direct or Inverse"
struct TrackGroup {
  std::vector<int> ListDI;
};

TrackGroup ProductTrack(TrackGroup const &tg1, TrackGroup const &tg2) {
  std::vector<int> ListDI = tg1.ListDI;
  ListDI.insert(ListDI.end(), tg2.ListDI.begin(), tg2.ListDI.end());
  return {ListDI};
}

TrackGroup InverseTrack(TrackGroup const &tg) {
  std::vector<int> ListDI_ret;
  size_t len = tg.ListDI.size();
  for (size_t u = 0; u < len; u++) {
    size_t v = len - 1 - u;
    ListDI_ret.push_back(-tg.ListDI[v]);
  }
  return {ListDI_ret};
}

template <typename T> using PairElt<T> = std::pair<TrackGroup, MyMatrix<T>>;

template <typename T>
PairElt<T> ProductPair(PairElt<T> const &p1, PairElt<T> const &p2) {
  return {ProductTrack(p1.first, p2.first), p1.second * p2.second};
}

PairElt<T> InversePair(PairElt<T> const &p) {
  return {InverseTrack(p.first), Inverse(p.second)};
}

template <typename T> struct StepEnum {
  std::vector<PairElt<T>> ListNeighbor;
}

template <typename T>
StepEnum<T> BuildInitialStepEnum(std::vector<MyMatrix<T>> const &ListMat) {
  std::vector<PairElt<T>> ListNeighbor;
  int i = 1;
  for (auto &eMat : ListMat) {
    ListNeighbor.push({i, eMat});
    i++;
  }
  return {ListNeighbor};
}

using TsingAdj = std::pair<size_t, Face>;

template <typename T> struct AdjacencyInfo {
  MyMatrix<T> EXT;
  StepEnum<T> se_red;
  std::vector<std::vector<TsingAdj>> ll_adj:
};

// The domain is defined originally as
// Tr(AX) <= Tr(PAP^T X)
// which we rewrite
// a.x <= phi(a).x        0 <= (phi(a) - a).x
//
// Matrixwise the scalar product a.x is rewritten as
// A X^T
// phi(a) is expressed as AQ
// 0 <= (AQ - A) X^T
// The mapping A ----> AQ maps to the adjacent domain.
// The mapping of the X ----> X c(Q)^T with c(Q) the
// contragredient representation.
template <typename T>
AdjacencyInfo<T> ComputeAdjacencyInfo(MyVector<T> const &x,
                                      StepEnum<T> const &se,
                                      std::string const &eCommand) {
  int n = x.size();
  int n_mat = se.ListNeighbor.size();
  MyMatrix<T> FAC(n_mat, n);
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> eMat = se.ListNeighbor[i_mat];
    MyVector<T> x_img = eMat.transpose() * x;
    MyVector<T> x_diff = x_img - x;
    AssignMatrixRow(FAC, i_mat, x_diff);
  }
  vectface vf = DualDescExternalProgram(FAC, eCommand, std::cerr);
  int n_ext = vf.size();
  MyMatrix<T> EXT(n_ext, n);
  int pos = 0;
  std::vector<Face> v_red(Face(n_ext), n_mat);
  for (auto &eInc : vf) {
    MyVector<T> eEXT = FindFacetInequality(FAC, eInc);
    AssignMatrixRow(EXT, pos, eEXT);
    for (int i_mat = 0; i_mat < n_mat; i_mat++) {
      v_red[i_mat][pos] = eInc[i_mat];
    }
    pos++;
  }
  std::vector<PairElt<T>> ListNeigborRed;
  int n_mat_red = 0;
  std::vector<int> l_i_mat;
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> EXT_red = SelectRow(EXT, v_red[i_mat]);
    int rnk = RankMat(EXT_red);
    if (rnk == n - 1) {
      ListNeighborRed.push_back(se.ListNeighbor[i_mat]);
      l_i_mat.push_back(n_mat_red);
      n_mat_red++;
    }
  }
  std::cerr << "n_mat_red=" << n_mat_red << "\n";
  std::vector<std::vector<TsingAdj>> ll_adj;
  for (int i_mat_red = 0; i_mat_red < n_mat_red; i_mat_red++) {
    int i_mat = l_i_mat[i_mat_red];
    Face const &f1 = v_red[i_mat];
    std::vector<TsingAdj> l_adj;
    for (int j_mat_red = 0; j_mat_red < n_mat_red; j_mat_red++) {
      if (i_mat_red != j_mat_red) {
        int j_mat = l_i_mat[j_mat_red];
        Face const &f2 = v_red[j_mat];
        Face f(n_ext);
        for (int i_ext = 0; i_ext < n_ext; i_ext++)
          if (f1[i_ext] == 1 && f2[i_ext] == 1)
            f[i_ext] = 1;
        MyMatrix<T> EXT_red = SelectRow(EXT, f);
        int rnk = RankMat(EXT_red);
        if (rnk == n - 2) {
          l_adj.push_back({j_mat_red, f});
        }
      }
    }
    ll_adj.push_back(l_adj);
  }
  StepEnum<T> se_red = {ListNeighborRed};
  return {EXT, se_red, ll_adj};
}

template <typename T>
std::option<StepEnum<T>> ComputeMissingNeighbors(AdjacencyInfo<T> const &ai) {}

struct RecOption {
  std::string eCommand;
  int n_iter_max;
};

template <typename T>
AdjacencyInfo<T> IterativePoincareRefinement(DataPoincare<T> const &dp,
                                             RecOption const &ro) {
  MyVector<T> x = dp.x;
  StepEnum<T> se = BuildInitialStepEnum(dp.ListGroupElt);
  int n_iter = 0;
  while (true) {
    std::cerr << "IterativePoincareRefinement n_iter=" << n_iter << "\n";
    AdjacencyInfo<T> ai = ComputeAdjacencyInfo(x, se, ro.eCommand);
    std::option<StepEnum<T>> opt = ComputeMissingNeighbors(ai);
    if (!opt) {
      return ai;
    }
    se = *opt;
    if (ro.n_iter_max > 0) {
      if (n_iter > ro.n_iter_max)
        throw TerminalException{1};
    }
    n_iter++;
  }
}

// clang-format off
#endif  // SRC_POINCARE_POLYHEDRON_TH_POINCARE_POLYHEDRON_H_
// clang-format on
