// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_FREEVECTORS_H_
#define SRC_DELAUNAY_FREEVECTORS_H_

// clang-format off
#include "LatticeDelaunay.h"
#include "LatticeStabEquiCan.h"
#include "POLY_cdd_graph.h"
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

// Port of the free-vector computation from
// GAPpackages/MyPolyhedral/lib/LatticeDelaunays.g (FreenessBeltListing /
// FuncDirectEnumerationFreeVector / FuncGetAllFreeVectors).
//
// The pipeline, starting from the Delaunay tesselation of the lattice:
//  1. The Voronoi relevant vectors = the edge vectors of the Delaunay
//     tesselation (differences of Delaunay-adjacent vertices), reduced to one
//     per +/- pair and closed under the point group.
//  2. The belts: triples {i,j,k} of relevant vectors with r_k = +/-(r_i +/- r_j)
//     (the three directions of a hexagonal 6-belt).
//  3. The free vectors: by FuncDirectEnumerationFreeVector, the configurations
//     where the span of selected relevant vectors (one per belt) leaves a
//     non-trivial orthogonal complement; a free vector is that complement
//     (1-dimensional case), mapped by the inverse Gram matrix.

#ifdef DEBUG
#define DEBUG_FREE_VECTORS
#endif

template <typename Tint>
MyVector<Tint> free_antipodal_canon(MyVector<Tint> const &v) {
  int len = v.size();
  for (int i = 0; i < len; i++) {
    if (v(i) != 0) {
      if (v(i) < 0) {
        return MyVector<Tint>(-v);
      }
      return v;
    }
  }
  return v;
}

template <typename Tint>
MyVector<Tint> free_affine_part(MyMatrix<Tint> const &EXT, int iRow) {
  int n = EXT.cols() - 1;
  MyVector<Tint> v(n);
  for (int i = 0; i < n; i++) {
    v(i) = EXT(iRow, i + 1);
  }
  return v;
}

template <typename Tint> struct FreeVectorOrbit {
  int subspace_dim;
  MyVector<Tint> free_vector; // valid when subspace_dim == 1
  int n_matched;              // number of relevant vectors spanned
  std::string orbit_size;     // group orbit size (decimal, can be large)
};

template <typename Tint> struct FreeVectorsResult {
  int n_relevant;                          // number of relevant vectors (mod +/-)
  std::vector<FreeVectorOrbit<Tint>> orbits;
};

template <typename T, typename Tint, typename Tgroup>
FreeVectorsResult<Tint>
compute_free_vectors(MyMatrix<T> const &GramMat,
                     DelaunayTesselation<T, Tgroup> const &DT,
                     std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int n = GramMat.rows();
  //
  // 1. The relevant vectors: edge vectors of the Delaunay tesselation, reduced
  //    to one per +/- pair and closed under the point group of the lattice.
  //
  std::vector<MyMatrix<Tint>> autom =
      ArithmeticAutomorphismGroup<T, Tint, Tgroup>(GramMat, os);
  std::unordered_set<MyVector<Tint>> set_vect;
  std::vector<MyVector<Tint>> ListVectRed;
  auto f_insert = [&](MyVector<Tint> const &v) -> void {
    MyVector<Tint> cv = free_antipodal_canon(v);
    if (set_vect.insert(cv).second) {
      ListVectRed.push_back(cv);
    }
  };
  for (auto &eDel : DT.l_dels) {
    MyMatrix<T> const &EXT_T = eDel.EXT;
    MyMatrix<Tint> EXT = UniversalMatrixConversion<Tint, T>(EXT_T);
    cdd::DDA<T> dda = cdd::DualDescriptionAdjacencies(EXT_T, os);
    int nbVert = EXT.rows();
    for (int i = 0; i < nbVert; i++) {
      for (int j = i + 1; j < nbVert; j++) {
        if (dda.SkelGraph.IsAdjacent(i, j)) {
          MyVector<Tint> r =
              free_affine_part(EXT, i) - free_affine_part(EXT, j);
          f_insert(r);
        }
      }
    }
  }
  // Closure under the point group (action v -> U^T v preserves the Gram norm
  // since U Gram U^T = Gram for the automorphisms).
  size_t pos = 0;
  while (pos < ListVectRed.size()) {
    MyVector<Tint> v = ListVectRed[pos];
    pos++;
    for (auto &U : autom) {
      MyVector<Tint> w = U.transpose() * v;
      f_insert(w);
    }
  }
  int nbVect = ListVectRed.size();
#ifdef DEBUG_FREE_VECTORS
  os << "FREEVECT: number of relevant vectors (mod +/-)=" << nbVect << "\n";
#endif
  std::unordered_map<MyVector<Tint>, int> map_vect;
  for (int i = 0; i < nbVect; i++) {
    map_vect[ListVectRed[i]] = i;
  }
  auto index_of = [&](MyVector<Tint> const &v) -> int {
    MyVector<Tint> cv = free_antipodal_canon(v);
    auto iter = map_vect.find(cv);
    if (iter == map_vect.end()) {
      return -1;
    }
    return iter->second;
  };
  //
  // The permutation group on the relevant vectors.
  //
  std::vector<Telt> ListPermGens;
  for (auto &U : autom) {
    std::vector<Tidx> eList(nbVect);
    for (int k = 0; k < nbVect; k++) {
      MyVector<Tint> w = U.transpose() * ListVectRed[k];
      int img = index_of(w);
      eList[k] = static_cast<Tidx>(img);
    }
    ListPermGens.push_back(Telt(eList));
  }
  Tgroup GRPperm(ListPermGens, nbVect);
  //
  // 2. The belts: triples {i,j,k} with r_k = +/-(r_i +/- r_j).
  //
  std::vector<std::vector<int>> ListBelt;
  std::set<std::vector<int>> SetBelt;
  auto try_belt = [&](int i, int j, MyVector<Tint> const &w) -> void {
    int k = index_of(w);
    if (k != -1 && k != i && k != j) {
      std::vector<int> triple{i, j, k};
      std::sort(triple.begin(), triple.end());
      if (SetBelt.insert(triple).second) {
        ListBelt.push_back(triple);
      }
    }
  };
  for (int i = 0; i < nbVect; i++) {
    for (int j = i + 1; j < nbVect; j++) {
      try_belt(i, j, ListVectRed[i] - ListVectRed[j]);
      try_belt(i, j, ListVectRed[i] + ListVectRed[j]);
    }
  }
  int nbBelt = ListBelt.size();
#ifdef DEBUG_FREE_VECTORS
  os << "FREEVECT: number of belts=" << nbBelt << "\n";
#endif
  //
  // 3. The free-vector enumeration (FuncDirectEnumerationFreeVector).
  //
  struct Sol {
    Face matchedBelt;
    Face matchedVectors;
    MyMatrix<Tint> basis;
  };
  Sol init{Face(nbBelt), Face(nbVect), IdentityMat<Tint>(n)};
  std::vector<Sol> ListSolution{init};
  std::vector<Sol> ListTotal;
  while (!ListSolution.empty()) {
    std::vector<Sol> NewList;
    auto FuncInsert = [&](Sol const &s) -> void {
      for (auto &e : NewList) {
        if (GRPperm.RepresentativeAction_OnSets(e.matchedVectors,
                                                s.matchedVectors)) {
          return;
        }
      }
      NewList.push_back(s);
    };
    for (auto &eSol : ListSolution) {
      int iBeltMin = -1;
      for (int b = 0; b < nbBelt; b++) {
        if (eSol.matchedBelt[b] == 0) {
          iBeltMin = b;
          break;
        }
      }
      if (iBeltMin == -1) {
        ListTotal.push_back(eSol);
        continue;
      }
      for (int idx : ListBelt[iBeltMin]) {
        // VectorSet = relevant vectors in matchedVectors, plus r_idx.
        std::vector<int> sel;
        for (auto b = eSol.matchedVectors.find_first(); b != Face::npos;
             b = eSol.matchedVectors.find_next(b)) {
          sel.push_back(static_cast<int>(b));
        }
        sel.push_back(idx);
        MyMatrix<Tint> VectorSet(sel.size(), n);
        for (size_t u = 0; u < sel.size(); u++) {
          for (int c = 0; c < n; c++) {
            VectorSet(u, c) = ListVectRed[sel[u]](c);
          }
        }
        MyMatrix<Tint> NSP = NullspaceTrMat(VectorSet);
        int dimNSP = NSP.rows();
        if (dimNSP > 0) {
          Face newMatched(nbVect);
          for (int iVect = 0; iVect < nbVect; iVect++) {
            bool ortho = true;
            for (int d = 0; d < dimNSP; d++) {
              Tint scal(0);
              for (int c = 0; c < n; c++) {
                scal += NSP(d, c) * ListVectRed[iVect](c);
              }
              if (scal != 0) {
                ortho = false;
                break;
              }
            }
            if (ortho) {
              newMatched[iVect] = 1;
            }
          }
          Face newBelt(nbBelt);
          for (int iBelt = 0; iBelt < nbBelt; iBelt++) {
            for (int idx2 : ListBelt[iBelt]) {
              if (newMatched[idx2] == 1) {
                newBelt[iBelt] = 1;
                break;
              }
            }
          }
          FuncInsert(Sol{newBelt, newMatched, NSP});
        }
      }
    }
    ListSolution = NewList;
  }
  //
  // Deduplicate the final solutions under the group and build the result.
  //
  MyMatrix<T> GramInv = Inverse(GramMat);
  FreeVectorsResult<Tint> result;
  result.n_relevant = nbVect;
  std::vector<Face> seen;
  for (auto &eSol : ListTotal) {
    bool is_new = true;
    for (auto &f : seen) {
      if (GRPperm.RepresentativeAction_OnSets(f, eSol.matchedVectors)) {
        is_new = false;
        break;
      }
    }
    if (!is_new) {
      continue;
    }
    seen.push_back(eSol.matchedVectors);
    FreeVectorOrbit<Tint> orb;
    orb.subspace_dim = eSol.basis.rows();
    orb.n_matched = eSol.matchedVectors.count();
    std::stringstream ss_orb;
    ss_orb << GRPperm.OrbitSize_OnSets(eSol.matchedVectors);
    orb.orbit_size = ss_orb.str();
    if (orb.subspace_dim == 1) {
      MyVector<T> u(n);
      for (int c = 0; c < n; c++) {
        u(c) = UniversalScalarConversion<T, Tint>(eSol.basis(0, c));
      }
      MyVector<T> w = GramInv * u;
      MyVector<T> wr = RemoveFractionVector(w);
      orb.free_vector = UniversalVectorConversion<Tint, T>(wr);
    } else {
      orb.free_vector = ZeroVector<Tint>(n);
    }
    result.orbits.push_back(orb);
  }
  return result;
}

template <typename Tint>
void WriteFreeVectorsGAP(std::ostream &os_out,
                         FreeVectorsResult<Tint> const &res) {
  os_out << "return rec(nbRelevantVector:=" << res.n_relevant << ",\n";
  os_out << "nbFreeVectorOrbit:=" << res.orbits.size() << ",\n";
  os_out << "ListFreeVectorOrbit:=[";
  bool IsFirst = true;
  for (auto &orb : res.orbits) {
    if (!IsFirst) {
      os_out << ",\n";
    }
    IsFirst = false;
    os_out << "rec(SubspaceDim:=" << orb.subspace_dim
           << ", OrbitSize:=" << orb.orbit_size
           << ", nMatched:=" << orb.n_matched
           << ", eVect:=" << StringVectorGAP(orb.free_vector) << ")";
  }
  os_out << "]);\n";
}

// clang-format off
#endif  // SRC_DELAUNAY_FREEVECTORS_H_
// clang-format on
