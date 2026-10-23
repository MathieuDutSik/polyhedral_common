// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_IGUSA_IGUSA_INITIAL_FACET_H_
#define SRC_IGUSA_IGUSA_INITIAL_FACET_H_

// clang-format off
#include "igusa_face_vertices.h"
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_IGUSA_INITIAL_FACET
#endif

/*
  An initial facet of a face Phi of P = conv(I_n) of rank r < n.

  Phi is the face of the valid pair (F, b), the lift through W of the
  bounded face Phi0 of P_r. For w primitive modulo W (W + Z w saturated),
  X[w] >= 1 is tight on Phi and Psi_w = Phi cap {X[w] = 1} is the face of
  the valid pair (F + w w^T, b + 1), of rank r + 1. All such w are
  equivalent under Stab(Phi). In a basis (B, w, ...) the points of Psi_w
  over Y in Inc(Phi0) have the blocks Y, b, 1 with 2b = c integral, and are
  positive definite iff Y^{-1}[c] < 4. Hence
      codim of Psi_w in Phi = 1 + r - rank C(Phi0),
      C(Phi0) = { c in Z^r : Y^{-1}[c] < 4 for some Y in Inc(Phi0) },
  and Psi_w is a facet of Phi iff C(Phi0) spans R^r.

  The vectors c are first taken over the vertices of Phi0, which gives a
  sufficient condition; if they do not span, all the integral points of
  Phi0 are computed, which gives the exact codimension.
 */

template <typename T, typename Tint> struct IgusaInitialFacet {
  // The rank r of the face and the vector w, primitive modulo W
  int rank;
  MyVector<Tint> w;
  // The rank of C(Phi0), and r independent vectors of it when it spans
  int rank_C;
  std::vector<MyVector<Tint>> l_certificate;
  // Whether all the integral points of Phi0 were needed
  bool used_all_points;
  // The codimension of Psi_w in Phi, 1 if it is a facet
  int codim;
  // The valid pair of Psi_w
  MyMatrix<T> Fnew;
  T bnew;
  // The vertices of Phi0, as found by igusa_face_vertices
  IgusaFaceVertices<T> face;
};

/*
  The vectors c of Z^r with Y^{-1}[c] < 4 for the forms Y of l_Y, added to
  l_c while they increase the rank. Returns the rank.
 */
template <typename T, typename Tint>
int igusa_insert_c_vectors(std::vector<MyMatrix<T>> const &l_Y, int r,
                           std::vector<MyVector<Tint>> &l_c,
                           std::ostream &os) {
  auto get_rank = [&]() -> int {
    if (l_c.empty()) {
      return 0;
    }
    MyMatrix<T> M(l_c.size(), r);
    for (size_t i = 0; i < l_c.size(); i++) {
      for (int j = 0; j < r; j++) {
        M(i, j) = UniversalScalarConversion<T, Tint>(l_c[i](j));
      }
    }
    return RankMat(M);
  };
  int rank = get_rank();
  for (auto &Y : l_Y) {
    if (rank == r) {
      break;
    }
    MyMatrix<T> Yinv = Inverse(Y);
    for (auto &c : computeLevel_GramMat<T, Tint>(Yinv, T(4), os)) {
      if (IsZeroVector(c) || EvaluationQuadForm<T, Tint>(Yinv, c) >= 4) {
        continue;
      }
      l_c.push_back(c);
      int rank_new = get_rank();
      if (rank_new == rank) {
        l_c.pop_back();
        continue;
      }
      rank = rank_new;
      if (rank == r) {
        break;
      }
    }
  }
  return rank;
}

template <typename T, typename Tint, typename Tgroup>
IgusaInitialFacet<T, Tint>
igusa_initial_facet(IgusaParameters<T> const &params, MyMatrix<T> const &F,
                    T const &b, std::vector<MyMatrix<T>> const &ListRepr,
                    std::ostream &os) {
  int n = F.rows();
  IgusaInitialFacet<T, Tint> ret;
  ret.face = igusa_face_vertices<T, Tint, Tgroup>(params, F, b, ListRepr, os);
  int r = ret.face.rank;
  ret.rank = r;
  if (r == n) {
    std::cerr << "IGUSA_INITIAL_FACET: the face is bounded, its facets are "
                 "the facets of the polytope of its vertices\n";
    throw TerminalException{1};
  }
  if (ret.face.l_vertex.empty()) {
    std::cerr << "IGUSA_INITIAL_FACET: no vertex found, (F, b) is not a "
                 "valid pair or the vertex representatives are incomplete\n";
    throw TerminalException{1};
  }
  // w: the first vector of a completion of the basis B of W
  MyMatrix<T> Compl = SubspaceCompletionInt(ret.face.B, n);
  MyVector<T> w_T = GetMatrixRow(Compl, 0);
  ret.w = UniversalVectorConversion<Tint, T>(w_T);
  ret.used_all_points = false;
  ret.rank_C = igusa_insert_c_vectors<T, Tint>(ret.face.l_vertex, r,
                                               ret.l_certificate, os);
  if (ret.rank_C < r) {
    // The exact rank needs all the integral points of Phi0
    LinSpaceMatrix<T> LinSpa_r = ComputeCanonicalSpace<T>(r);
    IgusaSpace<T, Tint> space_r =
        build_igusa_space<T, Tint>(LinSpa_r, params, os);
    igusa_insert_initial_cuts(space_r);
    IgusaFacetIncidence<T> inc =
        igusa_facet_incidence(space_r, ret.face.F0, b, os, false);
    ret.used_all_points = true;
    ret.rank_C = igusa_insert_c_vectors<T, Tint>(inc.l_incident, r,
                                                 ret.l_certificate, os);
  }
  ret.codim = 1 + r - ret.rank_C;
  MyMatrix<T> wwT = w_T * w_T.transpose();
  ret.Fnew = F + wwT;
  ret.bnew = b + 1;
#ifdef DEBUG_IGUSA_INITIAL_FACET
  os << "IGUSA_INITIAL_FACET: r=" << r << " w=" << StringVector(ret.w)
     << " rank C=" << ret.rank_C << " codim=" << ret.codim << "\n";
#endif
  return ret;
}

inline FullNamelist NAMELIST_GetStandard_IGUSA_INITIAL_FACET() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  {
    std::map<std::string, std::string> ListStringValues_doc;
    ListStringValues_doc["arithmetic"] =
        doc_rational_integer_arithmetics_by_library();
    ListStringValues_doc["IlpMethod"] = doc_IlpMethod();
    ListStringValues_doc["OutFile"] = doc_StderrStdoutFile("stderr");
    ListStringValues_doc["FileFacet"] = "The file containing the matrix F\n\
(n x n, positive semidefinite of rank r < n) of the face tr(F X) >= rhs";
    ListStringValues_doc["FacetRhs"] = "The right hand side rhs of the face,\n\
a rational number";
    ListStringValues_doc["FileListVertex"] = "The file containing the\n\
representatives of the orbits of vertices of P_r, in ListMatrix format";
    ListBlock["DATA"].setListStringValues_doc(ListStringValues_doc);
  }
  return FullNamelist(ListBlock);
}

// clang-format off
#endif  // SRC_IGUSA_IGUSA_INITIAL_FACET_H_
// clang-format on
