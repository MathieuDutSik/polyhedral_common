// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLYNORM_POLYNORM_PACKING_H_
#define SRC_POLYNORM_POLYNORM_PACKING_H_

// clang-format off
#include "PolyNorm_Basic.h"
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_POLYNORM_PACKING
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_POLYNORM_PACKING
#endif

/*
  The packing problem.

  The translates alpha P + v and alpha P + w have disjoint interiors if and
  only if w - v is not in int(alpha P) - int(alpha P) = int(alpha (P - P)),
  that is ||w - v||_D >= alpha for the difference body D = P - P. So the
  largest alpha for which alpha P + Z^n is a packing is
     alpha(P) = min_{z in Z^n, z != 0} ||z||_D,
  the first minimum of the (centrally symmetric) difference body. The
  lattice vectors attaining it are the contact vectors: the translates
  alpha P and alpha P + z touch exactly for those z.

  The computation is a bounded enumeration: any nonzero lattice vector
  gives an upper bound r0 on alpha, and the minimum is over the lattice
  points of r0 D, enumerated with the LLL-reduced iteration of
  GetListIntegralPoint.
 */

template <typename T, typename Tint> struct PolyNormPacking {
  T alpha;
  // The nonzero lattice vectors z with ||z||_D = alpha. They come in
  // opposite pairs since D is centrally symmetric.
  std::vector<MyVector<Tint>> ListContact;
};

template <typename T, typename Tint>
PolyNormPacking<T, Tint> ComputePolyNormPacking(MyMatrix<T> const &LmatDiff,
                                                std::ostream &os) {
  int n = LmatDiff.cols();
  // The upper bound from the basis vectors. Both signs give the same value
  // since D is symmetric, so only e_i is evaluated.
  std::optional<T> r0;
  for (int i = 0; i < n; i++) {
    MyVector<T> e = ZeroVector<T>(n);
    e(i) = 1;
    T val = PolyNorm_Gauge(LmatDiff, e);
    if (!r0 || val < *r0) {
      r0 = val;
    }
  }
#ifdef DEBUG_POLYNORM_PACKING
  os << "POLYNORM_PACKING: upper bound r0=" << *r0 << "\n";
#endif
  std::vector<MyVector<Tint>> ListPt =
      PolyNorm_LatticePoints<T, Tint>(LmatDiff, *r0, os);
#ifdef DEBUG_POLYNORM_PACKING
  os << "POLYNORM_PACKING: |Z^n cap r0 D|=" << ListPt.size() << "\n";
#endif
  std::optional<T> alpha;
  std::vector<MyVector<Tint>> ListContact;
  for (auto &z : ListPt) {
    if (!IsZeroVector(z)) {
      MyVector<T> z_T = UniversalVectorConversion<T, Tint>(z);
      T val = PolyNorm_Gauge(LmatDiff, z_T);
      if (!alpha || val < *alpha) {
        alpha = val;
        ListContact.clear();
      }
      if (val == *alpha) {
        ListContact.push_back(z);
      }
    }
  }
#ifdef SANITY_CHECK_POLYNORM_PACKING
  if (!alpha) {
    // e_i is in r0 D so the enumeration cannot be empty.
    std::cerr << "POLYNORM_PACKING: no nonzero lattice point found in r0 D\n";
    throw TerminalException{1};
  }
  if (*alpha > *r0) {
    std::cerr << "POLYNORM_PACKING: the minimum exceeds the upper bound\n";
    throw TerminalException{1};
  }
#endif
#ifdef DEBUG_POLYNORM_PACKING
  os << "POLYNORM_PACKING: alpha=" << *alpha
     << " |ListContact|=" << ListContact.size() << "\n";
#endif
  return {*alpha, std::move(ListContact)};
}

// clang-format off
#endif  // SRC_POLYNORM_POLYNORM_PACKING_H_
// clang-format on
