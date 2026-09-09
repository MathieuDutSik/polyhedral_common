// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEPLESKENSOUVIGNIER_H_
#define SRC_LATT_LATTICEPLESKENSOUVIGNIER_H_

// clang-format off
#include "PleskenSouvignier.h"
#include "Shvec_exact.h"
#include "ClassicLLL.h"
#include <optional>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_PLESKEN_SOUVIGNIER
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_PLESKEN_SOUVIGNIER
#endif

#ifdef TIMINGS
#define TIMINGS_PLESKEN_SOUVIGNIER
#endif

/*
  The lattice-level entry points of the Plesken-Souvignier engine of
  src_group/PleskenSouvignier.h: they LLL-reduce the Gram matrix, build
  the vector set of the paper -- the vectors of norm at most the largest
  diagonal entry of the reduced Gram matrix, which contains the standard
  basis and is invariant -- run the engine in the reduced basis, and
  transport the answers back to the original basis. The short vector
  enumeration and the LLL are what the engine header cannot depend on,
  which is why these wrappers live in src_latt.
 */

// The vector set of the paper for a Gram matrix, one vector per
// antipodal pair: the vectors of norm at most max_i GramMat(i, i).
template <typename T, typename Tint>
MyMatrix<Tint> PleskenSouvignierVectorFamily(MyMatrix<T> const &GramMat,
                                             T const &bound,
                                             std::ostream &os) {
  int n = GramMat.rows();
  std::vector<MyVector<Tint>> ListVect =
      computeLevel_GramMat<T, Tint>(GramMat, bound, os);
  return MatrixFromVectorFamilyDim(n, ListVect);
}

template <typename T> T PleskenSouvignierBound(MyMatrix<T> const &GramMat) {
  int n = GramMat.rows();
  T bound = GramMat(0, 0);
  for (int i = 1; i < n; i++) {
    if (GramMat(i, i) > bound) {
      bound = GramMat(i, i);
    }
  }
  return bound;
}

template <typename T, typename Tint>
std::vector<MyMatrix<Tint>>
ps_convert_integral_listmat(std::vector<MyMatrix<T>> const &ListMat) {
  std::vector<MyMatrix<Tint>> ListMatRet;
  ListMatRet.reserve(ListMat.size());
  for (auto &eMat : ListMat) {
    std::optional<MyMatrix<Tint>> opt =
        UniversalMatrixConversionCheck<Tint, T>(eMat);
    if (!opt) {
      std::cerr << "PS: the matrices of the configuration have to be "
                << "integral\n";
      throw TerminalException{1};
    }
    ListMatRet.push_back(*opt);
  }
  return ListMatRet;
}

/*
  The automorphism group of the configuration: generators over the ring
  in the ORIGINAL basis, satisfying g * M * g^T = M for every matrix of
  ListMat, and the orbit lengths along the stabilizer chain whose
  product is the group order. ListMat[0] has to be symmetric positive
  definite with integral entries, the others symmetric integral.
 */
template <typename T, typename Tint>
PleskenSouvignierAutomResult<Tint>
PleskenSouvignierLatticeAutomorphism(std::vector<MyMatrix<T>> const &ListMat,
                                     std::ostream &os) {
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(ListMat[0], os);
  MyMatrix<Tint> const &Pmat = rec.Pmat;
  MyMatrix<Tint> PmatInv = Inverse(Pmat);
  std::vector<MyMatrix<Tint>> ListMatInt =
      ps_convert_integral_listmat<T, Tint>(ListMat);
  std::vector<MyMatrix<Tint>> ListMatRed;
  ListMatRed.reserve(ListMatInt.size());
  for (auto &eMat : ListMatInt) {
    ListMatRed.push_back(Pmat * eMat * Pmat.transpose());
  }
  T bound = PleskenSouvignierBound(rec.GramMatRed);
  MyMatrix<Tint> SHVhalf =
      PleskenSouvignierVectorFamily<T, Tint>(rec.GramMatRed, bound, os);
#ifdef DEBUG_PLESKEN_SOUVIGNIER
  os << "PS: automorphism, bound " << bound << ", " << SHVhalf.rows()
     << " antipodal pairs\n";
#endif
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  os << "|PS: LatticeAutomorphism, preparation|=" << time << "\n";
#endif
  PleskenSouvignierAutomResult<Tint> result =
      PleskenSouvignierAutomorphism<Tint>(ListMatRed, SHVhalf, os);
  // g_red preserves P M P^T, so P^{-1} g_red P preserves M, and it is
  // integral, P being unimodular.
  for (auto &eGen : result.ListGen) {
    eGen = PmatInv * eGen * Pmat;
  }
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  for (auto &eGen : result.ListGen) {
    for (auto &eMat : ListMatInt) {
      if (eGen * eMat * eGen.transpose() != eMat) {
        std::cerr << "PS: a transported generator does not preserve the "
                  << "configuration\n";
        throw TerminalException{1};
      }
    }
  }
#endif
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  os << "|PS: LatticeAutomorphism|=" << time << "\n";
#endif
  return result;
}

/*
  The isometry test: P with P * ListMat1[i] * P^T = ListMat2[i] for
  every i, in the original bases, or nothing. The vector sets of the two
  sides are cut at the SAME bound, the one of the first lattice, which
  is what makes them correspond under any isometry. ListGenAut2 is any
  set of automorphisms of the second configuration in its original
  basis, typically the generators the caller already knows for a lattice
  it tests many candidates against; the orbits they generate prune the
  search, and an empty list is always allowed.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<Tint>> PleskenSouvignierLatticeIsometry(
    std::vector<MyMatrix<T>> const &ListMat1,
    std::vector<MyMatrix<T>> const &ListMat2,
    std::vector<MyMatrix<Tint>> const &ListGenAut2, std::ostream &os) {
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  LLLreduction<T, Tint> rec1 = LLLreducedBasis<T, Tint>(ListMat1[0], os);
  LLLreduction<T, Tint> rec2 = LLLreducedBasis<T, Tint>(ListMat2[0], os);
  MyMatrix<Tint> const &Pmat1 = rec1.Pmat;
  MyMatrix<Tint> const &Pmat2 = rec2.Pmat;
  std::vector<MyMatrix<Tint>> ListMatInt1 =
      ps_convert_integral_listmat<T, Tint>(ListMat1);
  std::vector<MyMatrix<Tint>> ListMatInt2 =
      ps_convert_integral_listmat<T, Tint>(ListMat2);
  std::vector<MyMatrix<Tint>> ListMatRed1;
  std::vector<MyMatrix<Tint>> ListMatRed2;
  for (auto &eMat : ListMatInt1) {
    ListMatRed1.push_back(Pmat1 * eMat * Pmat1.transpose());
  }
  for (auto &eMat : ListMatInt2) {
    ListMatRed2.push_back(Pmat2 * eMat * Pmat2.transpose());
  }
  T bound = PleskenSouvignierBound(rec1.GramMatRed);
  MyMatrix<Tint> SHVhalf1 =
      PleskenSouvignierVectorFamily<T, Tint>(rec1.GramMatRed, bound, os);
  MyMatrix<Tint> SHVhalf2 =
      PleskenSouvignierVectorFamily<T, Tint>(rec2.GramMatRed, bound, os);
#ifdef DEBUG_PLESKEN_SOUVIGNIER
  os << "PS: isometry, bound " << bound << ", pairs " << SHVhalf1.rows()
     << " / " << SHVhalf2.rows() << "\n";
#endif
  MyMatrix<Tint> Pmat2Inv = Inverse(Pmat2);
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  // A wrong generator would surface deep in the search as a family
  // invariance error, so catch the misuse at the boundary.
  for (auto &eGen : ListGenAut2) {
    for (auto &eMat : ListMatInt2) {
      if (eGen * eMat * eGen.transpose() != eMat) {
        std::cerr << "PS: ListGenAut2 has to consist of automorphisms of "
                  << "the SECOND configuration\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  std::vector<MyMatrix<Tint>> ListGenAut2Red;
  ListGenAut2Red.reserve(ListGenAut2.size());
  for (auto &eGen : ListGenAut2) {
    ListGenAut2Red.push_back(Pmat2 * eGen * Pmat2Inv);
  }
  std::optional<MyMatrix<Tint>> opt = PleskenSouvignierIsometry<Tint>(
      ListMatRed1, SHVhalf1, ListMatRed2, SHVhalf2, ListGenAut2Red, os);
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  os << "|PS: LatticeIsometry found=" << opt.has_value() << "|=" << time
     << "\n";
#endif
  if (!opt) {
    return {};
  }
  // P_red * (P1 M1 P1^T) * P_red^T = P2 M2 P2^T, so the original-basis
  // equivalence is P2^{-1} P_red P1.
  MyMatrix<Tint> P = Pmat2Inv * (*opt) * Pmat1;
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  for (size_t iMat = 0; iMat < ListMatInt1.size(); iMat++) {
    if (P * ListMatInt1[iMat] * P.transpose() != ListMatInt2[iMat]) {
      std::cerr << "PS: the transported isometry does not map ListMat1 to "
                << "ListMat2\n";
      throw TerminalException{1};
    }
  }
#endif
  return P;
}

// clang-format off
#endif  // SRC_LATT_LATTICEPLESKENSOUVIGNIER_H_
// clang-format on
