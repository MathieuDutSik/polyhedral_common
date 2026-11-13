// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEPLESKENSOUVIGNIER_H_
#define SRC_LATT_LATTICEPLESKENSOUVIGNIER_H_

// clang-format off
#include "PleskenSouvignier.h"
#include "Shvec_exact.h"
#include "ClassicLLL.h"
#include "InvariantVectorFamily.h"
#include <map>
#include <optional>
#include <unordered_set>
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
  src_group/PleskenSouvignier.h: they LLL-reduce the Gram matrix, take the
  invariant vector family of the kind asked for (InvariantFamilyHalf of
  InvariantVectorFamily.h), run the engine in the reduced basis, and
  transport the answers back to the original basis. The engine needs a
  family of full rank invariant under the group: it draws the searched
  basis from the family and discards the rational leaves when that basis
  does not span. Which family is best depends on the lattice (see the
  kinds); the choice is the caller's.
 */

/*
  The values (v M v^T)_M of the forms of ListMat on the rows v of SHV, in
  int64 when the bound n^2 max|v|^2 max|M| < 2^62 certifies every product
  and sum, and in Tint otherwise.
 */
template <typename Tint>
std::vector<std::vector<Tint>>
FormValuesOnFamily(std::vector<MyMatrix<Tint>> const &ListMat,
                   MyMatrix<Tint> const &SHV) {
  int m = SHV.rows();
  int n = SHV.cols();
  size_t nbMat = ListMat.size();
  std::vector<std::vector<Tint>> vals(m, std::vector<Tint>(nbMat));
  Tint maxV(0), maxF(0);
  for (int j = 0; j < m; j++) {
    for (int i = 0; i < n; i++) {
      Tint a = T_abs(SHV(j, i));
      if (a > maxV) {
        maxV = a;
      }
    }
  }
  for (auto const &F : ListMat) {
    for (int i = 0; i < n; i++) {
      for (int k = 0; k < n; k++) {
        Tint a = T_abs(F(i, k));
        if (a > maxF) {
          maxF = a;
        }
      }
    }
  }
  Tint bound = Tint(n) * Tint(n) * maxV * maxV * maxF;
  if (bound < UniversalScalarConversion<Tint, int64_t>(int64_t(1) << 62)) {
    MyMatrix<int64_t> SHV64 = UniversalMatrixConversion<int64_t, Tint>(SHV);
    for (size_t iMat = 0; iMat < nbMat; iMat++) {
      MyMatrix<int64_t> W =
          SHV64 * UniversalMatrixConversion<int64_t, Tint>(ListMat[iMat]);
      for (int j = 0; j < m; j++) {
        vals[j][iMat] =
            UniversalScalarConversion<Tint, int64_t>(W.row(j).dot(SHV64.row(j)));
      }
    }
    return vals;
  }
  for (size_t iMat = 0; iMat < nbMat; iMat++) {
    MyMatrix<Tint> W = SHV * ListMat[iMat];
    for (int j = 0; j < m; j++) {
      vals[j][iMat] = W.row(j).dot(SHV.row(j));
    }
  }
  return vals;
}

/*
  A smaller invariant family for a tuple of forms. Every common
  automorphism, and every isometry between two tuples, preserves the values
  (v A_i v^T)_i of a vector, so the family of the first form splits into
  the classes of equal values, permuted among themselves. The classes are
  taken in the order of their size and values, an invariant order, until
  their union spans the lattice: a family that need not be the whole one,
  but spans, as the searches and the canonical bases want it (a family of
  full rank only makes Plesken-Souvignier meet rational leaves). The whole
  family is returned when no proper union spans.
 */
template <typename Tint>
MyMatrix<Tint> RefineFamilyByForms(std::vector<MyMatrix<Tint>> const &ListMat,
                                   MyMatrix<Tint> const &SHVhalf,
                                   std::ostream &os) {
  int n = SHVhalf.cols();
  if (ListMat.size() <= 1 || SHVhalf.rows() == 0) {
    return SHVhalf;
  }
  std::vector<std::vector<Tint>> vals = FormValuesOnFamily(ListMat, SHVhalf);
  std::map<std::vector<Tint>, std::vector<int>> classes;
  for (int j = 0; j < SHVhalf.rows(); j++) {
    classes[vals[j]].push_back(j);
  }
  if (classes.size() == 1) {
    return SHVhalf;
  }
  std::vector<std::pair<size_t, std::vector<Tint>>> order;
  for (auto const &kv : classes) {
    order.push_back({kv.second.size(), kv.first});
  }
  std::sort(order.begin(), order.end());
  std::vector<int> rows;
  for (size_t i_cl = 0; i_cl + 1 < order.size(); i_cl++) {
    for (int j : classes[order[i_cl].second]) {
      rows.push_back(j);
    }
    if (static_cast<int>(rows.size()) < n) {
      continue;
    }
    MyMatrix<Tint> sub(rows.size(), n);
    for (size_t u = 0; u < rows.size(); u++) {
      sub.row(u) = SHVhalf.row(rows[u]);
    }
    if (IsCompleteSystem(ComputeFundamentalInvariant(sub, os), n)) {
#ifdef DEBUG_PLESKEN_SOUVIGNIER
      os << "PS: the family refined by the forms, " << sub.rows() << " of "
         << SHVhalf.rows() << " pairs\n";
#endif
      return sub;
    }
  }
  return SHVhalf;
}

/*
  The kind of a matrix for the engine, which compares the products of the
  vectors in one order only (see ps_cand in PleskenSouvignier.h): that is
  complete for a symmetric or an antisymmetric matrix, whose products in
  the other order are the same up to sign. Any other matrix M is replaced
  by the symmetric M + M^T and the antisymmetric M - M^T: an integral g
  preserves M if and only if it preserves both, their sum being 2 M.
 */
enum class PSFormKind { Symmetric, Antisymmetric, General };

template <typename Tint>
PSFormKind PleskenSouvignierFormKind(MyMatrix<Tint> const &M) {
  int n = M.rows();
  bool is_symm = true;
  bool is_antisymm = true;
  for (int i = 0; i < n; i++) {
    if (M(i, i) != 0) {
      is_antisymm = false;
    }
    for (int j = 0; j < i; j++) {
      if (M(i, j) != M(j, i)) {
        is_symm = false;
      }
      if (M(i, j) != -M(j, i)) {
        is_antisymm = false;
      }
    }
  }
  if (is_symm) {
    return PSFormKind::Symmetric;
  }
  if (is_antisymm) {
    return PSFormKind::Antisymmetric;
  }
  return PSFormKind::General;
}

template <typename Tint>
void PleskenSouvignierPushForm(std::vector<MyMatrix<Tint>> &ListMatInt,
                               MyMatrix<Tint> const &M, PSFormKind kind) {
  if (kind == PSFormKind::General) {
    MyMatrix<Tint> Mtr = M.transpose();
    ListMatInt.push_back(M + Mtr);
    ListMatInt.push_back(M - Mtr);
  } else {
    ListMatInt.push_back(M);
  }
}

/*
  The matrices over the ring, each scaled by the positive multiplier that
  makes it integral, and split as above when neither symmetric nor
  antisymmetric. The automorphisms of a form are those of any positive
  multiple of it, so the search can run on integral forms. Over an
  algebraic field, this needs forms with rational entries.
 */
template <typename T, typename Tint>
std::vector<MyMatrix<Tint>>
PleskenSouvignierIntegralForms(std::vector<MyMatrix<T>> const &ListMat) {
  std::vector<MyMatrix<Tint>> ListMatInt;
  for (auto &eMat : ListMat) {
    MyMatrix<T> eMatScal = RemoveFractionMatrixPlusCoeff(eMat).TheMat;
    if (!IsIntegralMatrix(eMatScal)) {
      std::cerr << "PleskenSouvignierIntegralForms: the Plesken-Souvignier "
                   "method needs forms with rational entries\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(eMatScal);
    PleskenSouvignierPushForm(ListMatInt, M, PleskenSouvignierFormKind(M));
  }
  return ListMatInt;
}

/*
  The same for a pair of configurations: ListMat1[i] and ListMat2[i] are
  scaled by the SAME multiplier, the one making ListMat1[i] integral, since
  P * M1 * P^T = M2 holds if and only if it holds for c M1 and c M2. When
  c ListMat2[i] is not integral, or not of the kind of ListMat1[i] (P M1
  P^T is symmetric, antisymmetric or neither as M1 is), no integral P can
  map the first form to the second, and nothing is returned. The split is
  decided on ListMat1[i] and applied to both.
 */
template <typename T, typename Tint>
std::optional<
    std::pair<std::vector<MyMatrix<Tint>>, std::vector<MyMatrix<Tint>>>>
PleskenSouvignierIntegralFormsPair(std::vector<MyMatrix<T>> const &ListMat1,
                                   std::vector<MyMatrix<T>> const &ListMat2) {
  std::vector<MyMatrix<Tint>> ListMatInt1;
  std::vector<MyMatrix<Tint>> ListMatInt2;
  for (size_t i = 0; i < ListMat1.size(); i++) {
    FractionMatrix<T> fr = RemoveFractionMatrixPlusCoeff(ListMat1[i]);
    if (!IsIntegralMatrix(fr.TheMat)) {
      std::cerr << "PleskenSouvignierIntegralFormsPair: the "
                   "Plesken-Souvignier method needs forms with rational "
                   "entries\n";
      throw TerminalException{1};
    }
    MyMatrix<T> eMat2Scal = fr.TheMult * ListMat2[i];
    if (!IsIntegralMatrix(eMat2Scal)) {
      return {};
    }
    MyMatrix<Tint> M1 = UniversalMatrixConversion<Tint, T>(fr.TheMat);
    MyMatrix<Tint> M2 = UniversalMatrixConversion<Tint, T>(eMat2Scal);
    PSFormKind kind = PleskenSouvignierFormKind(M1);
    if (PleskenSouvignierFormKind(M2) != kind) {
      return {};
    }
    PleskenSouvignierPushForm(ListMatInt1, M1, kind);
    PleskenSouvignierPushForm(ListMatInt2, M2, kind);
  }
  return std::make_pair(std::move(ListMatInt1), std::move(ListMatInt2));
}

/*
  The automorphism group of the configuration: generators over the ring
  in the ORIGINAL basis, satisfying g * M * g^T = M for every matrix of
  ListMat, and the orbit lengths along the stabilizer chain whose
  product is the group order. ListMat[0] has to be symmetric positive
  definite, the others are arbitrary. The matrices may be rational: they
  are made integral, and split when neither symmetric nor antisymmetric,
  by PleskenSouvignierIntegralForms.
 */
template <typename T, typename Tint, typename Tgroup>
PleskenSouvignierAutomResult<Tint, Tgroup>
PleskenSouvignierLatticeAutomorphism(std::vector<MyMatrix<T>> const &ListMat,
                                     InvariantFamilyKind kind,
                                     std::ostream &os) {
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  std::vector<MyMatrix<Tint>> ListMatInt =
      PleskenSouvignierIntegralForms<T, Tint>(ListMat);
  LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(ListMat[0], os);
  MyMatrix<Tint> const &Pmat = rec.Pmat;
  MyMatrix<Tint> PmatInv = Inverse(Pmat);
  std::vector<MyMatrix<Tint>> ListMatRed;
  ListMatRed.reserve(ListMatInt.size());
  for (auto &eMat : ListMatInt) {
    ListMatRed.push_back(Pmat * eMat * Pmat.transpose());
  }
  MyMatrix<Tint> SHVhalf =
      RefineFamilyByForms<Tint>(
          ListMatRed, InvariantFamilyHalf<T, Tint>(rec.GramMatRed, kind, os),
          os);
#ifdef DEBUG_PLESKEN_SOUVIGNIER
  os << "PS: automorphism, family " << NameInvariantFamilyKind(kind) << ", "
     << SHVhalf.rows() << " antipodal pairs\n";
#endif
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  os << "|PS: LatticeAutomorphism, preparation|=" << time << "\n";
#endif
  PleskenSouvignierAutomResult<Tint, Tgroup> result =
      PleskenSouvignierAutomorphism<Tint, Tgroup>(
          ListMatRed, SHVhalf, PSAutomorphismOptions<Tint>{}, os);
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
  every i, in the original bases, or nothing. The two sides build their
  families of the kind independently; the construction is intrinsic, which
  is what makes them correspond under any isometry. diagonal_bound is not
  intrinsic (its bound depends on the reduced basis) and is refused. ListGenAut2 is any
  set of automorphisms of the second configuration in its original
  basis, typically the generators the caller already knows for a lattice
  it tests many candidates against; the orbits they generate prune the
  search, and an empty list is always allowed.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<Tint>> PleskenSouvignierLatticeIsometry(
    std::vector<MyMatrix<T>> const &ListMat1,
    std::vector<MyMatrix<T>> const &ListMat2,
    std::vector<MyMatrix<Tint>> const &ListGenAut2, InvariantFamilyKind kind,
    std::ostream &os) {
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  if (kind == InvariantFamilyKind::DiagonalBound) {
    std::cerr << "PleskenSouvignierLatticeIsometry: the diagonal_bound family "
                 "depends on the basis, it cannot be used on two lattices\n";
    throw TerminalException{1};
  }
  if (ListMat1.size() != ListMat2.size()) {
    std::cerr << "PleskenSouvignierLatticeIsometry: the two lists of "
                 "matrices have different lengths\n";
    throw TerminalException{1};
  }
  // Before the LLL, which cannot run on forms with irrational entries.
  std::optional<
      std::pair<std::vector<MyMatrix<Tint>>, std::vector<MyMatrix<Tint>>>>
      opt_int = PleskenSouvignierIntegralFormsPair<T, Tint>(ListMat1, ListMat2);
  if (!opt_int) {
    return {};
  }
  LLLreduction<T, Tint> rec1 = LLLreducedBasis<T, Tint>(ListMat1[0], os);
  LLLreduction<T, Tint> rec2 = LLLreducedBasis<T, Tint>(ListMat2[0], os);
  MyMatrix<Tint> const &Pmat1 = rec1.Pmat;
  MyMatrix<Tint> const &Pmat2 = rec2.Pmat;
  std::vector<MyMatrix<Tint>> const &ListMatInt1 = opt_int->first;
  std::vector<MyMatrix<Tint>> const &ListMatInt2 = opt_int->second;
  std::vector<MyMatrix<Tint>> ListMatRed1;
  std::vector<MyMatrix<Tint>> ListMatRed2;
  for (auto &eMat : ListMatInt1) {
    ListMatRed1.push_back(Pmat1 * eMat * Pmat1.transpose());
  }
  for (auto &eMat : ListMatInt2) {
    ListMatRed2.push_back(Pmat2 * eMat * Pmat2.transpose());
  }
  MyMatrix<Tint> SHVhalf1 =
      RefineFamilyByForms<Tint>(
          ListMatRed1,
          InvariantFamilyHalf<T, Tint>(rec1.GramMatRed, kind, os), os);
  MyMatrix<Tint> SHVhalf2 =
      RefineFamilyByForms<Tint>(
          ListMatRed2,
          InvariantFamilyHalf<T, Tint>(rec2.GramMatRed, kind, os), os);
#ifdef DEBUG_PLESKEN_SOUVIGNIER
  os << "PS: isometry, families " << NameInvariantFamilyKind(kind)
     << ", pairs " << SHVhalf1.rows() << " / " << SHVhalf2.rows() << "\n";
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
      ListMatRed1, SHVhalf1, ListMatRed2, SHVhalf2, ListGenAut2Red,
      PSIsometryOptions<Tint>{}, os);
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
