// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_LATTICEPLESKENSOUVIGNIER_H_
#define SRC_LATT_LATTICEPLESKENSOUVIGNIER_H_

// clang-format off
#include "PleskenSouvignier.h"
#include "Shvec_exact.h"
#include "ClassicLLL.h"
#include "InvariantVectorFamily.h"
#include "ResidueArithmetic.h"
#include "factorizations.h"
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
  src_group/PleskenSouvignier.h: they LLL-reduce the Gram matrix, build a
  Z-spanning characteristic vector family of it
  (PleskenSouvignierZspanFamilyHalf), run the engine in the reduced basis,
  and transport the answers back to the original basis. The short vector
  enumeration and the LLL are what the engine header cannot depend on,
  which is why these wrappers live in src_latt.

  The vector set of the paper -- the vectors of norm at most the largest
  diagonal entry of the reduced Gram matrix, which contains the standard
  basis and is invariant -- is no longer what the engine is run on; it
  remains below (PleskenSouvignierVectorFamily) for
  LATT_GenerateCharacteristicVectorSet.
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

/*
  Whether the rows of SHV span Z^n. n independent rows span a sublattice of
  index D = |det|; the index of the span of all the rows divides D, and a
  prime p divides it if and only if the rows have rank below n modulo p.
  No Hermite normal form of the family is computed.
 */
template <typename Tint>
bool PleskenSouvignierRowsSpan(MyMatrix<Tint> const &SHV) {
  using Tfield = typename overlying_field<Tint>::field_type;
  int n = SHV.cols();
  int m = SHV.rows();
  // independent rows, by an echelon form over Q kept row by row
  std::vector<MyVector<Tfield>> ech;
  std::vector<int> piv;
  std::vector<int> sel;
  for (int j = 0; j < m && static_cast<int>(sel.size()) < n; j++) {
    MyVector<Tfield> v = UniversalVectorConversion<Tfield, Tint>(
        MyVector<Tint>(SHV.row(j).transpose()));
    for (size_t k = 0; k < ech.size(); k++) {
      if (v(piv[k]) != 0) {
        Tfield c = v(piv[k]) / ech[k](piv[k]);
        v -= c * ech[k];
      }
    }
    int pos = -1;
    for (int i = 0; i < n; i++) {
      if (v(i) != 0) {
        pos = i;
        break;
      }
    }
    if (pos >= 0) {
      ech.push_back(v);
      piv.push_back(pos);
      sel.push_back(j);
    }
  }
  if (static_cast<int>(sel.size()) < n) {
    return false;
  }
  MyMatrix<Tint> Bsel(n, n);
  for (int i = 0; i < n; i++) {
    Bsel.row(i) = SHV.row(sel[i]);
  }
  Tint D = T_abs(DeterminantMat(Bsel));
  if (D == 1) {
    return true;
  }
  for (auto const &p_T : FactorsInt(D)) {
    int64_t p = ResidueModulusOf(p_T);
    // the rank modulo p, by an echelon form kept row by row
    std::vector<std::vector<int64_t>> echp;
    std::vector<int> pivp;
    for (int j = 0; j < m && static_cast<int>(echp.size()) < n; j++) {
      std::vector<int64_t> v(n);
      for (int i = 0; i < n; i++) {
        v[i] = ResModOf(SHV(j, i), p);
      }
      for (size_t k = 0; k < echp.size(); k++) {
        int64_t c = v[pivp[k]];
        if (c != 0) {
          for (int i = 0; i < n; i++) {
            v[i] = ResSub(v[i], ResMul(c, echp[k][i], p), p);
          }
        }
      }
      int pos = -1;
      for (int i = 0; i < n; i++) {
        if (v[i] != 0) {
          pos = i;
          break;
        }
      }
      if (pos >= 0) {
        int64_t inv = ResInv(v[pos], p);
        for (int i = 0; i < n; i++) {
          v[i] = ResMul(v[i], inv, p);
        }
        echp.push_back(v);
        pivp.push_back(pos);
      }
    }
    if (static_cast<int>(echp.size()) < n) {
      return false;
    }
  }
  return true;
}

/*
  The vectors of norm at most B, one per antipodal pair, for the least B
  for which they span the lattice, and that B. It is an isometry invariant,
  so two isometric lattices get the same B and families of the same size,
  and at most the largest diagonal entry of the reduced Gram matrix GramMat,
  whose basis vectors have norm at most that: 2 for E8^3 and the Niemeier
  lattices whose roots span, 4 for the Leech lattice and for A2^12. The
  norms are multiples of g, the gcd of the G_ii and the 2 G_ij, and every
  level from the minimum up is tried.
 */
template <typename T, typename Tint>
std::pair<T, MyMatrix<Tint>>
PleskenSouvignierSpanningFamily(MyMatrix<T> const &GramMat, std::ostream &os) {
  if constexpr (!is_implementation_of_Q<T>::value &&
                !is_implementation_of_Z<T>::value) {
    // over a real field the norms are not discrete
    std::cerr << "PS: the spanning level needs a rational Gram matrix\n";
    throw TerminalException{1};
  } else {
    int n = GramMat.rows();
    T maxdiag = PleskenSouvignierBound(GramMat);
    // g = gcd / d, with d the lcm of the denominators
    Tint d(1);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j <= i; j++) {
        T x = (i == j) ? GramMat(i, i) : T(2) * GramMat(i, j);
        Tint den = UniversalScalarConversion<Tint, T>(GetDenominator(x));
        d = LCMpair(d, den);
      }
    }
    T d_T = UniversalScalarConversion<T, Tint>(d);
    Tint gnum(0);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j <= i; j++) {
        T x = (i == j) ? GramMat(i, i) : T(2) * GramMat(i, j);
        gnum = GcdPair(gnum, UniversalScalarConversion<Tint, T>(x * d_T));
      }
    }
    T step = UniversalScalarConversion<T, Tint>(T_abs(gnum)) / d_T;
    Tshortest<T, Tint> sv = T_ShortestVectorHalf<T, Tint>(GramMat, os);
    int prev_rows = -1;
    for (T B = sv.min; B < maxdiag; B += step) {
      MyMatrix<Tint> fam =
          PleskenSouvignierVectorFamily<T, Tint>(GramMat, B, os);
      if (fam.rows() != prev_rows && PleskenSouvignierRowsSpan(fam)) {
        return {B, fam};
      }
      prev_rows = fam.rows();
    }
    return {maxdiag,
            PleskenSouvignierVectorFamily<T, Tint>(GramMat, maxdiag, os)};
  }
}

/*
  The Z-spanning characteristic vector family fed to the Plesken-Souvignier
  backtrack, one vector per antipodal pair. A Z-spanning family is what the
  engine needs to work reliably: a searched basis drawn from it is unimodular,
  so the backtrack has no rational leaves. The construction is intrinsic (no
  norm bound), hence an isometry invariant, so the two sides of an isometry
  test build corresponding families independently. The root variant is used
  for a single form, where the norm-2 roots are a meaningful block; for
  several forms at once, where "root" is not a notion of the pair, the plain
  shortest-vector tower is used instead. GramMat is the LLL-reduced Gram
  matrix of the callers, which is not reduced again.
 */
template <typename T, typename Tint>
MyMatrix<Tint> PleskenSouvignierZspanFamilyHalf(MyMatrix<T> const &GramMat,
                                                bool single, std::ostream &os) {
  int n = GramMat.rows();
  MyMatrix<Tint> full =
      single ? inner_span_root_iterated_shortest<T, Tint>(
                   GramMat, GramReduction::AlreadyReduced, os)
             : inner_span_iterated_shortest<T, Tint>(
                   GramMat, GramReduction::AlreadyReduced, os);
  std::unordered_set<MyVector<Tint>> seen;
  std::vector<MyVector<Tint>> rows;
  for (int i = 0; i < full.rows(); i++) {
    MyVector<Tint> c = SignCanonicalizeVector(GetMatrixRow(full, i));
    if (seen.insert(c).second) {
      rows.push_back(c);
    }
  }
  return MatrixFromVectorFamilyDim(n, rows);
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
  The vector family of an automorphism computation. Both are invariant under
  the automorphisms and span the lattice: ZspanClosure the iterated
  shortest closure (PleskenSouvignierZspanFamilyHalf), SpanningLevel the
  vectors of norm at most the least level that spans
  (PleskenSouvignierSpanningFamily). The closure can be far larger when the
  short vectors span a sublattice of large index (12.7 10^6 vectors for
  A2^12, against the 1.9 10^5 of norm at most 4).
 */
enum class PSVectorFamily { ZspanClosure, SpanningLevel };

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
PleskenSouvignierLatticeAutomorphismFamily(
    std::vector<MyMatrix<T>> const &ListMat, PSVectorFamily family,
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
  bool single = (ListMat.size() == 1);
  MyMatrix<Tint> SHVhalf =
      (family == PSVectorFamily::ZspanClosure)
          ? PleskenSouvignierZspanFamilyHalf<T, Tint>(rec.GramMatRed, single,
                                                      os)
          : PleskenSouvignierSpanningFamily<T, Tint>(rec.GramMatRed, os).second;
#ifdef DEBUG_PLESKEN_SOUVIGNIER
  os << "PS: automorphism, family " << static_cast<int>(family)
     << " (single=" << single << "), " << SHVhalf.rows()
     << " antipodal pairs\n";
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

template <typename T, typename Tint, typename Tgroup>
PleskenSouvignierAutomResult<Tint, Tgroup>
PleskenSouvignierLatticeAutomorphism(std::vector<MyMatrix<T>> const &ListMat,
                                     std::ostream &os) {
  return PleskenSouvignierLatticeAutomorphismFamily<T, Tint, Tgroup>(
      ListMat, PSVectorFamily::ZspanClosure, os);
}

/*
  The isometry test: P with P * ListMat1[i] * P^T = ListMat2[i] for
  every i, in the original bases, or nothing. The two sides build their
  Z-spanning families independently; the construction is intrinsic, which
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
  bool single = (ListMat1.size() == 1);
  MyMatrix<Tint> SHVhalf1 =
      PleskenSouvignierZspanFamilyHalf<T, Tint>(rec1.GramMatRed, single, os);
  MyMatrix<Tint> SHVhalf2 =
      PleskenSouvignierZspanFamilyHalf<T, Tint>(rec2.GramMatRed, single, os);
#ifdef DEBUG_PLESKEN_SOUVIGNIER
  os << "PS: isometry, Z-spanning families (single=" << single << "), pairs "
     << SHVhalf1.rows() << " / " << SHVhalf2.rows() << "\n";
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
