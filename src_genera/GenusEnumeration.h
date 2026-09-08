// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSENUMERATION_H_
#define SRC_GENERA_GENUSENUMERATION_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include "LatticeStabEquiCan.h"
#include "InvariantVectorFamily.h"
#include "ClassicLLL.h"
#include "Positivity.h"
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_GENUS_ENUMERATION
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_GENUS_ENUMERATION
#endif

#ifdef TIMINGS
#define TIMINGS_GENUS_ENUMERATION
#endif

/*
  Enumeration of the isometry classes of a genus of even positive definite
  lattices, by Kneser p-neighbours.

  The scheme is NOT the adjacency-decomposition scheme of src_perfect and
  src_delaunay, and deliberately does not use their DataFunc machinery. There
  the orbits are joined by a geometric adjacency relation and the enumeration
  terminates when the orbit list is closed under it. Here every lattice of the
  genus is adjacent to a fixed lattice through some p-neighbour, so there is no
  useful notion of a facet, and the termination criterion is arithmetic: the
  Smith-Minkowski-Siegel mass

      mass(G) = sum over the classes L of the genus of 1 / |Aut(L)| .

  The enumeration stops when the accumulated sum reaches the mass supplied on
  input. That also certifies completeness independently of the enumeration, so
  a bug in the neighbour construction cannot silently return a short list: it
  shows up as a mass that never closes.

  Two facts govern the design.

  (1) The p-neighbour graph is connected on a SPINOR genus, not on a genus. A
  genus with several spinor genera therefore needs one seed lattice per spinor
  genus; with a single seed the enumeration converges to a proper subset and
  the mass falls short. The program reports that rather than pretending to be
  complete, and accepts several seeds on input.

  (2) The deduplication is by canonical form (src_latt/LatticeStabEquiCan.h),
  so testing whether a new neighbour is already known is a hash lookup rather
  than an isometry test against every class found so far. That is the point of
  doing this here: the pairwise scheme is quadratic in the class number.
 */

// The description of the genus. `prime` is the prime used for the neighbour
// step; 0 means "choose the smallest prime not dividing 2*det".
template <typename T> struct GenusSpec {
  int rank;
  T det;
  int prime;
  // The invariant vector family used for the canonical form and for the
  // automorphisms. See GenusInvariantVectorFamily.
  std::string method;
};

// The accumulated result of an enumeration.
template <typename T> struct GenusEnumerationResult {
  std::vector<MyMatrix<T>> ListGram;
  std::vector<T> ListAutOrder;
  T accumulated_mass;
  T target_mass;
  bool complete;
};

// Inverse of a modulo the prime p, by extended Euclid. p is small.
inline int InverseModP(int a, int p) {
  int a_red = ((a % p) + p) % p;
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  if (a_red == 0) {
    std::cerr << "GENUS: InverseModP called on a multiple of p\n";
    throw TerminalException{1};
  }
#endif
  int r0 = p, r1 = a_red, s0 = 0, s1 = 1;
  while (r1 != 0) {
    int q = r0 / r1;
    int r2 = r0 - q * r1;
    int s2 = s0 - q * s1;
    r0 = r1;
    r1 = r2;
    s0 = s1;
    s1 = s2;
  }
  return ((s0 % p) + p) % p;
}

template <typename Tint> int ResidueModP(Tint const &x, int p) {
  Tint p_T(p);
  Tint r = x % p_T;
  int r_int = UniversalScalarConversion<int, Tint>(r);
  return ((r_int % p) + p) % p;
}

// The p-adic valuation of a non-zero integer.
template <typename Tint> int ValuationP(Tint const &x, int p) {
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  if (x == 0) {
    std::cerr << "GENUS: ValuationP called on zero\n";
    throw TerminalException{1};
  }
#endif
  Tint p_T(p);
  Tint work = x;
  if (work < 0) {
    work = -work;
  }
  int val = 0;
  while (work % p_T == 0) {
    work /= p_T;
    val++;
  }
  return val;
}

/*
  The points of the projective space over F_p, as the vectors whose first
  non-zero coordinate is 1. There are (p^n - 1) / (p - 1) of them: 16383 for
  n = 14 and p = 2, which is why p = 2 or 3 is what one uses in rank 14.
 */
inline std::vector<std::vector<int>> EnumerateProjectiveLines(int n, int p) {
  std::vector<std::vector<int>> ListLines;
  for (int ipos = 0; ipos < n; ipos++) {
    int n_free = n - ipos - 1;
    int64_t n_case = 1;
    for (int i = 0; i < n_free; i++) {
      n_case *= p;
    }
    for (int64_t idx = 0; idx < n_case; idx++) {
      std::vector<int> v(n, 0);
      v[ipos] = 1;
      int64_t work = idx;
      for (int i = ipos + 1; i < n; i++) {
        v[i] = static_cast<int>(work % p);
        work /= p;
      }
      ListLines.push_back(v);
    }
  }
  return ListLines;
}

/*
  A key for a projective point: normalise so that the first non-zero
  coordinate is 1, then read the coordinates as a base-p number. Injective on
  projective points, which is all that is needed to index them.
 */
inline int64_t EncodeProjectiveLine(std::vector<int> const &v, int p) {
  size_t n = v.size();
  size_t ipos = 0;
  while (ipos < n && ((v[ipos] % p) + p) % p == 0) {
    ipos++;
  }
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  if (ipos == n) {
    std::cerr << "GENUS: EncodeProjectiveLine called on the zero vector\n";
    throw TerminalException{1};
  }
#endif
  int scal = InverseModP(v[ipos], p);
  int64_t code = 0;
  int64_t mult = 1;
  for (size_t i = 0; i < n; i++) {
    int val = (((v[i] % p) + p) % p) * scal;
    val = ((val % p) + p) % p;
    code += mult * val;
    mult *= p;
  }
  return code;
}

/*
  The p-neighbour of the lattice with Gram matrix G attached to an isotropic
  point v of P^{n-1}(F_p):

      N = { y in L : (y, v) = 0 mod p }  +  Z (v/p) .

  For N to be an even lattice one needs v_p((v,v)) >= e + 2 with e = v_p(2),
  that is (v,v) = 0 mod p^2 for odd p and mod 8 for p = 2. A point that is
  isotropic mod p only gives valuation e + 1 in general, and is corrected here
  by v <- v - a p u for a suitable u; that correction is the whole of the p = 2
  difficulty, and is invisible for odd p.

  Returns nothing when v lies in the radical mod p, or when the valuation is
  too small to be corrected (which can only happen for p = 2).
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
GenusNeighbor(MyMatrix<Tint> const &G, std::vector<int> const &v_line, int p,
              [[maybe_unused]] std::ostream &os) {
  int n = G.rows();
  MyVector<Tint> v(n);
  for (int i = 0; i < n; i++) {
    v(i) = Tint(v_line[i]);
  }
  // w = G v, the linear form y -> (y, v)
  auto compute_w = [&](MyVector<Tint> const &x) -> MyVector<Tint> {
    MyVector<Tint> w(n);
    for (int i = 0; i < n; i++) {
      Tint sum(0);
      for (int j = 0; j < n; j++) {
        sum += G(i, j) * x(j);
      }
      w(i) = sum;
    }
    return w;
  };
  MyVector<Tint> w = compute_w(v);
  int jpiv = -1;
  for (int i = 0; i < n; i++) {
    if (ResidueModP(w(i), p) != 0) {
      jpiv = i;
      break;
    }
  }
  if (jpiv == -1) {
    // v lies in the radical of the form mod p: no neighbour from it
    return {};
  }
  auto norm_of = [&](MyVector<Tint> const &x, MyVector<Tint> const &wx) -> Tint {
    Tint sum(0);
    for (int i = 0; i < n; i++) {
      sum += x(i) * wx(i);
    }
    return sum;
  };
  int e = (p == 2) ? 1 : 0;
  Tint nrm = norm_of(v, w);
  int val = ValuationP(nrm, p);
  if (val <= e) {
    return {};
  }
  if (val == e + 1) {
    // correct v so that the valuation rises to at least e + 2. With
    // a = nrm / (2 p (v, e_jpiv)) taken mod p, the two leading terms of
    // (v - a p e_jpiv, v - a p e_jpiv) = nrm - 2 a p (v, e_jpiv) + a^2 p^2 (..)
    // cancel. Both nrm and 2 p (v,e_jpiv) have valuation e + 1, so the ratio
    // is a p-adic unit and a is well defined.
    Tint den = Tint(2) * Tint(p) * w(jpiv);
    Tint p_T(p);
    Tint num_red = nrm;
    Tint den_red = den;
    for (int i = 0; i < e + 1; i++) {
      num_red /= p_T;
      den_red /= p_T;
    }
    int a = (ResidueModP(num_red, p) *
             InverseModP(ResidueModP(den_red, p), p)) % p;
    v(jpiv) -= Tint(a) * Tint(p);
    w = compute_w(v);
    nrm = norm_of(v, w);
#ifdef SANITY_CHECK_GENUS_ENUMERATION
    if (ValuationP(nrm, p) < e + 2) {
      std::cerr << "GENUS: the correction of the isotropic vector failed, "
                << "valuation is " << ValuationP(nrm, p) << " and should be at "
                << "least " << (e + 2) << "\n";
      throw TerminalException{1};
    }
#endif
  }
  /*
    A basis of L_v = { y : (y, v) = 0 mod p }: with c_i = w_i / w_jpiv mod p,
    the vectors e_i - c_i e_jpiv for i != jpiv, together with p e_jpiv. That is
    a triangular basis of determinant p, so the index is p as it must be.
   */
  int inv_pivot = InverseModP(ResidueModP(w(jpiv), p), p);
  MyMatrix<Tint> Bsub(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Bsub(i, j) = Tint(0);
    }
  }
  for (int i = 0; i < n; i++) {
    if (i == jpiv) {
      Bsub(i, jpiv) = Tint(p);
    } else {
      int c = (ResidueModP(w(i), p) * inv_pivot) % p;
      Bsub(i, i) = Tint(1);
      Bsub(i, jpiv) = -Tint(c);
    }
  }
  /*
    N = L_v + Z (v/p), so p N = p L_v + Z v. Take a Z-basis of the lattice
    generated by the rows of p * Bsub together with v, then divide by p.
   */
  MyMatrix<Tint> Agen(n + 1, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Agen(i, j) = Tint(p) * Bsub(i, j);
    }
  }
  for (int j = 0; j < n; j++) {
    Agen(n, j) = v(j);
  }
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time_zb;
#endif
  MyMatrix<Tint> Bp = GetZbasis(Agen);
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: neighbour GetZbasis|=" << time_zb << "\n";
#endif
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  if (Bp.rows() != n) {
    std::cerr << "GENUS: the neighbour basis has " << Bp.rows()
              << " rows instead of " << n << "\n";
    throw TerminalException{1};
  }
#endif
  // Gram of the neighbour: (Bp/p) G (Bp/p)^T = Bp G Bp^T / p^2
  MyMatrix<Tint> Prod = Bp * G * Bp.transpose();
  Tint p2 = Tint(p) * Tint(p);
  MyMatrix<T> GramN(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
#ifdef SANITY_CHECK_GENUS_ENUMERATION
      if (Prod(i, j) % p2 != 0) {
        std::cerr << "GENUS: the neighbour Gram matrix is not integral at ("
                  << i << "," << j << ")\n";
        throw TerminalException{1};
      }
#endif
      GramN(i, j) = UniversalScalarConversion<T, Tint>(Prod(i, j) / p2);
    }
  }
#ifdef SANITY_CHECK_GENUS_ENUMERATION
  for (int i = 0; i < n; i++) {
    T two(2);
    T q = GramN(i, i) / two;
    if (!IsInteger(q)) {
      std::cerr << "GENUS: the neighbour lattice is not even at " << i << "\n";
      throw TerminalException{1};
    }
  }
  if (DeterminantMat(GramN) !=
      UniversalScalarConversion<T, Tint>(DeterminantMat(G))) {
    std::cerr << "GENUS: the neighbour has a different determinant\n";
    throw TerminalException{1};
  }
#endif
  /*
    Reduce before returning. The basis coming out of GetZbasis and divided by p
    is arbitrary and can be very skewed, and everything downstream -- the
    canonical form and the automorphism group -- starts by enumerating short
    vectors, which is where a skewed basis is punishing. Without this the
    enumeration stalls on the first neighbour.
   */
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time_lll;
#endif
  LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(GramN, os);
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: neighbour LLL|=" << time_lll << "\n";
#endif
  return rec.GramMatRed;
}

/*
  The automorphism data of one lattice: the order, needed for the mass, and
  matrix generators, needed to cut the projective points into orbits. Both are
  derived from a single invariant vector family, since extracting it is the
  expensive part and computing them separately would do it twice.
 */
/*
  The invariant vector family the enumeration works from, together with
  whether it spans Z^n.

  Spanning Z^n is not required. The canonicalization handles a family that
  merely has full rank, by canonicalizing the position of Z^n relative to the
  span of the family under the group preserving it, and a smaller family is
  worth much more than the cheaper Hermite route that a spanning family
  allows: the weight matrix is quadratic in the size of the family and the
  canonical labelling of the graph is worse than quadratic, while the extra
  subspace step is not. So the choices below are ordered by size, not by
  whether they span.

  --- "fullrank": ExtractInvariantVectorFamilyFullRank, the shells of the
      lattice taken until they reach full rank. Cheap when the successive
      minima are close together, catastrophic when they are not, since a
      shell has to be taken whole: on TestData/SlowCanonic/slow_canonic_2 it
      returns 16942 vectors, 16816 of them in the single shell of norm 14.
  --- "cv_fullrank": the characteristic vector set V_cv of Section 2.2 of "A
      canonical form for positive definite matrices", with the coset unions
      of (2.2.5) dropped so that it stops at full rank. Asks for closest
      vectors to a few points rather than for whole shells, so it does not
      care how far apart the minima are, and it avoids the 2^k closest-vector
      computations that the spanning version pays when the minimal vectors
      span a sublattice of index 2^k. 238 and 324 vectors on the two
      SlowCanonic lattices, against 16942 and 2358 for the shells.
  --- "cv": the same but built so as to span Z^n, so 2^k times more expensive
      in the bad case, 8516 vectors on slow_canonic_1. Kept because the
      property is what Definition 1.2.1 asks for, not because it is the fast
      choice here. Note that "cv_fullrank" usually spans anyway: the two
      differ only when the minimal vectors fail to generate the lattice.
  --- "auto": both full rank families, keeping the smaller. Building a family
      costs a fraction of a percent of what is done with it -- measured on
      slow_canonic_1, 0.05 s and 1.23 s to build against 7.3 s and more for
      the weight matrix alone -- so the cheapest way to avoid the bad case of
      either is to build both and throw one away.

  The choice made by "auto" depends only on the isometry class, both family
  sizes being invariants, so the canonical forms it produces remain
  comparable across the enumeration.
 */
/*
  Whether the family generates Z^n and not merely a finite index subgroup.
  This has to be measured rather than assumed: a family built without asking
  for it often spans anyway, and every family that does can take the cheap
  Hermite canonicalization instead of the subspace one. Dropping the coset
  unions of V_cv, for instance, changes nothing at all on a lattice whose
  minimal vectors already span, and on the determinant-243 genus that is
  almost every class.

  The test is a Z-basis of the span followed by its determinant, which costs
  nothing next to the weight matrix it saves.
 */
template <typename Tint>
bool GenusFamilySpansLattice(MyMatrix<Tint> const &SHV) {
  int n = SHV.cols();
  if (RankMat(SHV) != n) {
    return false;
  }
  MyMatrix<Tint> Basis = GetZbasis(SHV);
  Tint det = DeterminantMat(Basis);
  return T_abs(det) == Tint(1);
}

template <typename T, typename Tint> struct GenusVectorFamily {
  MyMatrix<Tint> SHV;
  // Whether SHV generates Z^n and not merely a finite index subgroup. It
  // decides which canonicalization applies, the cheap Hermite one or the one
  // that has to place Z^n relative to the span of the family.
  bool spans_lattice;
};

template <typename T, typename Tint>
GenusVectorFamily<T, Tint>
GenusInvariantVectorFamily(MyMatrix<T> const &GramMat,
                           std::string const &method, std::ostream &os) {
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time;
#endif
  auto f_wrap = [&](MyMatrix<Tint> &&SHV) -> GenusVectorFamily<T, Tint> {
    bool spans = GenusFamilySpansLattice<Tint>(SHV);
    return {std::move(SHV), spans};
  };
  auto f_fullrank = [&]() -> GenusVectorFamily<T, Tint> {
    return f_wrap(ExtractInvariantVectorFamilyFullRank<T, Tint>(GramMat, os));
  };
  auto f_cv_fullrank = [&]() -> GenusVectorFamily<T, Tint> {
    return f_wrap(CharacteristicVectorSetCV<T, Tint>(GramMat, false, false, os));
  };
  auto f_cv = [&]() -> GenusVectorFamily<T, Tint> {
    return f_wrap(CharacteristicVectorSetCV<T, Tint>(GramMat, true, false, os));
  };
  auto f_get = [&]() -> GenusVectorFamily<T, Tint> {
    if (method == "fullrank") {
      return f_fullrank();
    }
    if (method == "cv_fullrank") {
      return f_cv_fullrank();
    }
    if (method == "cv") {
      return f_cv();
    }
    if (method == "auto") {
      GenusVectorFamily<T, Tint> fam_cv = f_cv_fullrank();
      GenusVectorFamily<T, Tint> fam_fr = f_fullrank();
      if (fam_cv.SHV.rows() <= fam_fr.SHV.rows()) {
        return fam_cv;
      }
      return fam_fr;
    }
    std::cerr << "GENUS: unknown invariant vector family method " << method
              << ", allowed are fullrank, cv_fullrank, cv and auto\n";
    throw TerminalException{1};
  };
  GenusVectorFamily<T, Tint> fam = f_get();
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: GenusInvariantVectorFamily|=" << time << "\n";
#endif
#ifdef DEBUG_GENUS_ENUMERATION
  os << "GENUS: family of " << fam.SHV.rows() << " vectors, spanning="
     << fam.spans_lattice << "\n";
#endif
  return fam;
}

template <typename T, typename Tint, typename Tgroup> struct LatticeAutInfo {
  T order;
  std::vector<MyMatrix<Tint>> ListGenMat;
};

template <typename T, typename Tint, typename Tgroup>
LatticeAutInfo<T, Tint, Tgroup>
GetLatticeAutInfo(MyMatrix<T> const &GramMat, std::string const &method,
                  std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time;
#endif
  /*
    The FULL RANK family, which is far smaller than the Z-basis one: on the
    rank-14 lattices here the two are 2358 and 24702 vectors, the gap coming
    from a single index-2 obstruction (see TestData/SlowCanonic).

    Taking the full-rank family is safe here only because the generators come
    from GetIntAutomorphism_ListMat_Vdiag, which restricts to the transformations
    preserving the LATTICE. Reading the order off the raw permutation group of a
    full-rank family instead would count rational isometries that do not
    preserve the lattice: that was measured on one class of the determinant-243
    genus as 713451110400 against the true 356725555200, a factor of two. The
    order is then recovered from the permutations induced by the integral
    generators, for which the family only has to be full rank -- a form
    preserving map fixing a full rank family pointwise is the identity.
   */
  MyMatrix<Tint> SHV =
      GenusInvariantVectorFamily<T, Tint>(GramMat, method, os).SHV;
#ifdef DEBUG_GENUS_ENUMERATION
  // The size of this family governs the cost of everything downstream, so it
  // is worth seeing when a lattice is expensive and why. Two things blow it
  // up: a family of full rank that does not span the lattice (see
  // TestData/SlowCanonic/slow_canonic_1.txt, index 2, 2358 against 24702
  // vectors) and successive minima that are spread out, which forces the
  // enumeration far up the shells before it reaches full rank
  // (slow_canonic_2.txt, rank 7 until norm 14, then 21540 vectors).
  os << "GENUS: automorphism family of " << SHV.rows() << " vectors\n";
#endif
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(SHV);
  int n_row = SHV_T.rows();
  std::vector<T> Vdiag(n_row, T(0));
  std::vector<MyMatrix<T>> ListMat{GramMat};
  /*
    ONE automorphism computation. Calling GetListGenAutomorphism_ListMat_Vdiag
    for the order and then ArithmeticAutomorphismGroup for the generators would
    run it twice, since the latter calls the former internally and adds the
    integrality work. Instead the matrix generators are obtained once and the
    order is recovered from the permutation they induce on SHV, which is a
    lookup per row.
   */
  std::vector<MyMatrix<T>> LGen_T =
      GetIntAutomorphism_ListMat_Vdiag<T, Tgroup>(SHV_T, ListMat, Vdiag, os);
  LatticeAutInfo<T, Tint, Tgroup> info;
  for (auto &M_T : LGen_T) {
    info.ListGenMat.push_back(UniversalMatrixConversion<Tint, T>(M_T));
  }
  // The generators satisfy g * Gram * g^T = Gram, so on the rows v of SHV the
  // form-preserving action is v -> v g, that is g^T v in column convention.
  std::unordered_map<MyVector<Tint>, Tidx> MapRow;
  for (int i = 0; i < n_row; i++) {
    MapRow[GetMatrixRow(SHV, i)] = static_cast<Tidx>(i);
  }
  std::vector<Telt> ListPermGens;
  for (auto &M : info.ListGenMat) {
    MyMatrix<Tint> Mtr = M.transpose();
    std::vector<Tidx> ePerm(n_row);
    for (int i = 0; i < n_row; i++) {
      MyVector<Tint> w = Mtr * GetMatrixRow(SHV, i);
      auto iter = MapRow.find(w);
#ifdef SANITY_CHECK_GENUS_ENUMERATION
      if (iter == MapRow.end()) {
        std::cerr << "GENUS: an automorphism does not permute the invariant "
                  << "vector family, which contradicts its invariance\n";
        throw TerminalException{1};
      }
#endif
      ePerm[i] = iter->second;
    }
    ListPermGens.push_back(Telt(ePerm));
  }
  Tgroup grp(ListPermGens, n_row);
  auto siz = grp.size();
  info.order = UniversalScalarConversion<T, decltype(siz)>(siz);
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: GetLatticeAutInfo|=" << time << "\n";
#endif
  return info;
}

/*
  The canonical Gram matrix, used as the key of the dictionary of classes.

  With method "fullrank" this is ComputeCanonicalFormFullRank, which works
  from the full rank invariant vector family rather than from one spanning
  Z^n. That distinction is what makes the canonical form usable at all on
  some lattices: the Z-basis family of slow_canonic_1 has 24702 members
  against 2358 for the full rank one, and the canonical labelling of a graph
  on 24702 vertices did not terminate in 110 s and exhausted memory inside
  nauty. The price is the subspace canonicalization needed to place Z^n
  relative to the lattice the family spans.

  With method "cv" this is ComputeCanonicalFormCV, which spans Z^n and so
  skips that subspace step entirely.

  The two canonical forms differ, as does that of ComputeCanonicalForm, so
  they must never be mixed within one enumeration.
 */
template <typename T, typename Tint, typename Tgroup>
MyMatrix<T> GenusCanonicalGram(MyMatrix<T> const &GramMat,
                               std::string const &method, std::ostream &os) {
  GenusVectorFamily<T, Tint> fam =
      GenusInvariantVectorFamily<T, Tint>(GramMat, method, os);
  auto get_basis = [&]() -> MyMatrix<Tint> {
    if (fam.spans_lattice) {
      // The cheap Hermite route, no subspace canonicalization needed.
      return ComputeCanonicalFormSpanning_family<T, Tint>(GramMat, fam.SHV, os);
    }
    return ComputeCanonicalFormFullRank_family<T, Tint, Tgroup>(GramMat,
                                                                fam.SHV, os);
  };
  MyMatrix<Tint> B = get_basis();
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  return B_T * GramMat * B_T.transpose();
}

/*
  The projective points to visit, one per orbit of the automorphism group.
  Without this the rank-14 enumeration would canonicalise 16383 neighbours per
  lattice at p = 2; the orbits cut that by the order of the group acting.
 */
template <typename Tint>
std::vector<std::vector<int>>
ProjectiveLineOrbitRepresentatives(int n, int p,
                                   std::vector<MyMatrix<Tint>> const &ListGen) {
  std::vector<std::vector<int>> ListLines = EnumerateProjectiveLines(n, p);
  size_t n_line = ListLines.size();
  if (ListGen.empty()) {
    return ListLines;
  }
  std::unordered_map<int64_t, size_t> MapIndex;
  for (size_t i_line = 0; i_line < n_line; i_line++) {
    MapIndex[EncodeProjectiveLine(ListLines[i_line], p)] = i_line;
  }
  // the permutation of the projective points induced by each generator
  std::vector<std::vector<size_t>> ListPermLine;
  for (auto &eGen : ListGen) {
    std::vector<size_t> ePerm(n_line);
    for (size_t i_line = 0; i_line < n_line; i_line++) {
      std::vector<int> const &v = ListLines[i_line];
      std::vector<int> w(n, 0);
      for (int j = 0; j < n; j++) {
        int64_t sum = 0;
        for (int i = 0; i < n; i++) {
          sum += static_cast<int64_t>(v[i]) *
                 UniversalScalarConversion<int64_t, Tint>(eGen(i, j));
        }
        w[j] = static_cast<int>(((sum % p) + p) % p);
      }
      ePerm[i_line] = MapIndex.at(EncodeProjectiveLine(w, p));
    }
    ListPermLine.push_back(ePerm);
  }
  std::vector<uint8_t> status(n_line, 0);
  std::vector<std::vector<int>> ListRepr;
  for (size_t i_line = 0; i_line < n_line; i_line++) {
    if (status[i_line] == 1) {
      continue;
    }
    ListRepr.push_back(ListLines[i_line]);
    std::vector<size_t> ListActive{i_line};
    status[i_line] = 1;
    while (!ListActive.empty()) {
      size_t pos = ListActive.back();
      ListActive.pop_back();
      for (auto &ePerm : ListPermLine) {
        size_t img = ePerm[pos];
        if (status[img] == 0) {
          status[img] = 1;
          ListActive.push_back(img);
        }
      }
    }
  }
  return ListRepr;
}

/*
  The smallest prime not dividing the determinant. Taking 2 whenever the
  determinant is odd matters a great deal: the number of projective points is
  (p^n - 1)/(p - 1), so in rank 16 the choice p = 2 gives 65535 of them while
  p = 3 gives 21523360, and in rank 14 it is 16383 against 2391484. The p = 2
  case needs the valuation correction in GenusNeighbor, which is exactly why
  that correction is there.

  Connectivity of the p-neighbour graph on a spinor genus is the reason for
  excluding the primes dividing the determinant. It is not relied upon blindly:
  if the graph fails to be connected the accumulated mass simply stops short of
  the target and the program reports the enumeration as incomplete.
 */
template <typename T> int ChooseNeighborPrime(T const &det) {
  std::vector<int> ListPrime{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
  for (auto &p : ListPrime) {
    T p_T(p);
    T quot = det / p_T;
    if (!IsInteger(quot)) {
      return p;
    }
  }
  std::cerr << "GENUS: no suitable neighbouring prime found for det=" << det
            << "\n";
  throw TerminalException{1};
}

/*
  The enumeration itself. Starts from the seed lattices, forms p-neighbours,
  keeps the new ones by canonical form, and stops when the accumulated mass
  reaches the target.
 */
template <typename T, typename Tint, typename Tgroup>
GenusEnumerationResult<T>
GenusEnumeration(std::vector<MyMatrix<T>> const &ListSeed, T const &TotalMass,
                 int prime, std::string const &method, std::ostream &os) {
#ifdef TIMINGS_GENUS_ENUMERATION
  MicrosecondTime time_total;
#endif
  int n = ListSeed[0].rows();
  GenusEnumerationResult<T> result;
  result.target_mass = TotalMass;
  result.accumulated_mass = T(0);
  result.complete = false;
  // The classes are held in canonical form, so recognising an already known
  // one is a dictionary lookup rather than a comparison against every class
  // found so far. That is the whole point of canonicalising: the pairwise
  // scheme is quadratic in the class number.
  std::map<MyMatrix<T>, size_t> MapCanonic;
  std::vector<MyMatrix<T>> ListGramWork;
  // The automorphism data is kept alongside: it is the expensive part, it is
  // needed twice (the order for the mass, the generators for the orbits of
  // projective points), and recomputing it in the neighbour loop doubles the
  // cost of the whole enumeration.
  std::vector<LatticeAutInfo<T, Tint, Tgroup>> ListAutInfo;

  auto f_insert = [&](MyMatrix<T> const &GramMat) -> bool {
    MyMatrix<T> GramCan =
        GenusCanonicalGram<T, Tint, Tgroup>(GramMat, method, os);
    if (MapCanonic.count(GramCan) == 1) {
      return false;
    }
    MapCanonic[GramCan] = result.ListGram.size();
    LatticeAutInfo<T, Tint, Tgroup> info =
        GetLatticeAutInfo<T, Tint, Tgroup>(GramCan, method, os);
    result.ListGram.push_back(GramCan);
    result.ListAutOrder.push_back(info.order);
    result.accumulated_mass += T(1) / info.order;
    ListGramWork.push_back(GramCan);
    ListAutInfo.push_back(info);
#ifdef DEBUG_GENUS_ENUMERATION
    os << "GENUS: class " << result.ListGram.size() << ", |Aut|=" << info.order
       << ", mass " << result.accumulated_mass << " / " << TotalMass << "\n";
#endif
    return true;
  };

  for (auto &eSeed : ListSeed) {
    if (result.accumulated_mass == TotalMass) {
      break;
    }
    (void)f_insert(eSeed);
  }
  size_t pos = 0;
  while (pos < ListGramWork.size() && result.accumulated_mass != TotalMass) {
    MyMatrix<T> GramMat = ListGramWork[pos];
    std::vector<MyMatrix<Tint>> const &ListGenMat = ListAutInfo[pos].ListGenMat;
    pos++;
    MyMatrix<Tint> G_int = UniversalMatrixConversion<Tint, T>(GramMat);
#ifdef TIMINGS_GENUS_ENUMERATION
    MicrosecondTime time_orb;
#endif
    std::vector<std::vector<int>> ListRepr =
        ProjectiveLineOrbitRepresentatives<Tint>(n, prime, ListGenMat);
#ifdef TIMINGS_GENUS_ENUMERATION
    os << "|GENUS: ProjectiveLineOrbits|=" << time_orb << "\n";
#endif
#ifdef DEBUG_GENUS_ENUMERATION
    os << "GENUS: lattice " << pos << " has " << ListRepr.size()
       << " orbits of projective points\n";
#endif
    [[maybe_unused]] size_t i_orbit = 0;
    for (auto &eLine : ListRepr) {
      i_orbit++;
#ifdef DEBUG_GENUS_ENUMERATION
      os << "GENUS: lattice " << pos << " orbit " << i_orbit << " / "
         << ListRepr.size() << "\n";
#endif
      std::optional<MyMatrix<T>> opt =
          GenusNeighbor<T, Tint>(G_int, eLine, prime, os);
      if (opt) {
#ifdef DEBUG_GENUS_ENUMERATION
        os << "GENUS:   neighbour built, testing against known classes\n";
#endif
        (void)f_insert(*opt);
        if (result.accumulated_mass == TotalMass) {
          break;
        }
      }
    }
  }
  result.complete = (result.accumulated_mass == TotalMass);
#ifdef TIMINGS_GENUS_ENUMERATION
  os << "|GENUS: GenusEnumeration|=" << time_total << "\n";
#endif
  return result;
}

//
// Input and output
//

/*
  The genus description, as key/value lines:

      rank 14
      det 243
      prime 3        (optional, 0 or absent means "choose automatically")

  It carries no information the enumeration needs beyond the prime -- the
  p-neighbour construction preserves the genus by itself -- but it lets the
  seed lattices be checked against the genus they are supposed to represent,
  which catches a mismatched input immediately rather than as a mass that
  never closes.
 */
template <typename T>
GenusSpec<T> ReadGenusSpecFile(std::string const &file_name) {
  if (!FILE_IsExistingFile(file_name)) {
    std::cerr << "GENUS: the file " << file_name << " does not exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(file_name);
  GenusSpec<T> spec;
  spec.rank = -1;
  spec.det = T(0);
  spec.prime = 0;
  spec.method = "auto";
  std::string key;
  while (is >> key) {
    if (key == "rank") {
      is >> spec.rank;
    } else if (key == "det") {
      is >> spec.det;
    } else if (key == "prime") {
      is >> spec.prime;
    } else if (key == "method") {
      is >> spec.method;
    } else {
      std::cerr << "GENUS: unrecognised key \"" << key << "\" in "
                << file_name << ". Allowed: rank, det, prime, method\n";
      throw TerminalException{1};
    }
  }
  if (spec.rank <= 0) {
    std::cerr << "GENUS: the rank is missing or not positive in " << file_name
              << "\n";
    throw TerminalException{1};
  }
  if (spec.det <= T(0)) {
    std::cerr << "GENUS: the determinant is missing or not positive in "
              << file_name << "\n";
    throw TerminalException{1};
  }
  if (spec.method != "fullrank" && spec.method != "cv" &&
      spec.method != "cv_fullrank" && spec.method != "auto") {
    std::cerr << "GENUS: the method in " << file_name << " is \""
              << spec.method
              << "\", allowed are fullrank, cv_fullrank, cv and auto\n";
    throw TerminalException{1};
  }
  return spec;
}

// The mass, as "num den" or as a single integer.
template <typename T> T ReadMassFile(std::string const &file_name) {
  if (!FILE_IsExistingFile(file_name)) {
    std::cerr << "GENUS: the file " << file_name << " does not exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(file_name);
  T num, den;
  if (!(is >> num)) {
    std::cerr << "GENUS: failed to read the mass from " << file_name << "\n";
    throw TerminalException{1};
  }
  if (!(is >> den)) {
    den = T(1);
  }
  if (den == T(0)) {
    std::cerr << "GENUS: the denominator of the mass is zero in " << file_name
              << "\n";
    throw TerminalException{1};
  }
  return num / den;
}

template <typename T>
void WriteGenusEnumerationResult(std::ostream &os,
                                 GenusEnumerationResult<T> const &result,
                                 std::string const &OutFormat, int prime) {
  if (OutFormat == "Summary") {
    os << "prime = " << prime << "\n";
    os << "class number = " << result.ListGram.size() << "\n";
    os << "accumulated mass = " << result.accumulated_mass << "\n";
    os << "target mass = " << result.target_mass << "\n";
    os << "complete = " << (result.complete ? "true" : "false") << "\n";
    if (!result.complete) {
      os << "The accumulated mass is short of the target. Either the genus "
         << "has several spinor genera and one seed per spinor genus is "
         << "needed, or the prime is not suitable.\n";
    }
    return;
  }
  if (OutFormat == "CPP") {
    os << result.ListGram.size() << "\n";
    for (auto &eGram : result.ListGram) {
      WriteMatrix(os, eGram);
    }
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(ListGram:=[";
    for (size_t i = 0; i < result.ListGram.size(); i++) {
      if (i > 0) {
        os << ",\n";
      }
      WriteMatrixGAP(os, result.ListGram[i]);
    }
    os << "],\n ListAutOrder:=[";
    for (size_t i = 0; i < result.ListAutOrder.size(); i++) {
      if (i > 0) {
        os << ", ";
      }
      os << result.ListAutOrder[i];
    }
    os << "],\n accumulated_mass:=" << result.accumulated_mass
       << ",\n target_mass:=" << result.target_mass
       << ",\n complete:=" << (result.complete ? "true" : "false")
       << ",\n prime:=" << prime << ");\n";
    return;
  }
  std::cerr << "GENUS: failed to find a matching entry for OutFormat="
            << OutFormat << "\n";
  throw TerminalException{1};
}

#endif  //  SRC_GENERA_GENUSENUMERATION_H_
