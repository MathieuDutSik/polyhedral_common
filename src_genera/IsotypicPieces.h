// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ISOTYPICPIECES_H_
#define SRC_GENERA_ISOTYPICPIECES_H_

// clang-format off
#include "MaximalLattices.h"
#include "GenusSymbol.h"
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <utility>
#include <vector>
// clang-format on

/*
  The isotypic pieces of a finite group G, abelian or not: the even
  positive definite G-lattices P whose rational G-module is rho^c, rho an
  irreducible rational representation, up to isometry of pairs (P, G).

  Three complete searches, each a proof that no piece is missed:
  - IsotypicPiecesCertified, rho absolutely irreducible (End_G(M) = Q):
    Lambda' = End(M), starts M (x) X / p^e, steps down;
  - IsotypicPiecesMaximalOrder, End_G(M) with self-adjoint part Q (Schur
    index 2, K^+ = Q): a maximal order Lambda'_p, starts moved up to the
    maximal lattices (unique), steps down;
  - IsotypicPiecesTwisted, End_G(M) a field K with self-adjoint part larger
    than Q: the twists of M'' under the Hermite bound, steps up, steps down.
  A representation none applies to is an error (IsometryGlueFactors.h):
  there is no search that could miss pieces, and thus genera, silently.
  The G-lattices and their steps are in GroupLattice.h, the maximal orders
  and maximal lattices in MaximalLattices.h.
 */
/*
  The key of a start X (genus symbol gs, determinant d): its local symbols
  at the odd primes outside those of G, its 2-adic symbol when 2 is not a
  prime of G, and at each prime p of G the rational class of X over Q_p,
  the square class of d (valuation mod 2, the unit part modulo p, or
  modulo 8 at 2) and the Hasse invariant.
 */
template <typename T, typename Tint>
std::vector<int64_t> StartKey(GenusSymbol const &gs, Tint const &d,
                              std::vector<Tint> const &group_primes,
                              std::ostream &os) {
  std::set<int64_t> gp;
  for (auto const &p : group_primes) {
    gp.insert(UniversalScalarConversion<int64_t, Tint>(p));
  }
  std::vector<int64_t> key{gs.rank, gs.even ? 1 : 0};
  for (auto const &[p, comps] : gs.sym_odd) {
    if (gp.count(p) == 0) {
      key.push_back(p);
      for (auto const &cc : comps) {
        key.insert(key.end(), cc.begin(), cc.end());
      }
    }
  }
  key.push_back(-1);
  if (gp.count(2) == 0) {
    for (auto const &cc : gs.sym2) {
      key.insert(key.end(), cc.begin(), cc.end());
    }
  }
  key.push_back(-1);
  std::vector<T> diag =
      RationalDiagonalForm(GenusSymbolRepresentative<T, Tint>(gs, os));
  for (int64_t p : gp) {
    Tint u = d;
    int v = 0;
    Tint p_i(static_cast<long>(p));
    while (ResInt(u, p_i) == 0) {
      u /= p_i;
      v++;
    }
    int64_t m = (p == 2) ? 8 : p;
    int64_t u_m = UniversalScalarConversion<int64_t, Tint>(ResInt(u, Tint(static_cast<long>(m))));
    int64_t cls = u_m;
    if (p != 2) {
      // the Legendre symbol of u modulo p
      int64_t x = 1;
      for (int64_t k = 0; k < (p - 1) / 2; k++) {
        x = (x * u_m) % p;
      }
      cls = (x == 1) ? 1 : -1;
    }
    key.push_back(p);
    key.push_back(v % 2);
    key.push_back(cls);
    key.push_back(HasseInvariantDiagonal<T, Tint>(diag, p));
  }
  return key;
}

/*
  The certified search, for rho of Schur index 1 (ListGens(Z[G]) of full
  rank n^2 in End(M) = M_n(Z)).

  Let Lambda be the image of Z_p[G] in End(V_rho) and Lambda' = End(M_p) =
  M_n(Z_p), a maximal order containing it, p^a the exponent of Lambda' /
  Lambda (ConductorExponent). For a G-lattice L, N = Lambda' L is a
  Lambda'-lattice and p^a N subset Lambda L = L subset N. By Morita, N = M
  (x) X with the form b_M (x) q, and as p^a N is in L, integral, q has
  values in p^{-2a - v} (v the valuation of the content of b_M): N is an M
  (x) X with its form divided by p^e, X integral, 0 <= e <= 2a + v, and
  v_p det N >= v_p det L - 2 a rank L. L is reached from N by a chain of
  maximal G-sublattices, whose determinants grow from det N to det L. So
  every piece of determinant D is found by the steps down from the
  starts M (x) X / prod p^e_p, X over the genera (odd and even) of rank c
  whose determinant puts that of the start between v_p(D) - 2 a_p rank and
  v_p(D) at each prime p of G (and right outside them), cutting the lattices
  whose determinant goes beyond the largest admissible valuation: a finite
  and complete search, up to the completeness of the genus symbols and the
  Morita equivalence. One pair per class, by the canonical form (of the
  form scaled to be integral).
 */

// v_p of the exponent of M_n(Z) / span(g, g in G), or nullopt when the
// span is not of full rank n^2.
template <typename T, typename Tint>
std::optional<int> ConductorExponent(GroupLattice<T> const &M, Tint const &p) {
  int n = M.gram.rows();
  std::vector<MyMatrix<T>> ListP;
  for (auto const &g : M.gens) {
    ListP.push_back(g.transpose());
  }
  std::vector<MyMatrix<T>> elts =
      GroupElementsAndInverses(ListP, 100000).first;
  MyMatrix<T> Span(static_cast<int>(elts.size()), n * n);
  for (size_t k = 0; k < elts.size(); k++) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        Span(static_cast<int>(k), i * n + j) = elts[k](i, j);
      }
    }
  }
  MyMatrix<T> B = GetZbasis(Span);
  if (B.rows() != n * n) {
    return {};
  }
  MyMatrix<T> Binv = Inverse(B);
  Tint e(1);
  for (int i = 0; i < n * n; i++) {
    for (int j = 0; j < n * n; j++) {
      e = LCMpair(e, UniversalScalarConversion<Tint, T>(GetDenominator(Binv(i, j))));
    }
  }
  int a = 0;
  while (ResInt(e, p) == 0) {
    e /= p;
    a++;
  }
  return a;
}

// v_p of a nonzero rational.
template <typename T, typename Tint> int RationalValuation(T const &x, Tint const &p) {
  T den = GetDenominator(x);
  T num = x * den;
  Tint num_i = UniversalScalarConversion<Tint, T>(T_abs(num));
  Tint den_i = UniversalScalarConversion<Tint, T>(den);
  int v = 0;
  while (ResInt(num_i, p) == 0) {
    num_i /= p;
    v++;
  }
  while (ResInt(den_i, p) == 0) {
    den_i /= p;
    v--;
  }
  return v;
}

template <typename T, typename Tint, typename Tgroup>
std::optional<std::vector<GroupLattice<T>>>
IsotypicPiecesCertified(GroupLattice<T> const &M, int c,
                        std::vector<Tint> const &dets,
                        std::vector<Tint> const &group_primes,
                        std::ostream &os) {
  int n = M.gram.rows();
  int rank = n * c;
  // per prime of G: a_p, the content and det valuations of M, the range
  std::vector<int> a, v_cont, v_detM, hi, lo;
  for (auto const &p : group_primes) {
    std::optional<int> opt = ConductorExponent<T, Tint>(M, p);
    if (!opt) {
      return {};
    }
    a.push_back(*opt);
    int vc = std::numeric_limits<int>::max();
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (M.gram(i, j) != 0) {
          vc = std::min(vc, RationalValuation<T, Tint>(M.gram(i, j), p));
        }
      }
    }
    v_cont.push_back(vc);
    v_detM.push_back(RationalValuation<T, Tint>(DeterminantMat(M.gram), p));
    int h = 0, l = std::numeric_limits<int>::max();
    for (auto const &D : dets) {
      int vD = RationalValuation<T, Tint>(UniversalScalarConversion<T, Tint>(D), p);
      h = std::max(h, vD);
      l = std::min(l, vD);
    }
    hi.push_back(h);
    lo.push_back(l - 2 * a.back() * rank);
  }
  // the determinants of X outside the primes of G
  // M is integral: its determinant is an integer
  Tint det_M = UniversalScalarConversion<Tint, T>(DeterminantMat(M.gram));
  auto out_part = [&](Tint x) -> Tint {
    for (auto const &p : group_primes) {
      while (ResInt(x, p) == 0) {
        x /= p;
      }
    }
    return x;
  };
  Tint det_M_out = out_part(det_M);
  std::set<Tint> dets_X_out;
  for (auto const &D : dets) {
    Tint D_out = out_part(D);
    Tint M_part = MyPow<Tint>(det_M_out, c);
    if (ResInt(D_out, M_part) != 0) {
      continue;
    }
    Tint rest = D_out / M_part;
    Tint root(1);
    bool ok = true;
    if (rest > 1) {
      for (auto const &[q, e] : FactorsIntMap(rest)) {
        if (static_cast<int>(e) % n != 0) {
          ok = false;
          break;
        }
        root *= MyPow<Tint>(q, static_cast<int>(e) / n);
      }
    }
    if (ok) {
      dets_X_out.insert(root);
    }
  }
  LatticeMethodChoice choice = DefaultPairMethods().canonical;
  std::set<std::vector<MyMatrix<T>>> known;
  std::vector<GroupLattice<T>> queue, res;
  auto insert = [&](GroupLattice<T> const &L) -> void {
    T det = DeterminantMat(L.gram);
    for (size_t ip = 0; ip < group_primes.size(); ip++) {
      if (RationalValuation<T, Tint>(det, group_primes[ip]) > hi[ip]) {
        return;
      }
    }
    // the canonical form of the form scaled to be integral
    T den(1);
    for (int i = 0; i < L.gram.rows(); i++) {
      for (int j = 0; j < L.gram.cols(); j++) {
        den = LCMpair(den, GetDenominator(L.gram(i, j)));
      }
    }
    std::vector<MyMatrix<T>> ListP;
    for (auto const &g : L.gens) {
      ListP.push_back(g.transpose());
    }
    MyMatrix<T> A = den * L.gram;
    MyMatrix<Tint> B = PairCanonicalBasis<T, Tint, Tgroup>(A, ListP, choice, os);
    GroupLattice<T> Lc =
        GroupLatticeInBasis(L, UniversalMatrixConversion<T, Tint>(B));
    std::vector<MyMatrix<T>> key{Lc.gram};
    for (auto const &g : Lc.gens) {
      key.push_back(g);
    }
    if (known.insert(key).second) {
      queue.push_back(Lc);
    }
  };
  // the starts: per prime of G, the scaling e and the valuation of det X
  std::vector<std::pair<Tint, Tint>> choices{{Tint(1), Tint(1)}};
  for (size_t ip = 0; ip < group_primes.size(); ip++) {
    Tint const &p = group_primes[ip];
    std::vector<std::pair<Tint, Tint>> next;
    for (int e = 0; e <= 2 * a[ip] + v_cont[ip]; e++) {
      // v(det N) = c v(det M) + n v(det X) - rank e in [lo, hi]
      for (int vx = 0; c * v_detM[ip] + n * vx - rank * e <= hi[ip]; vx++) {
        if (c * v_detM[ip] + n * vx - rank * e < lo[ip]) {
          continue;
        }
        for (auto const &[sc, dx] : choices) {
          next.push_back({sc * MyPow<Tint>(p, e), dx * MyPow<Tint>(p, vx)});
        }
      }
    }
    choices = next;
  }
  size_t n_X = 0;
  for (auto const &d_out : dets_X_out) {
    for (auto const &[sc, dx] : choices) {
      for (bool even : {true, false}) {
        for (auto const &gs : EnumerateGenusSymbols<Tint>(c, d_out * dx, even)) {
          MyMatrix<T> X = GenusSymbolRepresentative<T, Tint>(gs, os);
          GroupLattice<T> N = TensorWithLattice(M, X);
          N.gram /= UniversalScalarConversion<T, Tint>(sc);
          insert(N);
          n_X++;
        }
      }
    }
  }
  size_t n_start = queue.size();
  for (size_t pos = 0; pos < queue.size(); pos++) {
    GroupLattice<T> L = queue[pos];
    for (auto const &p : group_primes) {
      T_small_prime p_s = UniversalScalarConversion<T_small_prime, Tint>(p);
      size_t n_down = 0;
      std::vector<GroupLattice<T>> neigh = FlipNeighbors(L, p_s, n_down, true, false);
      for (size_t i = 0; i < n_down; i++) {
        insert(neigh[i]);
      }
    }
  }
  for (auto const &L : queue) {
    if (!IsEvenIntegralGram(L.gram)) {
      continue;
    }
    Tint det_L = UniversalScalarConversion<Tint, T>(DeterminantMat(L.gram));
    if (std::binary_search(dets.begin(), dets.end(), det_L)) {
      res.push_back(L);
    }
  }
  os << "ISOM_PIECES: certified, " << n_X << " X, " << n_start << " starts, "
     << queue.size() << " lattices, " << res.size()
     << " pieces of admissible determinant\n";
  return res;
}

/*
  The complete search for any Schur index, through a maximal order.

  B = Q rho(G) in End(V_rho), Lambda = Z rho(G). At a prime p of G,
  MaximalOrderAt finds a maximal order Lambda'_p of B_p containing Lambda:
  if Lambda is not maximal, Lambda + p^{a-1} Lambda' (Lambda' maximal
  containing it, p^a Lambda' in Lambda, a >= 1 minimal) is a ring between
  Lambda and p^{-1} Lambda, so some x in p^{-1} Lambda generates with
  Lambda a ring inside p^{-1} Lambda, larger than Lambda; and none exists
  once Lambda is maximal.

  For a G-lattice L, N = Lambda'_p L satisfies p^a N in L in N, and as p^a
  N in L is integral, the form on N has values in p^{-2a}: N lies in a
  maximal p^{-2a}-integral Lambda'_p-lattice. Those are unique up to the
  isometries of the hermitian space, the centralizer of G (Eichler for the
  split Morita cases, Shimura for the hermitian spaces over a local
  division algebra). So a start whose completion at each prime p of G is
  such a maximal lattice, and whose completions at the other primes are
  those of the pieces (M (x) X, X over the local symbols there, Morita at
  p not dividing |G|), contains a lattice of the genus of every piece,
  reached by the steps down with the determinant cut. For rho of Schur
  index 2 (Albert type III, a definite quaternion algebra with its
  canonical involution) the positive hermitian forms of rank c are all
  isometric, so a single rational class.
 */

/*
  The steps down from the starts with the determinant cut (valuations at
  the primes of G at most hi), one pair per class by the canonical form of
  the form scaled to be integral; the even lattices of determinant in dets.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<GroupLattice<T>>
StepsDownWithCut(std::vector<GroupLattice<T>> const &starts,
                 std::vector<Tint> const &dets,
                 std::vector<Tint> const &group_primes,
                 std::vector<int> const &hi, std::ostream &os) {
  LatticeMethodChoice choice = DefaultPairMethods().canonical;
  std::set<std::vector<MyMatrix<T>>> known;
  std::vector<GroupLattice<T>> queue, res;
  auto insert = [&](GroupLattice<T> const &L) -> void {
    T det = DeterminantMat(L.gram);
    for (size_t ip = 0; ip < group_primes.size(); ip++) {
      if (RationalValuation<T, Tint>(det, group_primes[ip]) > hi[ip]) {
        return;
      }
    }
    T den(1);
    for (int i = 0; i < L.gram.rows(); i++) {
      for (int j = 0; j < L.gram.cols(); j++) {
        den = LCMpair(den, GetDenominator(L.gram(i, j)));
      }
    }
    std::vector<MyMatrix<T>> ListP;
    for (auto const &g : L.gens) {
      ListP.push_back(g.transpose());
    }
    MyMatrix<T> A = den * L.gram;
    MyMatrix<Tint> B = PairCanonicalBasis<T, Tint, Tgroup>(A, ListP, choice, os);
    GroupLattice<T> Lc =
        GroupLatticeInBasis(L, UniversalMatrixConversion<T, Tint>(B));
    std::vector<MyMatrix<T>> key{Lc.gram};
    for (auto const &g : Lc.gens) {
      key.push_back(g);
    }
    if (known.insert(key).second) {
      queue.push_back(Lc);
    }
  };
  for (auto const &N : starts) {
    insert(N);
  }
  size_t n_start = queue.size();
  for (size_t pos = 0; pos < queue.size(); pos++) {
    GroupLattice<T> L = queue[pos];
    for (auto const &p : group_primes) {
      T_small_prime p_s = UniversalScalarConversion<T_small_prime, Tint>(p);
      size_t n_down = 0;
      std::vector<GroupLattice<T>> neigh = FlipNeighbors(L, p_s, n_down, true, false);
      for (size_t i = 0; i < n_down; i++) {
        insert(neigh[i]);
      }
    }
  }
  for (auto const &L : queue) {
    if (!IsEvenIntegralGram(L.gram)) {
      continue;
    }
    Tint det_L = UniversalScalarConversion<Tint, T>(DeterminantMat(L.gram));
    if (std::binary_search(dets.begin(), dets.end(), det_L)) {
      res.push_back(L);
    }
  }
  os << "ISOM_PIECES: " << n_start << " starts, " << queue.size()
     << " lattices, " << res.size() << " pieces of admissible determinant\n";
  return res;
}

/*
  The twists of M: for a in the self-adjoint part S of End_G(M) (K when
  End_G(rho) = K is a totally real field, K^+ for a CM field, Q for a
  quaternion algebra with its canonical involution), M(a) is M with the
  form b_M(a x, y), positive when a is totally positive. A form over K of
  rank c is diagonal, <a_1, ..., a_c>, so the sums of twists M(a_i) reach
  every rational class, which M (x) X, X over Z, does not when S != Q.

  The a are the integral elements of S (a Z-basis of S cap M_n(Z)) with
  Tr(a^2) <= bound, totally positive: a search by a box, the bound making
  room for the norms asked for (up to the unit squares, a reduced a has a
  trace of the order of n N(a)^{1/n}).
 */
// The totally positive a in the integral self-adjoint commutant with
// Tr(a^2) <= bound2, as the twisted Gram matrices a^T A (one per form).
template <typename T, typename Tint>
std::vector<MyMatrix<T>> TwistedGrams(GroupLattice<T> const &M,
                                      std::vector<MyMatrix<T>> const &S,
                                      T const &bound2) {
  int k = S.size();
  MyMatrix<T> Q(k, k);
  for (int i = 0; i < k; i++) {
    for (int j = 0; j < k; j++) {
      MyMatrix<T> P = S[i] * S[j];
      T tr(0);
      for (int u = 0; u < P.rows(); u++) {
        tr += P(u, u);
      }
      Q(i, j) = tr;
    }
  }
  MyMatrix<T> Qinv = Inverse(Q);
  std::vector<int64_t> box(k);
  for (int i = 0; i < k; i++) {
    T r2 = bound2 * Qinv(i, i);
    int64_t r = 0;
    while (UniversalScalarConversion<T, int64_t>((r + 1) * (r + 1)) <= r2) {
      r++;
    }
    box[i] = r;
  }
  std::set<MyMatrix<T>> seen;
  std::vector<MyMatrix<T>> res;
  std::vector<int64_t> v(k);
  std::function<void(int)> rec = [&](int i) -> void {
    if (i == k) {
      MyMatrix<T> X = ZeroMatrix<T>(M.gram.rows(), M.gram.rows());
      for (int u = 0; u < k; u++) {
        X += UniversalScalarConversion<T, int64_t>(v[u]) * S[u];
      }
      T q(0);
      MyMatrix<T> X2 = X * X;
      for (int u = 0; u < X2.rows(); u++) {
        q += X2(u, u);
      }
      if (q == 0 || q > bound2) {
        return;
      }
      MyMatrix<T> Aa = X.transpose() * M.gram;
      if (IsPositiveDefinite(Aa, std::cerr) && seen.insert(Aa).second) {
        res.push_back(Aa);
      }
      return;
    }
    for (int64_t x = -box[i]; x <= box[i]; x++) {
      v[i] = x;
      rec(i + 1);
    }
  };
  rec(0);
  return res;
}

/*
  The pieces by the maximal orders (any Schur index, see the comment
  above): nullopt when a maximal order is beyond the bounds.
 */
template <typename T, typename Tint, typename Tgroup>
std::optional<std::vector<GroupLattice<T>>>
IsotypicPiecesMaximalOrder(GroupLattice<T> const &M, int c,
                           std::vector<Tint> const &dets,
                           std::vector<Tint> const &group_primes,
                           std::ostream &os) {
  int n = M.gram.rows();
  std::vector<MaximalOrderData<T, Tint>> orders;
  std::vector<int> hi;
  for (auto const &p : group_primes) {
    // the maximal lattices are unique at p (MaximalLattices.h), or this
    // search is no proof
    if (!MaximalLatticesUniqueAt<T, Tint>(M, p)) {
      os << "ISOM_PIECES: maximal orders, the maximal lattices are not "
            "unique at p="
         << p << "\n";
      return {};
    }
    std::optional<MaximalOrderData<T, Tint>> opt =
        MaximalOrderAt<T, Tint>(M, p, 1000000);
    if (!opt) {
      return {};
    }
    orders.push_back(*opt);
    int h = 0;
    for (auto const &D : dets) {
      h = std::max(h, RationalValuation<T, Tint>(
                          UniversalScalarConversion<T, Tint>(D), p));
    }
    hi.push_back(h);
    os << "ISOM_PIECES: maximal order at p=" << p << ", a=" << opt->a
       << ", dim " << opt->basis.size() << "\n";
  }
  // the determinants of X outside the primes of G
  auto out_part = [&](Tint x) -> Tint {
    for (auto const &p : group_primes) {
      while (ResInt(x, p) == 0) {
        x /= p;
      }
    }
    return x;
  };
  Tint det_M_out =
      out_part(UniversalScalarConversion<Tint, T>(DeterminantMat(M.gram)));
  std::set<Tint> dets_X_out;
  for (auto const &D : dets) {
    Tint D_out = out_part(D);
    Tint M_part = MyPow<Tint>(det_M_out, c);
    if (ResInt(D_out, M_part) != 0) {
      continue;
    }
    Tint rest = D_out / M_part;
    Tint root(1);
    bool ok = true;
    if (rest > 1) {
      for (auto const &[q, e] : FactorsIntMap(rest)) {
        if (static_cast<int>(e) % n != 0) {
          ok = false;
          break;
        }
        root *= MyPow<Tint>(q, static_cast<int>(e) / n);
      }
    }
    if (ok) {
      dets_X_out.insert(root);
    }
  }
  // one X per local symbol outside the primes of G and rational class at
  // them
  std::vector<Tint> extra{Tint(1)};
  for (auto const &p : group_primes) {
    std::vector<Tint> next;
    for (auto const &x : extra) {
      Tint pe(1);
      for (int e = 0; e <= 3; e++) {
        next.push_back(x * pe);
        pe *= p;
      }
    }
    extra = next;
  }
  std::map<std::vector<int64_t>, std::pair<Tint, GenusSymbol>> starts_X;
  for (auto const &d_out : dets_X_out) {
    for (auto const &e : extra) {
      for (bool even : {true, false}) {
        for (auto const &gs : EnumerateGenusSymbols<Tint>(c, d_out * e, even)) {
          // the local symbols outside the primes of G and the rational
          // class at them (several when End_G(rho) is a field K != Q: a
          // hermitian form <a> over K / K^+ depends on a modulo the local
          // norms)
          std::vector<int64_t> key = StartKey<T, Tint>(gs, d_out * e,
                                                       group_primes, os);
          auto iter = starts_X.find(key);
          if (iter == starts_X.end() || d_out * e < iter->second.first) {
            starts_X[key] = {d_out * e, gs};
          }
        }
      }
    }
  }
  // the starts before the steps at the primes of G, with their layout: M
  // (x) X (coordinates i c + k) or a direct sum of twists (u n + i)
  std::vector<std::pair<GroupLattice<T>, bool>> before;
  for (auto const &[key, val] : starts_X) {
    MyMatrix<T> X = GenusSymbolRepresentative<T, Tint>(val.second, os);
    before.push_back({TensorWithLattice(M, X), true});
  }
  // when S is larger than Q (and the commutant not a field, else the
  // twists are used), the M (x) X do not reach every rational class: no
  // complete search
  if (SelfAdjointCommutantBasis<T, Tint>(M).size() > 1) {
    return {};
  }
  std::vector<GroupLattice<T>> starts;
  for (auto [L, tensor_layout] : before) {
    // L is M (x) Z^c up to the form (UpToMaximalLattice)
    for (size_t ip = 0; ip < group_primes.size(); ip++) {
      L = UpToMaximalLattice<T, Tint>(L, orders[ip], group_primes[ip], n,
                                      tensor_layout);
    }
    starts.push_back(L);
  }
  os << "ISOM_PIECES: maximal orders, " << starts_X.size() << " X\n";
  return StepsDownWithCut<T, Tint, Tgroup>(starts, dets, group_primes, hi, os);
}

/*
  The complete search through the twists, when End_G(M) is a field K (rho
  of Schur index 1) whose self-adjoint part S is larger than Q (K totally
  real != Q, or K^+ != Q).

  M is first replaced by M' = O_K M (O_K the maximal order of the
  commutant, CommutativeMaximalOrder), locally free over O_K of rank n / d, d
  = [K : Q], with b_M' = Tr_{K/Q} h' for a K-valued form h'. Let P be a
  piece with c copies, of determinant at most D. X = Hom_G(M', P) is an
  O_K-module of rank c with the form h, b_P(x m, y m') = Tr h'(h(x, y) m,
  m'), and Q(x) = Tr_{K/Q} h(x, x) is positive on X, of rank R = c d over
  Z. Y = X(M') = X (x)_{O_K} M' has det(Y) = det_Q(X)^{n/d} det(M')^c /
  |d_K|^{c n / d}, and [P : Y] <= f = prod_{p | |G|} p^{2 a_p rank P}: with
  Lambda' a maximal order containing Lambda = Z rho(G), p^a Lambda' in
  Lambda, N = Lambda' P, the Lambda'-maps Lambda' M' -> N restricted to M'
  and multiplied by p^a lie in X, and they cover p^a N (Lambda' M' is a
  progenerator), so Y contains p^{2a} N, thus p^{2a} P. Hence det_Q(X) <=
  (D f^2 |d_K|^{c n/d} / det(M')^c)^{d/n}, and by Hermite some x != 0 has
  Q(x) <= gamma_R det_Q(X)^{1/R}: x(M') is the twist M'(phi), phi = h(x, x),
  with Tr(phi) = (n / d) Q(x) as a matrix. P contains x(M') (+) (P cap
  x(M')^perp) with an index at most det x(M') <= det(M') (Tr phi / n)^n,
  the complement a piece of c - 1 copies of determinant at most D det x(M').
  By induction P is an integral G-overlattice of an orthogonal sum of c
  twists under the bound, reached by steps up at the primes of its
  determinant. The starts, the sums of twists of integral Gram matrices
  A_M' phi under the bound, and their integral overlattices by steps up are
  a complete finite search, with no unit or class group of K.
 */

// gamma_R^R, the Hermite constants known exactly for R <= 8, and
// Hermite's bound (4/3)^{R (R - 1) / 2} beyond.
inline double HermiteConstantPower(int R) {
  static const double known[9] = {1, 1, 4.0 / 3.0, 2, 4, 8, 64.0 / 3.0, 64, 256};
  if (R <= 8) {
    return known[R];
  }
  return std::pow(4.0 / 3.0, R * (R - 1) / 2.0);
}

template <typename T, typename Tint, typename Tgroup>
std::optional<std::vector<GroupLattice<T>>>
IsotypicPiecesTwisted(GroupLattice<T> const &M_in, int c,
                      std::vector<Tint> const &dets, std::ostream &os) {
  if (SelfAdjointCommutantBasis<T, Tint>(M_in).size() <= 1) {
    return {};
  }
  int n = M_in.gram.rows();
  // the commutant: a field K (commutative), else not this search
  std::vector<MyMatrix<T>> C = IntegralCommutantBasis<T, Tint>(M_in);
  for (auto const &x : C) {
    for (auto const &y : C) {
      if (x * y != y * x) {
        return {};
      }
    }
  }
  int d = static_cast<int>(C.size());
  T nd = UniversalScalarConversion<T, int>(n / d);
  auto trace_form = [&](std::vector<MyMatrix<T>> const &B) -> MyMatrix<T> {
    int k = B.size();
    MyMatrix<T> F(k, k);
    for (int a = 0; a < k; a++) {
      for (int b = 0; b < k; b++) {
        MyMatrix<T> P = B[a] * B[b];
        T tr(0);
        for (int u = 0; u < n; u++) {
          tr += P(u, u);
        }
        F(a, b) = tr / nd;
      }
    }
    return F;
  };
  // the lattice spanned by the images of L by the operators (on the
  // columns, in the coordinates of L)
  auto saturate = [&](GroupLattice<T> const &L,
                      std::vector<MyMatrix<T>> const &ops) -> GroupLattice<T> {
    MyMatrix<T> Gen = IdentityMat<T>(L.gram.rows());
    for (auto const &u : ops) {
      Gen = Concatenate(Gen, MyMatrix<T>(u.transpose()));
    }
    return GroupLatticeOfRows(L, Gen);
  };
  // O_K (Round 2), then M' = O_K M
  std::vector<MyMatrix<T>> OK = CommutativeMaximalOrder<T, Tint>(C, n);
  GroupLattice<T> M = saturate(M_in, OK);
  // then M'' = Lambda'_p M' at the primes of G, with the a_p
  std::vector<Tint> gp = GroupLatticePrimes<T, Tint>(M_in);
  T s0(1);
  std::vector<int> a_p, hi;
  for (auto const &p : gp) {
    std::optional<MaximalOrderData<T, Tint>> opt =
        MaximalOrderAt<T, Tint>(M, p, 100000000);
    if (!opt) {
      return {};
    }
    M = saturate(M, opt->basis);
    a_p.push_back(opt->a);
    T p_T = UniversalScalarConversion<T, Tint>(p);
    for (int e = 0; e < 2 * opt->a; e++) {
      s0 *= p_T;
    }
    int h = 0;
    for (auto const &D : dets) {
      h = std::max(h, RationalValuation<T, Tint>(
                          UniversalScalarConversion<T, Tint>(D), p));
    }
    hi.push_back(h);
  }
  std::vector<MyMatrix<T>> S = SelfAdjointCommutantBasis<T, Tint>(M);
  MyMatrix<T> const &A = M.gram;
  MyMatrix<T> Ainv = Inverse(A);
  double dK = UniversalScalarConversion<double, T>(
      T_abs(DeterminantMat(trace_form(IntegralCommutantBasis<T, Tint>(M)))));
  double log_detM =
      std::log(UniversalScalarConversion<double, T>(DeterminantMat(A)));
  double log_s0 = std::log(UniversalScalarConversion<double, T>(s0));
  // the bounds along the induction, in logarithms: det(N) <= D, N =
  // Lambda' P, then det of the complement <= det(N) det(x M'') s0^{2n}
  double logD = 0;
  for (auto const &dd : dets) {
    logD =
        std::max(logD, std::log(UniversalScalarConversion<double, Tint>(dd)));
  }
  double bound_tr = 0;
  double bound_det = 0;
  for (int i = 0; i < c; i++) {
    int ci = c - i;
    int R = ci * d;
    double log_detQX = (static_cast<double>(d) / n) *
                       (logD + (static_cast<double>(ci) * n / d) * std::log(dK) -
                        ci * log_detM);
    double minQ =
        std::exp((std::log(HermiteConstantPower(R)) + log_detQX) / R);
    double tr = (static_cast<double>(n) / d) * minQ;
    bound_tr = std::max(bound_tr, tr);
    double log_det_x = log_detM + n * std::log(tr / n);
    bound_det = std::max(bound_det, std::exp(log_det_x));
    logD += log_det_x + 2 * n * log_s0;
  }
  os << "ISOM_PIECES: twists, [K:Q] = " << d << ", |d_K| = " << dK
     << ", trace bound " << bound_tr << "\n";
  T Btr = UniversalScalarConversion<T, int64_t>(
      static_cast<int64_t>(std::ceil(bound_tr * (1 + 1e-9))) + 1);
  T bound2 = Btr * Btr;
  T det_max = UniversalScalarConversion<T, int64_t>(
      static_cast<int64_t>(std::min(bound_det * (1 + 1e-9) + 1, 9e18)));
  // the Gram matrices A phi, phi in S, of x(M'') in N: in (1/s0) Z
  std::vector<MyMatrix<T>> gramsS;
  for (auto const &Si : S) {
    gramsS.push_back(s0 * A * Si);
  }
  std::vector<MyMatrix<T>> GB =
      UnflattenMatrices(IntegralSpaceSaturation(FlattenMatrices(gramsS)), n);
  for (auto &G : GB) {
    G /= s0;
  }
  int k = static_cast<int>(GB.size());
  MyMatrix<T> Q(k, k);
  for (int a = 0; a < k; a++) {
    for (int b = 0; b < k; b++) {
      MyMatrix<T> P = Ainv * GB[a] * Ainv * GB[b];
      T tr(0);
      for (int u = 0; u < n; u++) {
        tr += P(u, u);
      }
      Q(a, b) = tr;
    }
  }
  MyMatrix<T> Qinv = Inverse(Q);
  std::vector<int64_t> box(k);
  for (int a = 0; a < k; a++) {
    T r2 = bound2 * Qinv(a, a);
    int64_t r = 0;
    while (UniversalScalarConversion<T, int64_t>((r + 1) * (r + 1)) <= r2) {
      r++;
    }
    box[a] = r;
  }
  LatticeMethodChoice choice = DefaultPairMethods().canonical;
  // the canonical form of the form scaled to be integral
  auto key_of = [&](GroupLattice<T> const &L)
      -> std::pair<std::vector<MyMatrix<T>>, GroupLattice<T>> {
    T den(1);
    for (int i = 0; i < L.gram.rows(); i++) {
      for (int j = 0; j < L.gram.cols(); j++) {
        den = LCMpair(den, GetDenominator(L.gram(i, j)));
      }
    }
    std::vector<MyMatrix<T>> LP;
    for (auto const &g : L.gens) {
      LP.push_back(g.transpose());
    }
    MyMatrix<T> As = den * L.gram;
    MyMatrix<Tint> B = PairCanonicalBasis<T, Tint, Tgroup>(As, LP, choice, os);
    GroupLattice<T> Lc =
        GroupLatticeInBasis(L, UniversalMatrixConversion<T, Tint>(B));
    std::vector<MyMatrix<T>> key{Lc.gram};
    for (auto const &g : Lc.gens) {
      key.push_back(g);
    }
    return {key, Lc};
  };
  std::set<std::vector<MyMatrix<T>>> seen_tw;
  std::vector<GroupLattice<T>> twists;
  std::vector<int64_t> v(k);
  std::function<void(int)> rec = [&](int i) -> void {
    if (i == k) {
      MyMatrix<T> G = ZeroMatrix<T>(n, n);
      for (int u = 0; u < k; u++) {
        G += UniversalScalarConversion<T, int64_t>(v[u]) * GB[u];
      }
      MyMatrix<T> phi = Ainv * G;
      MyMatrix<T> phi2 = phi * phi;
      T tr(0), tr2(0);
      for (int u = 0; u < n; u++) {
        tr += phi(u, u);
        tr2 += phi2(u, u);
      }
      if (tr <= 0 || tr > Btr || tr2 > bound2) {
        return;
      }
      T detG = DeterminantMat(G);
      if (detG <= 0 || detG > det_max || !IsPositiveDefinite(G, std::cerr)) {
        return;
      }
      auto [key, Lc] = key_of(GroupLattice<T>{G, M.gens});
      if (seen_tw.insert(key).second) {
        twists.push_back(Lc);
      }
      return;
    }
    for (int64_t x = -box[i]; x <= box[i]; x++) {
      v[i] = x;
      rec(i + 1);
    }
  };
  rec(0);
  os << "ISOM_PIECES: twists, trace bound " << Btr << ", " << twists.size()
     << " twists\n";
  // the starts: the sums of c twists; at each prime of G up to a maximal
  // p^{-2a}-integral Lambda'_p-lattice
  std::set<std::vector<MyMatrix<T>>> known;
  std::vector<GroupLattice<T>> queue;
  auto insert = [&](GroupLattice<T> const &L) -> void {
    auto [key, Lc] = key_of(L);
    if (known.insert(key).second) {
      queue.push_back(Lc);
    }
  };
  // the maximal orders Lambda'_p, on M''
  std::vector<MaximalOrderData<T, Tint>> orders;
  for (auto const &p : gp) {
    orders.push_back(*MaximalOrderAt<T, Tint>(M, p, 100000000));
  }
  // up to the maximal lattices at the primes of G where they are unique;
  // at the others (bad primes, MaximalLattices.h) the pieces are reached
  // by all the steps up from the sums of twists, as at the primes outside
  // G
  std::vector<bool> unique_at;
  for (auto const &p : gp) {
    unique_at.push_back(MaximalLatticesUniqueAt<T, Tint>(M, p));
    if (!unique_at.back()) {
      os << "ISOM_PIECES: twists, the maximal lattices are not unique at p="
         << p << ", all the steps up there\n";
    }
  }
  auto up_to_maximal = [&](GroupLattice<T> L) -> GroupLattice<T> {
    for (size_t ip = 0; ip < gp.size(); ip++) {
      if (unique_at[ip]) {
        L = UpToMaximalLattice<T, Tint>(L, orders[ip], gp[ip], n, false);
      }
    }
    return L;
  };
  std::vector<size_t> idx(c, 0);
  std::function<void(int, size_t)> rec2 = [&](int i, size_t from) -> void {
    if (i == c) {
      GroupLattice<T> L = twists[idx[0]];
      for (int u = 1; u < c; u++) {
        L = DirectSumGroupLattice(L, twists[idx[u]]);
      }
      insert(up_to_maximal(L));
      return;
    }
    for (size_t t = from; t < twists.size(); t++) {
      idx[i] = t;
      rec2(i + 1, t);
    }
  };
  if (!twists.empty()) {
    rec2(0, 0);
  }
  size_t n_start = queue.size();
  // the integral overlattices at the primes outside G and at the primes of
  // G where the maximal lattices are not unique
  std::set<Tint> gp_set;
  for (size_t ip = 0; ip < gp.size(); ip++) {
    if (unique_at[ip]) {
      gp_set.insert(gp[ip]);
    }
  }
  for (size_t pos = 0; pos < queue.size(); pos++) {
    GroupLattice<T> L = queue[pos];
    T det = DeterminantMat(L.gram);
    Tint num = UniversalScalarConversion<Tint, T>(T_abs(det) * GetDenominator(det));
    if (num <= 1) {
      continue;
    }
    T den(1);
    for (int i = 0; i < L.gram.rows(); i++) {
      for (int j = 0; j < L.gram.cols(); j++) {
        den = LCMpair(den, GetDenominator(L.gram(i, j)));
      }
    }
    for (auto const &[q, e] : FactorsIntMap(num)) {
      if (e < 2 || gp_set.count(q) > 0) {
        continue;
      }
      GroupLattice<T> Ls = L;
      Ls.gram *= den;
      T_small_prime q_s = UniversalScalarConversion<T_small_prime, Tint>(q);
      size_t n_down = 0;
      for (auto N : FlipNeighbors(Ls, q_s, n_down, false, true)) {
        N.gram /= den;
        insert(N);
      }
    }
  }
  os << "ISOM_PIECES: twists, " << n_start << " starts, " << queue.size()
     << " lattices before the steps down\n";
  return StepsDownWithCut<T, Tint, Tgroup>(queue, dets, gp, hi, os);
}

// clang-format off
#endif  // SRC_GENERA_ISOTYPICPIECES_H_
// clang-format on
