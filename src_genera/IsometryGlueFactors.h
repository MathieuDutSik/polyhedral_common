// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ISOMETRYGLUEFACTORS_H_
#define SRC_GENERA_ISOMETRYGLUEFACTORS_H_

// clang-format off
#include "IsometryGlue.h"
#include "GenusSymbol.h"
#include <algorithm>
#include <fstream>
#include <map>
#include <optional>
#include <sstream>
#include <set>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  The seeds of every genus of pairs (L, sigma), L in the genus of a given
  lattice and sigma of characteristic polynomial prod_i mu_i^{b_i}, the
  mu_i pairwise coprime (cyclotomic for sigma of finite order).

  The pieces P_i = L cap ker mu_i(sigma) are primitive and sigma-stable,
  and so are M_j = L cap (V_1 + ... + V_j): M_j is a sigma-stable primitive
  extension of M_{j-1} (+) P_j, so L is reached by gluing one piece at a
  time (IsometryGlue.h), and gluing a representative of every genus of
  pairs at every step gives a seed in every genus of pairs of L.

  The bound. Let X be a sigma-stable primitive sublattice of L, Y its
  orthogonal complement, sigma acting on them with coprime minimal
  polynomials f and g. The glue H = L / (X (+) Y) is the graph of a
  sigma-equivariant isomorphism between subgroups of D(X) and D(Y): it is
  killed by f(sigma) and by g(sigma), so by Res(f, g) = u f + v g (u, v in
  Z[x]), and it has at most min(rank X, rank Y) generators, so its order s
  divides S = Res(f, g)^min(rank X, rank Y). Better, when X is a piece
  P_i with f = Phi_d: P_i is a module of rank b_i = rank P_i / phi(d) over
  O = Z[zeta_d], D(P_i) a sum of b_i cyclic O-modules O / I, and H, killed
  by g(zeta_d) in O, lies in the sum of the (O / I)[g(zeta_d)], each of
  order dividing N(g(zeta_d)) = |Res(f, g)|: s divides Res(f, g)^b_i (b_i
  = rank P_i for f = x -+ 1). So the exponent is min(c_X, c_Y), c = b_i
  for a single piece, the rank otherwise. As H embeds in D(X) and in
  D(Y) and |D(X)| |D(Y)| = d s^2: s divides det X and det X divides d s,
  that is det X = s e with s | S and e | d.

  Polynomials are the vectors of their coefficients, of increasing degree
  (PolyProductCoeffs, MinimalPolynomialCoeffs and PolyResultant are in
  IsometryGlue.h).
 */

// Phi_d, by dividing x^d - 1 by the Phi_e, e | d, e < d.
template <typename T> std::vector<T> CyclotomicPolynomialCoeffs(int d) {
  std::vector<T> num(d + 1, T(0));
  num[0] = T(-1);
  num[d] = T(1);
  for (int e = 1; e < d; e++) {
    if (d % e != 0) {
      continue;
    }
    std::vector<T> den = CyclotomicPolynomialCoeffs<T>(e);
    // exact division by the monic den
    int dn = static_cast<int>(num.size()) - 1;
    int dd = static_cast<int>(den.size()) - 1;
    std::vector<T> quo(dn - dd + 1, T(0));
    for (int k = dn - dd; k >= 0; k--) {
      quo[k] = num[k + dd];
      for (int i = 0; i <= dd; i++) {
        num[k + i] -= quo[k] * den[i];
      }
    }
    num = quo;
  }
  return num;
}

// The positive divisors of N > 0.
template <typename Tint> std::vector<Tint> DivisorsOf(Tint const &N) {
  std::vector<Tint> divisors{Tint(1)};
  if (N > 1) {
    for (auto const &[q, mult] : FactorsIntMap(N)) {
      std::vector<Tint> next;
      for (auto const &x : divisors) {
        Tint y = x;
        for (size_t i = 0; i <= mult; i++) {
          next.push_back(y);
          y *= q;
        }
      }
      divisors = next;
    }
  }
  std::sort(divisors.begin(), divisors.end());
  return divisors;
}

/*
  The possible determinants s e (s | S, e | d) of a sigma-stable primitive
  sublattice X of L, det L = d, with S = Res^min(rank X, rank Y) as above.
 */
template <typename Tint>
std::vector<Tint> SublatticeDeterminants(Tint const &d, Tint const &res,
                                         int max_gen) {
  Tint S = MyPow<Tint>(T_abs(res), max_gen);
  std::set<Tint> dets;
  for (auto const &s : DivisorsOf(S)) {
    for (auto const &e : DivisorsOf(d)) {
      dets.insert(s * e);
    }
  }
  return std::vector<Tint>(dets.begin(), dets.end());
}

/*
  One factor: the pairs (P, sigma_P) of a list, one per genus of pairs, all
  of minimal polynomial mu and of rank rank.
  - orth: sigma_P = sign Id, the pieces are the even genera of rank r of
    the admissible determinants, made here;
  - file: the pieces of a file (Gram matrix, then isometry, in turn), e.g.
    the hermitian genera of HermitianPieces.jl.
 */
template <typename T> struct IsometryFactor {
  std::string desc;
  bool orth;
  int sign;
  int rank;
  std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> pieces;
  std::vector<T> mu;
  // the rank over Z[x] / mu: rank / deg mu
  int b;
};

/*
  A factor from one line of a spec file:
      orth 1 r
      orth -1 r
      file <path>
 */
template <typename T>
IsometryFactor<T> ParseIsometryFactor(std::string const &line) {
  std::istringstream ls(line);
  std::string kind;
  ls >> kind;
  IsometryFactor<T> f;
  f.desc = line;
  if (kind == "orth") {
    f.orth = true;
    if (!(ls >> f.sign >> f.rank) || (f.sign != 1 && f.sign != -1) ||
        f.rank <= 0) {
      std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
      throw TerminalException{1};
    }
    f.mu = {T(-f.sign), T(1)};
    f.b = f.rank;
    return f;
  }
  if (kind != "file") {
    std::cerr << "ISOM_FACTORS: unknown kind " << kind << "\n";
    throw TerminalException{1};
  }
  std::string path;
  if (!(ls >> path)) {
    std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
    throw TerminalException{1};
  }
  f.orth = false;
  f.sign = 0;
  std::vector<MyMatrix<T>> lst = ReadListMatrixFile<T>(path);
  if (lst.size() % 2 != 0) {
    std::cerr << "ISOM_FACTORS: " << path
              << " has an odd number of matrices\n";
    throw TerminalException{1};
  }
  for (size_t i = 0; i + 1 < lst.size(); i += 2) {
    f.pieces.push_back({lst[i], lst[i + 1]});
  }
  if (f.pieces.empty()) {
    // no piece of that kind: no pair at all
    f.rank = -1;
    return f;
  }
  f.rank = f.pieces[0].first.rows();
  f.mu = MinimalPolynomialCoeffs(f.pieces[0].second);
  f.b = f.rank / (static_cast<int>(f.mu.size()) - 1);
  for (auto const &[A, g] : f.pieces) {
    CheckIsometryInput(A, g);
    if (A.rows() != f.rank || MinimalPolynomialCoeffs(g) != f.mu) {
      std::cerr << "ISOM_FACTORS: the pieces of " << path
                << " differ in rank or minimal polynomial\n";
      throw TerminalException{1};
    }
  }
  return f;
}

// The factors of a spec file, one per line, # for the comments.
template <typename T>
std::vector<IsometryFactor<T>> ReadIsometryFactors(std::string const &file) {
  std::ifstream is(file);
  if (!is) {
    std::cerr << "ISOM_FACTORS: cannot open " << file << "\n";
    throw TerminalException{1};
  }
  std::vector<IsometryFactor<T>> res;
  std::string line;
  while (std::getline(is, line)) {
    std::istringstream ls(line);
    std::string kind;
    if (!(ls >> kind) || kind[0] == '#') {
      continue;
    }
    res.push_back(ParseIsometryFactor<T>(line));
  }
  return res;
}

// The product of the minimal polynomials of the factors in idx.
template <typename T>
std::vector<T> FactorsProduct(std::vector<IsometryFactor<T>> const &factors,
                              std::vector<size_t> const &idx) {
  std::vector<T> p{T(1)};
  for (size_t i : idx) {
    p = PolyProductCoeffs(p, factors[i].mu);
  }
  return p;
}

/*
  The bound on the number of generators of the glue, over Z[x] / mu for a
  single factor: b, the rank otherwise.
 */
template <typename T>
int GlueGenerators(std::vector<IsometryFactor<T>> const &factors,
                   std::vector<size_t> const &idx) {
  if (idx.size() == 1) {
    return factors[idx[0]].b;
  }
  int rank = 0;
  for (size_t i : idx) {
    rank += factors[i].rank;
  }
  return rank;
}

/*
  The possible determinants of the sublattice of L on the factors in idx,
  the others making its complement.
 */
template <typename T, typename Tint>
std::vector<Tint>
FactorDeterminants(std::vector<IsometryFactor<T>> const &factors,
                   std::vector<size_t> const &idx, Tint const &d) {
  std::vector<size_t> rest;
  for (size_t i = 0; i < factors.size(); i++) {
    if (std::find(idx.begin(), idx.end(), i) == idx.end()) {
      rest.push_back(i);
    }
  }
  T res = PolyResultant(FactorsProduct(factors, idx),
                        FactorsProduct(factors, rest));
  if (res == 0) {
    std::cerr << "ISOM_FACTORS: the minimal polynomials are not coprime\n";
    throw TerminalException{1};
  }
  return SublatticeDeterminants(
      d, UniversalScalarConversion<Tint, T>(res),
      std::min(GlueGenerators(factors, idx), GlueGenerators(factors, rest)));
}

/*
  Whether X (+) Y, of determinants det1 and det2, can be glued into a
  lattice of determinant det by a glue whose order divides S: s^2 = det1
  det2 / det, s | S, s | det1 and s | det2.
 */
template <typename T, typename Tint>
bool IsAdmissibleGlue(T const &det1, T const &det2, T const &det,
                      Tint const &S) {
  T sq = det1 * det2 / det;
  if (!IsInteger(sq)) {
    return false;
  }
  Tint sq_i = UniversalScalarConversion<Tint, T>(sq);
  Tint s = IntegerSquareRoot(sq_i);
  if (s * s != sq_i) {
    return false;
  }
  Tint d1 = UniversalScalarConversion<Tint, T>(det1);
  Tint d2 = UniversalScalarConversion<Tint, T>(det2);
  return ResInt(S, s) == 0 && ResInt(d1, s) == 0 && ResInt(d2, s) == 0;
}

template <typename T, typename Tint>
bool IsAdmissibleDeterminant(std::vector<Tint> const &dets, T const &det) {
  Tint det_i = UniversalScalarConversion<Tint, T>(det);
  return std::binary_search(dets.begin(), dets.end(), det_i);
}

/*
  The seeds: one glued pair at least in every genus of pairs (L, sigma), L
  in the genus of A, sigma acting on the factors as given.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<GluedPair<T>>
FactorSeeds(MyMatrix<T> const &A, std::vector<IsometryFactor<T>> factors,
            std::ostream &os) {
  int n = A.rows();
  T det = DeterminantMat(A);
  Tint d = UniversalScalarConversion<Tint, T>(det);
  GenusSymbol S = GenusSymbolOfGram<T, Tint>(A);
  int sum_rank = 0;
  for (auto const &f : factors) {
    if (f.rank < 0) {
      std::cerr << "ISOM_FACTORS: " << f.desc << ": no piece\n";
      return {};
    }
    sum_rank += f.rank;
  }
  if (factors.empty() || sum_rank != n) {
    std::cerr << "ISOM_FACTORS: the ranks of the factors add up to "
              << sum_rank << ", not " << n << "\n";
    throw TerminalException{1};
  }
  // the pieces of each factor, of an admissible determinant
  for (size_t i = 0; i < factors.size(); i++) {
    IsometryFactor<T> &f = factors[i];
    std::vector<Tint> dets = FactorDeterminants(factors, {i}, d);
    if (f.orth) {
      MyMatrix<T> g = T(f.sign) * IdentityMat<T>(f.rank);
      for (auto const &D : dets) {
        for (auto const &gs : EnumerateEvenGenusSymbols<Tint>(f.rank, D)) {
          f.pieces.push_back({GenusSymbolRepresentative<T, Tint>(gs, os), g});
        }
      }
    } else {
      std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> kept;
      for (auto const &pc : f.pieces) {
        if (IsAdmissibleDeterminant(dets, DeterminantMat(pc.first))) {
          kept.push_back(pc);
        }
      }
      f.pieces = kept;
    }
    std::cerr << "ISOM_FACTORS: " << f.desc << ": " << f.pieces.size()
              << " pieces of admissible determinant\n";
    if (f.pieces.empty()) {
      return {};
    }
  }
  // glue the factors with the fewest pieces first
  std::sort(factors.begin(), factors.end(),
            [](IsometryFactor<T> const &x, IsometryFactor<T> const &y) {
              return x.pieces.size() < y.pieces.size();
            });
  IsometryGlueBounds bounds;
  // the glued pairs so far, glue_order the product of the glue orders
  std::vector<GluedPair<T>> current;
  for (auto const &[A1, g1] : factors[0].pieces) {
    current.push_back({A1, g1, DeterminantMat(A1), T(1)});
  }
  std::vector<size_t> idx{0};
  for (size_t j = 1; j < factors.size(); j++) {
    idx.push_back(j);
    bool last = (j + 1 == factors.size());
    std::vector<Tint> dets;
    if (!last) {
      dets = FactorDeterminants(factors, idx, d);
    }
    std::optional<T> target;
    if (last) {
      target = det;
    }
    // the glue between M_{j-1} and P_j: killed by the resultant of their
    // minimal polynomials, its order dividing Res^min(c_prev, b_j)
    std::vector<size_t> idx_prev(idx.begin(), idx.end() - 1);
    T res = PolyResultant(FactorsProduct(factors, idx_prev), factors[j].mu);
    Tint S_glue =
        MyPow<Tint>(UniversalScalarConversion<Tint, T>(T_abs(res)),
                    std::min(GlueGenerators(factors, idx_prev),
                             GlueGenerators(factors, {j})));
    std::set<std::vector<MyMatrix<T>>> known;
    std::vector<GluedPair<T>> next;
    for (auto const &cur : current) {
      for (auto const &[A2, g2] : factors[j].pieces) {
        T det2 = DeterminantMat(A2);
        if (last && !IsAdmissibleGlue(cur.det, det2, det, S_glue)) {
          continue;
        }
        if (!last) {
          // some admissible determinant of M_j must be reachable
          bool reachable = false;
          for (auto const &D : dets) {
            T D_T = UniversalScalarConversion<T, Tint>(D);
            if (IsAdmissibleGlue(cur.det, det2, D_T, S_glue)) {
              reachable = true;
              break;
            }
          }
          if (!reachable) {
            continue;
          }
        }
        for (auto const &gl : GlueIsometryPairs<T, Tint, Tgroup>(
                 cur.gram, cur.isom, A2, g2, target, bounds, os)) {
          if (!last && !IsAdmissibleDeterminant(dets, gl.det)) {
            continue;
          }
          if (last && GenusSymbolOfGram<T, Tint>(gl.gram) != S) {
            continue;
          }
          // GlueIsometryPairs returns canonical pairs
          if (known.insert(PairForms(gl.gram, std::vector<MyMatrix<T>>{
                                                  gl.isom.transpose()}))
                  .second) {
            next.push_back(
                {gl.gram, gl.isom, gl.det, cur.glue_order * gl.glue_order});
          }
        }
      }
    }
    std::cerr << "ISOM_FACTORS: step " << j << ", " << next.size()
              << " glued pairs\n";
    current = next;
  }
  std::vector<GluedPair<T>> res;
  for (auto const &gp : current) {
    if (factors.size() == 1 &&
        (gp.det != det || GenusSymbolOfGram<T, Tint>(gp.gram) != S)) {
      continue;
    }
    res.push_back(gp);
  }
  std::cerr << "ISOM_FACTORS: " << res.size()
            << " glued pairs in the genus of L\n";
  return res;
}

/*
  Odd lattices. The even sublattice L_0 = {x : x.x even} of an odd lattice L
  has index 2 (x -> x.x mod 2 is linear) and is stable under every
  isometry of L. A pair (L, sigma) is the triple (L_0, sigma, w), with L =
  L_0 + w and w in D(L_0) of order 2, fixed by sigma, of odd norm w.w; the
  even sublattice of L_0 + w is L_0 back, and the centralizer of sigma in
  O(L) is the stabilizer of w in the centralizer C of sigma in O(L_0). So
  the classes of the pairs (L, sigma) with L of a given genus are, for the
  classes of the pairs (L_0, sigma) with L_0 in the genus of the even
  sublattices, the C-orbits of these w with L_0 + w in that genus, of mass
  |orbit| / |C|.
 */

// A basis (rows, in the coordinates of A) of the even sublattice of the
// odd lattice of Gram matrix A.
template <typename T> MyMatrix<T> EvenSublatticeBasis(MyMatrix<T> const &A) {
  int n = A.rows();
  std::vector<int> odd;
  for (int i = 0; i < n; i++) {
    if (!IsInteger(A(i, i) / T(2))) {
      odd.push_back(i);
    }
  }
  if (odd.empty()) {
    std::cerr << "ISOM_FACTORS: the lattice is even\n";
    throw TerminalException{1};
  }
  int i0 = odd[0];
  MyMatrix<T> B = IdentityMat<T>(n);
  for (int j : odd) {
    if (j == i0) {
      B(j, j) = T(2);
    } else {
      B(j, i0) = T(1);
    }
  }
  return B;
}

template <typename T> struct OddOverlattice {
  // a basis of L_0 + w, rows in the coordinates of L_0
  MyMatrix<T> basis;
  size_t orbit_size;
};

/*
  The orbits of the group generated by ListGen (isometries of A0 acting on
  the rows, h A0 h^T = A0, commuting with P) on the w = v / 2 of D(A0), v
  in {0, 1}^n, with v A0 = 0 mod 2, v P = v mod 2 and v A0 v^T = 4 mod 8
  (w.w odd), one overlattice A0 + w per orbit with the orbit size.
 */
template <typename T, typename Tint>
std::vector<OddOverlattice<T>>
OddOverlatticeOrbits(MyMatrix<T> const &A0, MyMatrix<T> const &P,
                     std::vector<MyMatrix<Tint>> const &ListGen) {
  int n = A0.rows();
  auto mod2 = [](T const &x) -> int {
    Tint xi = UniversalScalarConversion<Tint, T>(x);
    return (ResInt(xi, Tint(2)) == 0) ? 0 : 1;
  };
  // the kernel of A0 mod 2, by Gaussian elimination over F_2
  std::vector<std::vector<int>> M(n, std::vector<int>(n));
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      M[i][j] = mod2(A0(i, j));
    }
  }
  std::vector<int> pivot_col;
  int rank = 0;
  for (int c = 0; c < n && rank < n; c++) {
    int piv = -1;
    for (int r = rank; r < n; r++) {
      if (M[r][c]) {
        piv = r;
        break;
      }
    }
    if (piv < 0) {
      continue;
    }
    std::swap(M[piv], M[rank]);
    for (int r = 0; r < n; r++) {
      if (r != rank && M[r][c]) {
        for (int k = 0; k < n; k++) {
          M[r][k] ^= M[rank][k];
        }
      }
    }
    pivot_col.push_back(c);
    rank++;
  }
  // A0 is symmetric: the kernel of v -> v A0 is that of x -> A0 x
  std::vector<std::vector<int>> kernel;
  for (int f = 0; f < n; f++) {
    if (std::find(pivot_col.begin(), pivot_col.end(), f) != pivot_col.end()) {
      continue;
    }
    std::vector<int> v(n, 0);
    v[f] = 1;
    for (int r = 0; r < rank; r++) {
      v[pivot_col[r]] = M[r][f];
    }
    kernel.push_back(v);
  }
  if (kernel.size() > 24) {
    std::cerr << "ISOM_FACTORS: D(L_0)[2] is too large\n";
    throw TerminalException{1};
  }
  auto key = [](std::vector<int> const &v) -> uint64_t {
    uint64_t k = 0;
    for (size_t i = 0; i < v.size(); i++) {
      k |= static_cast<uint64_t>(v[i]) << i;
    }
    return k;
  };
  auto act = [&](std::vector<int> const &v,
                 MyMatrix<T> const &h) -> std::vector<int> {
    std::vector<int> w(n, 0);
    for (int j = 0; j < n; j++) {
      int s = 0;
      for (int i = 0; i < n; i++) {
        if (v[i]) {
          s ^= mod2(h(i, j));
        }
      }
      w[j] = s;
    }
    return w;
  };
  std::vector<MyMatrix<T>> gens;
  for (auto const &h : ListGen) {
    gens.push_back(UniversalMatrixConversion<T, Tint>(h));
  }
  // the admissible w
  std::set<uint64_t> admissible;
  std::vector<std::vector<int>> list;
  size_t n_elt = size_t(1) << kernel.size();
  for (size_t u = 1; u < n_elt; u++) {
    std::vector<int> v(n, 0);
    for (size_t b = 0; b < kernel.size(); b++) {
      if ((u >> b) & 1) {
        for (int i = 0; i < n; i++) {
          v[i] ^= kernel[b][i];
        }
      }
    }
    if (act(v, P) != v) {
      continue;
    }
    T norm(0);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        if (v[i] && v[j]) {
          norm += A0(i, j);
        }
      }
    }
    Tint norm_i = UniversalScalarConversion<Tint, T>(norm);
    if (ResInt(norm_i, Tint(8)) != 4) {
      continue;
    }
    admissible.insert(key(v));
    list.push_back(v);
  }
  std::vector<OddOverlattice<T>> res;
  std::set<uint64_t> seen;
  for (auto const &v0 : list) {
    if (seen.count(key(v0)) == 1) {
      continue;
    }
    std::vector<std::vector<int>> queue{v0};
    seen.insert(key(v0));
    size_t orbit_size = 0;
    while (!queue.empty()) {
      std::vector<int> v = queue.back();
      queue.pop_back();
      orbit_size++;
      for (auto const &h : gens) {
        std::vector<int> w = act(v, h);
        if (admissible.count(key(w)) == 0) {
          std::cerr << "ISOM_FACTORS: the centralizer does not preserve the "
                       "admissible classes\n";
          throw TerminalException{1};
        }
        if (seen.insert(key(w)).second) {
          queue.push_back(w);
        }
      }
    }
    // A0 + v / 2, a Z-basis from the integer generators 2 Id and v
    MyMatrix<Tint> Gen_i = ZeroMatrix<Tint>(n + 1, n);
    for (int i = 0; i < n; i++) {
      Gen_i(i, i) = Tint(2);
      Gen_i(n, i) = Tint(v0[i]);
    }
    MyMatrix<T> B =
        UniversalMatrixConversion<T, Tint>(GetZbasis(Gen_i)) / T(2);
    res.push_back({B, orbit_size});
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_ISOMETRYGLUEFACTORS_H_
// clang-format on
