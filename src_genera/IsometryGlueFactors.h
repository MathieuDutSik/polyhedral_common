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
  One factor: the pairs (P, G_P) of a list, one per genus of pairs, G acting
  through its generators g_1, ..., g_k (k = 1 for a single isometry), all of
  rank rank and with the same minimal polynomial mus[j] for g_j.
  - orth: g_j = signs[j] Id, the pieces are the even genera of rank r of the
    admissible determinants, made here;
  - file: the pieces of a file (Gram matrix, then isometry g, in turn),
    e.g. the hermitian genera of GENUS_HermitianPieces, with g_j = g^e_j.
  For an abelian G the pieces of L are its isotypic parts, on which G acts
  through a character: the factors must have distinct characters.
 */
template <typename T> struct IsometryFactor {
  std::string desc;
  bool orth;
  std::vector<int> signs;
  std::vector<int64_t> powers;
  int rank;
  // the Gram matrix and the actions of the generators (on the columns)
  std::vector<std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>> pieces;
  std::vector<std::vector<T>> mus;
};

template <typename T>
MyMatrix<T> FactorMatrixPower(MyMatrix<T> const &g, int64_t e) {
  MyMatrix<T> res = IdentityMat<T>(g.rows());
  for (int64_t i = 0; i < e; i++) {
    res = res * g;
  }
  return res;
}

/*
  A factor from one line of a spec file, for a group of n_gen generators:
      orth s_1 ... s_k r
      file <path> e_1 ... e_k   (e_1 = 1 by default when k = 1)
 */
template <typename T>
IsometryFactor<T> ParseIsometryFactor(std::string const &line, int n_gen) {
  std::istringstream ls(line);
  std::string kind;
  ls >> kind;
  IsometryFactor<T> f;
  f.desc = line;
  if (kind == "orth") {
    f.orth = true;
    for (int j = 0; j < n_gen; j++) {
      int s = 0;
      if (!(ls >> s) || (s != 1 && s != -1)) {
        std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
        throw TerminalException{1};
      }
      f.signs.push_back(s);
      f.mus.push_back({T(-s), T(1)});
    }
    if (!(ls >> f.rank) || f.rank <= 0) {
      std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
      throw TerminalException{1};
    }
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
  for (int j = 0; j < n_gen; j++) {
    int64_t e = 1;
    if (!(ls >> e)) {
      if (n_gen > 1) {
        std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
        throw TerminalException{1};
      }
      e = 1;
    }
    if (e < 0) {
      std::cerr << "ISOM_FACTORS: the powers must be non negative\n";
      throw TerminalException{1};
    }
    f.powers.push_back(e);
  }
  std::vector<MyMatrix<T>> lst = ReadListMatrixFile<T>(path);
  if (lst.size() % 2 != 0) {
    std::cerr << "ISOM_FACTORS: " << path
              << " has an odd number of matrices\n";
    throw TerminalException{1};
  }
  for (size_t i = 0; i + 1 < lst.size(); i += 2) {
    CheckIsometryInput(lst[i], lst[i + 1]);
    std::vector<MyMatrix<T>> gens;
    for (int64_t e : f.powers) {
      gens.push_back(FactorMatrixPower(lst[i + 1], e));
    }
    f.pieces.push_back({lst[i], gens});
  }
  if (f.pieces.empty()) {
    // no piece of that kind: no pair at all
    f.rank = -1;
    return f;
  }
  f.rank = f.pieces[0].first.rows();
  for (auto const &g : f.pieces[0].second) {
    f.mus.push_back(MinimalPolynomialCoeffs(g));
  }
  for (auto const &[A, gens] : f.pieces) {
    bool same = A.rows() == f.rank;
    for (size_t j = 0; same && j < gens.size(); j++) {
      same = MinimalPolynomialCoeffs(gens[j]) == f.mus[j];
    }
    if (!same) {
      std::cerr << "ISOM_FACTORS: the pieces of " << path
                << " differ in rank or minimal polynomial\n";
      throw TerminalException{1};
    }
  }
  return f;
}

// One generator unless stated by a first line "group k".
template <typename T>
IsometryFactor<T> ParseIsometryFactor(std::string const &line) {
  return ParseIsometryFactor<T>(line, 1);
}

// The factors of a spec file, one per line, # for the comments, a first
// line "group k" for a group of k generators.
template <typename T>
std::vector<IsometryFactor<T>> ReadIsometryFactors(std::string const &file) {
  std::ifstream is(file);
  if (!is) {
    std::cerr << "ISOM_FACTORS: cannot open " << file << "\n";
    throw TerminalException{1};
  }
  std::vector<IsometryFactor<T>> res;
  std::string line;
  int n_gen = 1;
  while (std::getline(is, line)) {
    std::istringstream ls(line);
    std::string kind;
    if (!(ls >> kind) || kind[0] == '#') {
      continue;
    }
    if (kind == "group") {
      if (!res.empty() || !(ls >> n_gen) || n_gen < 1) {
        std::cerr << "ISOM_FACTORS: bad group line: " << line << "\n";
        throw TerminalException{1};
      }
      continue;
    }
    res.push_back(ParseIsometryFactor<T>(line, n_gen));
  }
  return res;
}

/*
  The minimal polynomial of g_j on the sum of the factors in idx: the lcm
  of theirs, the product of the distinct ones (cyclotomic polynomials).
 */
template <typename T>
std::vector<T> FactorsProduct(std::vector<IsometryFactor<T>> const &factors,
                              std::vector<size_t> const &idx, size_t j) {
  std::vector<T> p{T(1)};
  std::set<std::vector<T>> seen;
  for (size_t i : idx) {
    if (seen.insert(factors[i].mus[j]).second) {
      p = PolyProductCoeffs(p, factors[i].mus[j]);
    }
  }
  return p;
}

/*
  The bound on the number of generators of the glue for g_j, over Z[x] /
  mu_j for a single factor (rank / deg mu_j), the rank otherwise.
 */
template <typename T>
int GlueGenerators(std::vector<IsometryFactor<T>> const &factors,
                   std::vector<size_t> const &idx, size_t j) {
  if (idx.size() == 1) {
    IsometryFactor<T> const &f = factors[idx[0]];
    return f.rank / (static_cast<int>(f.mus[j].size()) - 1);
  }
  int rank = 0;
  for (size_t i : idx) {
    rank += factors[i].rank;
  }
  return rank;
}

// n with Phi_n = mu (mu cyclotomic), by trying the n up to a bound.
template <typename T> int64_t CyclotomicIndexOf(std::vector<T> const &mu) {
  int64_t deg = static_cast<int64_t>(mu.size()) - 1;
  for (int64_t n = 1; n <= 4 * deg * deg + 6; n++) {
    if (CyclotomicPolynomialCoeffs<T>(static_cast<int>(n)) == mu) {
      return n;
    }
  }
  std::cerr << "ISOM_FACTORS: a minimal polynomial is not cyclotomic\n";
  throw TerminalException{1};
}

/*
  Bounds on the glue H between the sums X and Y of the factors in idx_in
  and in idx_out: an integer killing it (exponent, 0 when none is known)
  and a multiple of its order.

  One generator: H is killed by the resultant of the minimal polynomials
  and has at most min(c_X, c_Y) generators (see above).

  A group (abelian, the factors its isotypic parts), X a single factor of
  character chi_X. A factor mu(g_j), mu the minimal polynomial of g_j on a
  character chi of Y, kills the characters psi of Y with the same minimal
  polynomial for g_j (mu_{psi, j} = mu); a set of such factors covering the
  characters of Y gives a product a killing Y, so H, and acting on X as an
  algebraic integer alpha of Q(chi_X) whose norm E = prod |Res(mu_{X,j},
  mu)|^[Q(chi_X) : Q(zeta_{X,j})] kills H (the cover of least E is taken,
  by a search over the subsets). H, a module over Z[chi_X] with at most c_X
  = rank X / [Q(chi_X) : Q] generators, has order dividing E^min(c_X, c_Y)
  (c_Y likewise when Y is a single factor, rank Y otherwise).

  X and Y both sums of several factors: the glue of L over the sum of all
  the pieces P_i embeds in the sum of its projections, the glues H_i of the
  P_i with their complements, of orders dividing the bounds S_i above; H is
  a quotient of it, so its order divides prod_i S_i. And by the resultant
  bound of any generator separating X and Y (coprime minimal polynomials on
  them), as for a single isometry, which bounds the exponent too.
 */
template <typename Tint> struct GlueBound {
  Tint exponent;
  Tint order;
};

// phi of the lcm of the orders of the eigenvalues: [Q(chi) : Q].
template <typename T>
int64_t CharacterFieldDegree(IsometryFactor<T> const &f) {
  int64_t ord = 1;
  for (auto const &mu : f.mus) {
    int64_t o = CyclotomicIndexOf(mu);
    ord = ord / std::gcd(ord, o) * o;
  }
  int64_t deg = 0;
  for (int64_t k = 1; k <= ord; k++) {
    if (std::gcd(k, ord) == 1) {
      deg++;
    }
  }
  return deg;
}

template <typename T, typename Tint>
GlueBound<Tint> GlueBounds(std::vector<IsometryFactor<T>> const &factors,
                           std::vector<size_t> const &idx_in,
                           std::vector<size_t> const &idx_out) {
  size_t n_gen = factors[0].mus.size();
  // no glue with an empty side
  if (idx_in.empty() || idx_out.empty()) {
    return {Tint(1), Tint(1)};
  }
  if (n_gen == 1) {
    T res = T_abs(PolyResultant(FactorsProduct(factors, idx_in, 0),
                                FactorsProduct(factors, idx_out, 0)));
    if (res == 0) {
      std::cerr << "ISOM_FACTORS: the minimal polynomials are not coprime\n";
      throw TerminalException{1};
    }
    Tint res_i = UniversalScalarConversion<Tint, T>(res);
    return {res_i, MyPow<Tint>(res_i, std::min(GlueGenerators(factors, idx_in, 0),
                                               GlueGenerators(factors, idx_out, 0)))};
  }
  int rank_in = 0, rank_out = 0;
  for (size_t i : idx_in) {
    rank_in += factors[i].rank;
  }
  for (size_t i : idx_out) {
    rank_out += factors[i].rank;
  }
  if (idx_in.size() != 1 && idx_out.size() != 1) {
    Tint prod(1);
    for (size_t i = 0; i < factors.size(); i++) {
      std::vector<size_t> rest;
      for (size_t k = 0; k < factors.size(); k++) {
        if (k != i) {
          rest.push_back(k);
        }
      }
      prod *= GlueBounds<T, Tint>(factors, {i}, rest).order;
    }
    Tint expo(0);
    for (size_t j = 0; j < n_gen; j++) {
      T res = T_abs(PolyResultant(FactorsProduct(factors, idx_in, j),
                                  FactorsProduct(factors, idx_out, j)));
      if (res != 0) {
        Tint res_i = UniversalScalarConversion<Tint, T>(res);
        expo = (expo == 0) ? res_i : GcdPair(expo, res_i);
        prod = GcdPair(prod, MyPow<Tint>(res_i, std::min(rank_in, rank_out)));
      }
    }
    if (expo == 0) {
      expo = prod;
    }
    return {expo, prod};
  }
  std::vector<size_t> const &iX = (idx_in.size() == 1) ? idx_in : idx_out;
  std::vector<size_t> const &iY = (idx_in.size() == 1) ? idx_out : idx_in;
  IsometryFactor<T> const &X = factors[iX[0]];
  int64_t deg_X = CharacterFieldDegree(X);
  // the candidate factors (j, mu): their cost and the characters of Y
  // they kill
  struct Cand {
    Tint cost;
    std::vector<size_t> kills;
  };
  std::vector<Cand> cands;
  for (size_t j = 0; j < n_gen; j++) {
    std::set<std::vector<T>> done;
    for (size_t i : iY) {
      std::vector<T> const &mu = factors[i].mus[j];
      if (mu == X.mus[j] || !done.insert(mu).second) {
        continue;
      }
      T res = T_abs(PolyResultant(X.mus[j], mu));
      int64_t e = deg_X / (static_cast<int64_t>(X.mus[j].size()) - 1);
      Cand c{MyPow<Tint>(UniversalScalarConversion<Tint, T>(res), e), {}};
      for (size_t k = 0; k < iY.size(); k++) {
        if (factors[iY[k]].mus[j] == mu) {
          c.kills.push_back(k);
        }
      }
      cands.push_back(c);
    }
  }
  if (cands.size() > 20) {
    std::cerr << "ISOM_FACTORS: too many candidate factors for the glue "
                 "bound\n";
    throw TerminalException{1};
  }
  std::optional<Tint> best;
  for (size_t mask = 1; mask < (size_t(1) << cands.size()); mask++) {
    std::vector<uint8_t> covered(iY.size(), 0);
    Tint cost(1);
    for (size_t c = 0; c < cands.size(); c++) {
      if ((mask >> c) & 1) {
        cost *= cands[c].cost;
        for (size_t k : cands[c].kills) {
          covered[k] = 1;
        }
      }
    }
    bool all = true;
    for (auto x : covered) {
      all = all && x;
    }
    if (all && (!best || cost < *best)) {
      best = cost;
    }
  }
  if (!best) {
    std::cerr << "ISOM_FACTORS: two factors have the same character\n";
    throw TerminalException{1};
  }
  Tint E = *best;
  int c_X = X.rank / static_cast<int>(deg_X);
  int c_Y = rank_in + rank_out - X.rank;
  if (iY.size() == 1) {
    c_Y = factors[iY[0]].rank /
          static_cast<int>(CharacterFieldDegree(factors[iY[0]]));
  }
  return {E, MyPow<Tint>(E, std::min(c_X, c_Y))};
}

// The possible determinants s e, s | S, e | d (see above).
template <typename Tint>
std::vector<Tint> SublatticeDeterminantsOfBound(Tint const &d, Tint const &S) {
  std::set<Tint> dets;
  for (auto const &s : DivisorsOf(S)) {
    for (auto const &e : DivisorsOf(d)) {
      dets.insert(s * e);
    }
  }
  return std::vector<Tint>(dets.begin(), dets.end());
}

/*
  The possible determinants of the sublattice of L on the factors in idx,
  the others making its complement.
 */
template <typename T, typename Tint>
std::optional<std::vector<Tint>>
FactorDeterminants(std::vector<IsometryFactor<T>> const &factors,
                   std::vector<size_t> const &idx, Tint const &d) {
  std::vector<size_t> rest;
  for (size_t i = 0; i < factors.size(); i++) {
    if (std::find(idx.begin(), idx.end(), i) == idx.end()) {
      rest.push_back(i);
    }
  }
  return SublatticeDeterminantsOfBound(
      d, GlueBounds<T, Tint>(factors, idx, rest).order);
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

// All the generators of a glued pair, acting on the columns.
template <typename T> std::vector<MyMatrix<T>> GluedGens(GluedPair<T> const &gp) {
  std::vector<MyMatrix<T>> gens{gp.isom};
  for (auto const &g : gp.extra_isom) {
    gens.push_back(g);
  }
  return gens;
}

/*
  The seeds: one glued pair at least in every genus of pairs (L, G), L in
  the genus of A, G acting on the factors as given.
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
  size_t n_gen = factors[0].mus.size();
  for (auto const &f : factors) {
    if (f.mus.size() != n_gen) {
      std::cerr << "ISOM_FACTORS: the factors have different numbers of "
                   "generators\n";
      throw TerminalException{1};
    }
  }
  // the pieces of each factor, of an admissible determinant
  for (size_t i = 0; i < factors.size(); i++) {
    IsometryFactor<T> &f = factors[i];
    // a single factor against the rest: always a bound
    std::vector<Tint> dets = *FactorDeterminants(factors, {i}, d);
    if (f.orth) {
      std::vector<MyMatrix<T>> gens;
      for (int s : f.signs) {
        gens.push_back(T(s) * IdentityMat<T>(f.rank));
      }
      for (auto const &D : dets) {
        for (auto const &gs : EnumerateEvenGenusSymbols<Tint>(f.rank, D)) {
          f.pieces.push_back({GenusSymbolRepresentative<T, Tint>(gs, os), gens});
        }
      }
    } else {
      std::vector<std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>> kept;
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
  auto split_gens = [](std::vector<MyMatrix<T>> const &gens)
      -> std::vector<MyMatrix<T>> {
    return std::vector<MyMatrix<T>>(gens.begin() + 1, gens.end());
  };
  // the glued pairs so far, glue_order the product of the glue orders
  std::vector<GluedPair<T>> current;
  for (auto const &[A1, gens1] : factors[0].pieces) {
    current.push_back(
        {A1, gens1[0], DeterminantMat(A1), T(1), split_gens(gens1)});
  }
  std::vector<size_t> idx{0};
  for (size_t j = 1; j < factors.size(); j++) {
    idx.push_back(j);
    bool last = (j + 1 == factors.size());
    // the determinants of M_j, when bounded
    std::optional<std::vector<Tint>> dets;
    if (!last) {
      dets = FactorDeterminants(factors, idx, d);
    }
    std::optional<T> target;
    if (last) {
      target = det;
    }
    // the glue between M_{j-1} and P_j
    std::vector<size_t> idx_prev(idx.begin(), idx.end() - 1);
    GlueBound<Tint> gb = GlueBounds<T, Tint>(factors, idx_prev, {j});
    Tint const &S_glue = gb.order;
    std::set<std::vector<MyMatrix<T>>> known;
    std::vector<GluedPair<T>> next;
    for (auto const &cur : current) {
      for (auto const &[A2, gens2] : factors[j].pieces) {
        T det2 = DeterminantMat(A2);
        if (last && !IsAdmissibleGlue(cur.det, det2, det, S_glue)) {
          continue;
        }
        if (!last && dets) {
          // some admissible determinant of M_j must be reachable
          bool reachable = false;
          for (auto const &D : *dets) {
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
        for (auto const &gl : GlueGroupPairs<T, Tint, Tgroup>(
                 cur.gram, GluedGens(cur), A2, gens2, target,
                 std::optional<Tint>(gb.exponent), bounds, os)) {
          if (!last && dets && !IsAdmissibleDeterminant(*dets, gl.det)) {
            continue;
          }
          if (last && GenusSymbolOfGram<T, Tint>(gl.gram) != S) {
            continue;
          }
          std::vector<MyMatrix<T>> ListP;
          for (auto const &g : GluedGens(gl)) {
            ListP.push_back(g.transpose());
          }
          if (known.insert(PairForms(gl.gram, ListP)).second) {
            GluedPair<T> gp = gl;
            gp.glue_order = cur.glue_order * gl.glue_order;
            next.push_back(gp);
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
  in {0, 1}^n, with v A0 = 0 mod 2, v P = v mod 2 for every P of ListP
  (the generators of the group) and v A0 v^T = 4 mod 8 (w.w odd), one
  overlattice A0 + w per orbit with the orbit size.
 */
template <typename T, typename Tint>
std::vector<OddOverlattice<T>>
OddOverlatticeOrbits(MyMatrix<T> const &A0, std::vector<MyMatrix<T>> const &ListP,
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
    bool fixed = true;
    for (auto const &P : ListP) {
      fixed = fixed && (act(v, P) == v);
    }
    if (!fixed) {
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
