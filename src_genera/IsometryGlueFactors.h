// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_ISOMETRYGLUEFACTORS_H_
#define SRC_GENERA_ISOMETRYGLUEFACTORS_H_

// clang-format off
#include "IsometryGlue.h"
#include "IsotypicPieces.h"
#include "HermitianLattice.h"
#include "CyclotomicArithmetic.h"
#include "MAT_MatrixModPrime.h"
#include "GenusSymbol.h"
#include <algorithm>
#include <fstream>
#include <functional>
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

// Phi_d with coefficients in T (CyclotomicPolynomial of
// CyclotomicArithmetic.h).
template <typename T> std::vector<T> CyclotomicPolynomialCoeffs(int d) {
  std::vector<T> res;
  for (int64_t c : CyclotomicPolynomial(d)) {
    res.push_back(UniversalScalarConversion<T, int64_t>(c));
  }
  return res;
}

/*
  The possible determinants s e (s | S, e | d) of a sigma-stable primitive
  sublattice X of L, det L = d, S a multiple of the order of the glue (e.g.
  Res^min(rank X, rank Y) as above).
 */
template <typename Tint>
std::vector<Tint> SublatticeDeterminants(Tint const &d, Tint const &S) {
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
    e.g. the hermitian genera of GENUS_HermitianPieces, with g_j = g^e_j;
  - irrep: the pieces rho^c, rho an absolutely irreducible rational
    representation of G given by a G-lattice M (IsotypicPieces.h), made
    here;
  - herm: the hermitian genera over Z[zeta_d] of rank b (HermitianPieces of
    HermitianLattice.h), made here, with g_j = zeta^e_j: G acting through a
    linear character of values in Q(zeta_d) and its Galois conjugates.
  For an abelian G the pieces of L are its isotypic parts, on which G acts
  through a character (character = true): the factors must have distinct
  characters. For an irrep factor G does not act through a character, and
  the glue is bounded by the resultants of the generators alone.
 */
template <typename T> struct IsometryFactor {
  std::string desc;
  bool orth;
  bool character = true;
  std::optional<GroupLattice<T>> irrep;
  int multiplicity = 0;
  int64_t herm_d = 0;
  std::vector<int> signs;
  std::vector<int64_t> powers;
  int rank;
  // the Gram matrix and the actions of the generators (on the columns)
  std::vector<std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>> pieces;
  std::vector<std::vector<T>> mus;
};


// The factor "orth s_1 ... s_k r": g_j = s_j Id on the even genera of rank r.
template <typename T>
IsometryFactor<T> MakeOrthFactor(std::vector<int> const &signs, int rank) {
  IsometryFactor<T> f;
  f.desc = "orth";
  f.orth = true;
  for (int s : signs) {
    if (s != 1 && s != -1) {
      std::cerr << "ISOM_FACTORS: a sign must be 1 or -1\n";
      throw TerminalException{1};
    }
    f.signs.push_back(s);
    f.mus.push_back({T(-s), T(1)});
    f.desc += " " + std::to_string(s);
  }
  if (rank <= 0) {
    std::cerr << "ISOM_FACTORS: the rank of an orth factor must be positive\n";
    throw TerminalException{1};
  }
  f.rank = rank;
  f.desc += " " + std::to_string(rank);
  return f;
}

// The factor "file path e_1 ... e_k": g_j = g^e_j on the pieces of the file.
template <typename T>
IsometryFactor<T> MakeFileFactor(std::string const &path,
                                 std::vector<int64_t> const &powers) {
  IsometryFactor<T> f;
  f.desc = "file " + path;
  f.orth = false;
  for (int64_t e : powers) {
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
      gens.push_back(MatrixPower(lst[i + 1], e));
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

/*
  The factor "irrep path c": the file holds the Gram matrix of a G-lattice
  M in rho, then the actions of the n_gen generators on it (on the
  columns); the pieces are made in FactorSeeds.
 */
template <typename T>
IsometryFactor<T> MakeIrrepFactor(std::string const &path, int c, int n_gen) {
  std::vector<MyMatrix<T>> lst = ReadListMatrixFile<T>(path);
  if (static_cast<int>(lst.size()) != n_gen + 1 || c <= 0) {
    std::cerr << "ISOM_FACTORS: " << path << " must hold a Gram matrix and "
              << n_gen << " generators, and c must be positive\n";
    throw TerminalException{1};
  }
  IsometryFactor<T> f;
  f.desc = "irrep " + path + " " + std::to_string(c);
  f.orth = false;
  f.character = false;
  f.multiplicity = c;
  GroupLattice<T> M{lst[0], {}};
  for (int j = 0; j < n_gen; j++) {
    CheckIsometryInput(lst[0], lst[j + 1]);
    M.gens.push_back(lst[j + 1]);
    f.mus.push_back(MinimalPolynomialCoeffs(lst[j + 1]));
  }
  f.irrep = M;
  f.rank = M.gram.rows() * c;
  return f;
}

/*
  The factor "herm d b e_1 ... e_k": the hermitian genera of rank b over
  Z[zeta_d], g_j = zeta_d^e_j, made in FactorSeeds; mu_j = Phi_{d / gcd(d,
  e_j)}.
 */
template <typename T>
IsometryFactor<T> MakeHermFactor(int64_t d, int b,
                                 std::vector<int64_t> const &powers) {
  if (d < 3 || b <= 0) {
    std::cerr << "ISOM_FACTORS: herm needs d >= 3 and b > 0\n";
    throw TerminalException{1};
  }
  IsometryFactor<T> f;
  f.desc = "herm " + std::to_string(d) + " " + std::to_string(b);
  f.orth = false;
  f.herm_d = d;
  f.multiplicity = b;
  for (int64_t e : powers) {
    if (e < 0) {
      std::cerr << "ISOM_FACTORS: the powers must be non negative\n";
      throw TerminalException{1};
    }
    f.powers.push_back(e);
    f.mus.push_back(CyclotomicPolynomialCoeffs<T>(
        static_cast<int>(d / std::gcd(d, e))));
    f.desc += " " + std::to_string(e);
  }
  f.rank = b * static_cast<int>(EulerPhi(d));
  return f;
}

/*
  A factor from one line of a spec file, for a group of n_gen generators:
      orth s_1 ... s_k r
      file <path> e_1 ... e_k   (e_1 = 1 by default when k = 1)
      irrep <path> c
      herm d b e_1 ... e_k
 */
template <typename T>
IsometryFactor<T> ParseIsometryFactor(std::string const &line, int n_gen) {
  std::istringstream ls(line);
  std::string kind;
  ls >> kind;
  if (kind == "orth") {
    std::vector<int> signs(n_gen);
    int rank = 0;
    for (auto &s : signs) {
      if (!(ls >> s)) {
        std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
        throw TerminalException{1};
      }
    }
    if (!(ls >> rank)) {
      std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
      throw TerminalException{1};
    }
    return MakeOrthFactor<T>(signs, rank);
  }
  if (kind == "herm") {
    int64_t d = 0;
    int b = 0;
    std::vector<int64_t> powers(n_gen, 1);
    if (!(ls >> d >> b)) {
      std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
      throw TerminalException{1};
    }
    for (int j = 0; j < n_gen; j++) {
      if (!(ls >> powers[j])) {
        std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
        throw TerminalException{1};
      }
    }
    return MakeHermFactor<T>(d, b, powers);
  }
  if (kind == "irrep") {
    std::string path;
    int c = 0;
    if (!(ls >> path >> c)) {
      std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
      throw TerminalException{1};
    }
    return MakeIrrepFactor<T>(path, c, n_gen);
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
  std::vector<int64_t> powers(n_gen, 1);
  for (int j = 0; j < n_gen; j++) {
    if (!(ls >> powers[j]) && n_gen > 1) {
      std::cerr << "ISOM_FACTORS: bad line: " << line << "\n";
      throw TerminalException{1};
    }
  }
  return MakeFileFactor<T>(path, powers);
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
  The cyclotomic factors of the minimal polynomial mu of an element of
  finite order (a product of distinct cyclotomic polynomials): the n with
  Phi_n | mu, by trial division.
 */
template <typename T> std::set<int> CyclotomicFactorsOf(std::vector<T> mu) {
  std::set<int> res;
  int deg = static_cast<int>(mu.size()) - 1;
  for (int n = 1; deg > 0 && n <= 4 * deg * deg + 6; n++) {
    std::vector<T> phi = CyclotomicPolynomialCoeffs<T>(n);
    int dphi = static_cast<int>(phi.size()) - 1;
    if (dphi > deg) {
      continue;
    }
    // the division of mu by the monic phi
    std::vector<T> rem = mu;
    std::vector<T> quo(deg - dphi + 1, T(0));
    for (int k = deg - dphi; k >= 0; k--) {
      T q = rem[k + dphi];
      quo[k] = q;
      for (int i = 0; i <= dphi; i++) {
        rem[k + i] -= q * phi[i];
      }
    }
    bool zero = true;
    for (int i = 0; i < dphi; i++) {
      zero = zero && rem[i] == 0;
    }
    if (zero) {
      res.insert(n);
      mu = quo;
      deg -= dphi;
    }
  }
  if (deg != 0) {
    std::cerr << "ISOM_FACTORS: a minimal polynomial is not a product of "
                 "cyclotomic polynomials\n";
    throw TerminalException{1};
  }
  return res;
}

/*
  The minimal polynomial of g_j on the sum of the factors in idx: the lcm
  of theirs, the product of the distinct cyclotomic polynomials dividing
  one of them.
 */
template <typename T>
std::vector<T> FactorsProduct(std::vector<IsometryFactor<T>> const &factors,
                              std::vector<size_t> const &idx, size_t j) {
  std::set<int> cyc;
  for (size_t i : idx) {
    for (int n : CyclotomicFactorsOf(factors[i].mus[j])) {
      cyc.insert(n);
    }
  }
  std::vector<T> p{T(1)};
  for (int n : cyc) {
    p = PolyProductCoeffs(p, CyclotomicPolynomialCoeffs<T>(n));
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
  return EulerPhi(ord);
}

/*
  The order of the image of G acting on the sum of the factors, from a
  G-lattice of each (the signs of an orth factor, M for an irrep factor,
  zeta_d for a herm factor, the first piece otherwise), by closure of the
  block diagonal generators.
 */
template <typename T, typename Tint>
Tint FactorsGroupOrder(std::vector<IsometryFactor<T>> const &factors) {
  size_t n_gen = factors[0].mus.size();
  std::vector<std::vector<MyMatrix<T>>> blocks(n_gen);
  for (auto const &f : factors) {
    for (size_t j = 0; j < n_gen; j++) {
      if (f.orth) {
        blocks[j].push_back(T(f.signs[j]) * IdentityMat<T>(1));
      } else if (f.irrep) {
        blocks[j].push_back(f.irrep->gens[j]);
      } else if (f.herm_d != 0) {
        // zeta_d acting on Z[zeta_d] by the companion matrix of Phi_d
        std::vector<T> phi = CyclotomicPolynomialCoeffs<T>(
            static_cast<int>(f.herm_d));
        int k = static_cast<int>(phi.size()) - 1;
        MyMatrix<T> C = ZeroMatrix<T>(k, k);
        for (int i = 0; i + 1 < k; i++) {
          C(i + 1, i) = T(1);
        }
        for (int i = 0; i < k; i++) {
          C(i, k - 1) = -phi[i];
        }
        blocks[j].push_back(MatrixPower(C, f.powers[j]));
      } else {
        blocks[j].push_back(f.pieces[0].second[j]);
      }
    }
  }
  std::vector<MyMatrix<T>> ListP;
  for (auto const &bl : blocks) {
    MyMatrix<T> P = bl[0];
    for (size_t k = 1; k < bl.size(); k++) {
      P = BlockDiagonal(P, bl[k]);
    }
    ListP.push_back(P);
  }
  size_t order = GroupElementsAndInverses(ListP, 1000000).first.size();
  return Tint(static_cast<long>(order));
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
  bool all_character = true;
  for (auto const *idx : {&idx_in, &idx_out}) {
    for (size_t i : *idx) {
      all_character = all_character && factors[i].character;
    }
  }
  if (!all_character) {
    // the glue is killed by the resultant of any generator separating the
    // two sides, and is generated by min(rank X, rank Y) elements
    Tint expo(0), order(0);
    for (size_t j = 0; j < n_gen; j++) {
      T res = T_abs(PolyResultant(FactorsProduct(factors, idx_in, j),
                                  FactorsProduct(factors, idx_out, j)));
      if (res != 0) {
        Tint res_i = UniversalScalarConversion<Tint, T>(res);
        Tint ord_j = MyPow<Tint>(res_i, std::min(rank_in, rank_out));
        expo = (expo == 0) ? res_i : GcdPair(expo, res_i);
        order = (order == 0) ? ord_j : GcdPair(order, ord_j);
      }
    }
    if (expo == 0) {
      // no generator separates them: the glue is killed by the order N of
      // the image H of G, as N e_X is in Z[H] for the central idempotent
      // e_X, acting as N on X and 0 on Y
      Tint N = FactorsGroupOrder<T, Tint>(factors);
      return {N, MyPow<Tint>(N, std::min(rank_in, rank_out))};
    }
    return {expo, order};
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
  return SublatticeDeterminants(
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
    std::vector<size_t> rest;
    for (size_t k = 0; k < factors.size(); k++) {
      if (k != i) {
        rest.push_back(k);
      }
    }
    Tint S_i = GlueBounds<T, Tint>(factors, {i}, rest).order;
    std::vector<Tint> dets = SublatticeDeterminants(d, S_i);
    // at an odd prime outside the glue the piece is an orthogonal summand
    // of L locally: its odd symbol fits in that of L
    std::set<int64_t> glue_primes;
    if (S_i > 1) {
      for (auto const &[q, mult] : FactorsIntMap(S_i)) {
        glue_primes.insert(UniversalScalarConversion<int64_t, Tint>(q));
      }
    }
    if (f.orth) {
      std::vector<MyMatrix<T>> gens;
      for (int s : f.signs) {
        gens.push_back(T(s) * IdentityMat<T>(f.rank));
      }
      for (auto const &D : dets) {
        for (auto const &gs : EnumerateEvenGenusSymbols<Tint>(f.rank, D)) {
          if (OddSymbolFitsIn<Tint>(gs, S, glue_primes)) {
            f.pieces.push_back(
                {GenusSymbolRepresentative<T, Tint>(gs, os), gens});
          }
        }
      }
    } else if (f.irrep) {
      // the complete searches: the twists when End_G(M) is a field with a
      // self-adjoint part larger than Q, End(M) for Schur index 1, the
      // maximal orders otherwise; none applies (a maximal order beyond the
      // bounds, a quaternion algebra over K^+ != Q) is an error, as a
      // search that may miss pieces would miss genera without notice
      std::vector<Tint> gp = GroupLatticePrimes<T, Tint>(*f.irrep);
      std::optional<std::vector<GroupLattice<T>>> cert =
          IsotypicPiecesTwisted<T, Tint, Tgroup>(*f.irrep, f.multiplicity,
                                                 dets, os);
      if (!cert) {
        cert = IsotypicPiecesCertified<T, Tint, Tgroup>(
            *f.irrep, f.multiplicity, dets, gp, os);
      }
      if (!cert) {
        cert = IsotypicPiecesMaximalOrder<T, Tint, Tgroup>(
            *f.irrep, f.multiplicity, dets, gp, os);
      }
      if (!cert) {
        std::cerr << "ISOM_FACTORS: " << f.desc
                  << ": no complete search of the pieces for this "
                     "representation\n";
        throw TerminalException{1};
      }
      std::vector<GroupLattice<T>> const &list_P = *cert;
      for (auto const &P : list_P) {
        if (OddSymbolFitsIn<Tint>(GenusSymbolOfGram<T, Tint>(P.gram), S,
                                  glue_primes)) {
          f.pieces.push_back({P.gram, P.gens});
        }
      }
    } else {
      if (f.herm_d != 0) {
        Tint N(1);
        for (auto const &D : dets) {
          N = LCMpair(N, D);
        }
        for (auto const &[G, P] :
             HermitianPieces<T, Tint>(f.herm_d, f.multiplicity, N, os)) {
          MyMatrix<T> g = P.transpose();
          std::vector<MyMatrix<T>> gens;
          for (int64_t e : f.powers) {
            gens.push_back(MatrixPower(g, e));
          }
          f.pieces.push_back({G, gens});
        }
      }
      std::vector<std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>> kept;
      for (auto const &pc : f.pieces) {
        if (IsAdmissibleDeterminant(dets, DeterminantMat(pc.first)) &&
            OddSymbolFitsIn<Tint>(GenusSymbolOfGram<T, Tint>(pc.first), S,
                                  glue_primes)) {
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
  /*
    The gluing, as a tree: the sum of the factors of a set is glued from
    the sums of two parts.
    - A group: the parts are separated by a generator g_j (the factors of
      one minimal polynomial of g_j against the others), so that the glue
      between them is bounded by the resultant for g_j; the most balanced
      split by rank is taken. The large factors then do not meet every
      intermediate pair (a case of Z/2 x Z/2 on A2^4 was beyond 10 minutes
      as a chain).
    - A single isometry: a chain by increasing number of pieces, except
      that the whole lattice glues last the factor with the fewest pieces:
      the last step, with the target determinant and genus, meets every
      intermediate pair (an I9 case: 638 seconds to 2).
    Every intermediate sum is primitive in L: its determinant is
    admissible and its odd symbol fits in that of L outside the primes of
    its glue with the complement.
   */
  IsometryGlueBounds bounds;
  auto split_gens = [](std::vector<MyMatrix<T>> const &gens)
      -> std::vector<MyMatrix<T>> {
    return std::vector<MyMatrix<T>>(gens.begin() + 1, gens.end());
  };
  auto primes_of = [](Tint const &x) -> std::set<int64_t> {
    std::set<int64_t> res;
    if (x > 1) {
      for (auto const &[q, mult] : FactorsIntMap(x)) {
        res.insert(UniversalScalarConversion<int64_t, Tint>(q));
      }
    }
    return res;
  };
  auto rank_of = [&](std::vector<size_t> const &idx) -> int {
    int r = 0;
    for (size_t i : idx) {
      r += factors[i].rank;
    }
    return r;
  };
  std::vector<size_t> all(factors.size());
  for (size_t i = 0; i < factors.size(); i++) {
    all[i] = i;
  }
  // the split of idx into two parts
  auto split = [&](std::vector<size_t> const &idx, bool top)
      -> std::pair<std::vector<size_t>, std::vector<size_t>> {
    if (n_gen > 1) {
      std::optional<std::pair<std::vector<size_t>, std::vector<size_t>>> best;
      int best_gap = 0;
      int rank_all = rank_of(idx);
      for (size_t j = 0; j < n_gen; j++) {
        std::map<std::vector<T>, std::vector<size_t>> classes;
        for (size_t i : idx) {
          classes[factors[i].mus[j]].push_back(i);
        }
        if (classes.size() < 2) {
          continue;
        }
        for (auto const &[mu, part] : classes) {
          std::vector<size_t> other;
          for (size_t i : idx) {
            if (std::find(part.begin(), part.end(), i) == part.end()) {
              other.push_back(i);
            }
          }
          int gap = std::abs(rank_all - 2 * rank_of(part));
          if (!best || gap < best_gap) {
            best = std::make_pair(part, other);
            best_gap = gap;
          }
        }
      }
      if (best) {
        return *best;
      }
    }
    // a chain: the factor of the fewest pieces last at the top, the one of
    // the most pieces last below
    auto cmp = [&](size_t x, size_t y) -> bool {
      return factors[x].pieces.size() < factors[y].pieces.size();
    };
    size_t last = top ? *std::min_element(idx.begin(), idx.end(), cmp)
                      : *std::max_element(idx.begin(), idx.end(), cmp);
    std::vector<size_t> first;
    for (size_t i : idx) {
      if (i != last) {
        first.push_back(i);
      }
    }
    return {first, {last}};
  };
  std::function<std::vector<GluedPair<T>>(std::vector<size_t> const &, bool)>
      seeds = [&](std::vector<size_t> const &idx,
                  bool top) -> std::vector<GluedPair<T>> {
    if (idx.size() == 1) {
      std::vector<GluedPair<T>> res;
      for (auto const &[A1, gens1] : factors[idx[0]].pieces) {
        res.push_back(
            {A1, gens1[0], DeterminantMat(A1), T(1), split_gens(gens1)});
      }
      return res;
    }
    auto [idx_a, idx_b] = split(idx, top);
    std::vector<GluedPair<T>> seeds_a = seeds(idx_a, false);
    std::vector<GluedPair<T>> seeds_b = seeds(idx_b, false);
    GlueBound<Tint> gb = GlueBounds<T, Tint>(factors, idx_a, idx_b);
    std::set<int64_t> glue_primes = primes_of(gb.order);
    // the conditions on the sum of idx, primitive in L
    std::optional<std::vector<Tint>> dets;
    std::set<int64_t> sum_primes;
    if (!top) {
      dets = FactorDeterminants(factors, idx, d);
      std::vector<size_t> rest;
      for (size_t i : all) {
        if (std::find(idx.begin(), idx.end(), i) == idx.end()) {
          rest.push_back(i);
        }
      }
      sum_primes = primes_of(GlueBounds<T, Tint>(factors, idx, rest).order);
    }
    std::optional<T> target;
    if (top) {
      target = det;
    }
    std::vector<GenusSymbol> sym_a, sym_b;
    if (top) {
      for (auto const &x : seeds_a) {
        sym_a.push_back(GenusSymbolOfGram<T, Tint>(x.gram));
      }
      for (auto const &x : seeds_b) {
        sym_b.push_back(GenusSymbolOfGram<T, Tint>(x.gram));
      }
    }
    std::set<std::vector<MyMatrix<T>>> known;
    std::vector<GluedPair<T>> res;
    for (size_t ia = 0; ia < seeds_a.size(); ia++) {
      GluedPair<T> const &a = seeds_a[ia];
      for (size_t ib = 0; ib < seeds_b.size(); ib++) {
        GluedPair<T> const &b = seeds_b[ib];
        if (top) {
          if (!IsAdmissibleGlue(a.det, b.det, det, gb.order) ||
              !OddSymbolsAddUp<Tint>(sym_a[ia], sym_b[ib], S, glue_primes)) {
            continue;
          }
        } else if (dets) {
          // some admissible determinant of the sum must be reachable
          bool reachable = false;
          for (auto const &D : *dets) {
            if (IsAdmissibleGlue(a.det, b.det,
                                 UniversalScalarConversion<T, Tint>(D),
                                 gb.order)) {
              reachable = true;
              break;
            }
          }
          if (!reachable) {
            continue;
          }
        }
        for (auto const &gl : GlueGroupPairs<T, Tint, Tgroup>(
                 a.gram, GluedGens(a), b.gram, GluedGens(b), target,
                 std::optional<Tint>(gb.exponent), bounds, os)) {
          if (top) {
            if (GenusSymbolOfGram<T, Tint>(gl.gram) != S) {
              continue;
            }
          } else {
            if (dets && !IsAdmissibleDeterminant(*dets, gl.det)) {
              continue;
            }
            if (!OddSymbolFitsIn<Tint>(GenusSymbolOfGram<T, Tint>(gl.gram), S,
                                       sum_primes)) {
              continue;
            }
          }
          std::vector<MyMatrix<T>> ListP;
          for (auto const &g : GluedGens(gl)) {
            ListP.push_back(g.transpose());
          }
          if (known.insert(PairForms(gl.gram, ListP)).second) {
            GluedPair<T> gp = gl;
            gp.glue_order = a.glue_order * b.glue_order * gl.glue_order;
            res.push_back(gp);
          }
        }
      }
    }
    std::cerr << "ISOM_FACTORS: " << idx.size() << " factors (rank "
              << rank_of(idx) << "), " << res.size() << " glued pairs\n";
    return res;
  };
  std::vector<GluedPair<T>> res;
  for (auto const &gp : seeds(all, true)) {
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
  // the kernel of A0 mod 2 (A0 symmetric: x A0 = 0 iff A0 x = 0)
  ModMatrix K = LeftNullspaceModP(ModMatrixOfMatrix(A0, 2), 2);
  std::vector<std::vector<int>> kernel;
  for (int b = 0; b < K.rows(); b++) {
    std::vector<int> v(n);
    for (int i = 0; i < n; i++) {
      v[i] = static_cast<int>(K(b, i));
    }
    kernel.push_back(v);
  }
  // the vectors are keyed by 64-bit masks
  if (n > 64) {
    std::cerr << "ISOM_FACTORS: the odd overlattices are for rank <= 64\n";
    throw TerminalException{1};
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
  // v -> v h mod 2, for the matrices reduced mod 2
  auto act = [&](std::vector<int> const &v,
                 ModMatrix const &h) -> std::vector<int> {
    std::vector<int> w(n, 0);
    for (int j = 0; j < n; j++) {
      T_residue s = 0;
      for (int i = 0; i < n; i++) {
        if (v[i]) {
          s ^= h(i, j);
        }
      }
      w[j] = static_cast<int>(s);
    }
    return w;
  };
  std::vector<ModMatrix> ListP2;
  for (auto const &P : ListP) {
    ListP2.push_back(ModMatrixOfMatrix(P, 2));
  }
  std::vector<ModMatrix> gens;
  for (auto const &h : ListGen) {
    gens.push_back(ModMatrixOfMatrix(h, 2));
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
    for (auto const &P2 : ListP2) {
      fixed = fixed && (act(v, P2) == v);
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
