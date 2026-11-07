// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSENUMERATIONGROUP_H_
#define SRC_GENERA_GENUSENUMERATIONGROUP_H_

// clang-format off
#include "GenusEnumeration.h"
#include <algorithm>
#include <map>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_GENUS_ENUMERATION_GROUP
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_GENUS_ENUMERATION_GROUP
#endif

#ifdef TIMINGS
#define TIMINGS_GENUS_ENUMERATION_GROUP
#endif

/*
  Enumeration of the classes of pairs (L, sigma) in the genus of a pair:
  L even positive definite, sigma an isometry of L of finite order m, two
  pairs isometric when an isometry of the lattices conjugates the
  isometries. The termination is by the mass

      mass = sum over the classes of 1 / |C_{O(L)}(sigma)| ,

  which GENUS_IsometryMass computes. This is the cyclic case of the
  enumeration of the lattices of a genus with a prescribed symmetry
  (Section 6 of the paper on the genus enumeration).

  Conventions. Internally sigma acts on the row vectors by x -> x P, with P
  A P^T = A; on input and output it is the matrix g = P^T acting on the
  columns, g^T A g = A.

  The neighbours. At a good prime p (p not dividing m det), L (x) Z_p is
  the direct sum of the eigenlattices of sigma for the eigenvalues in Z_p
  and of the other primary components. The sigma-stable p-neighbours are
  the neighbours N = L_v + Z v/p of the eigenvectors v:
  - for lambda != +-1 every eigenvector is isotropic;
  - for lambda = +-1 the isotropic ones.
  The vector has to be an eigenvector modulo p^2, not only modulo p. It is
  obtained by the projector E_lambda = R(P) / R(lambda'), lambda' the
  Teichmueller lift of lambda to Z/p^2 and R = (x^m - 1) / (x - lambda'):
  x^m - 1 kills P, and R(lambda') = m lambda'^{m-1} is a unit since p does
  not divide m. For lambda = +-1 the correction to an isotropic vector
  modulo p^2 is done inside the eigenlattice, which keeps it an
  eigenvector. The centralizer C_{O(L)}(sigma) acts on each eigenspace
  modulo p, and the neighbours are taken over its orbits on the projective
  space of the eigenspace, by ProjectiveLineOrbitRepresentatives.

  The classes. A pair is told apart by the canonical form of the pair of
  forms (A, A P^T) (the second one not symmetric), and C_{O(L)}(sigma) is
  the automorphism group of that pair of forms.
 */

template <typename T, typename Tgroup> struct GenusEnumerationGroupResult {
  using TintGroup = typename Tgroup::Tint;
  using Tmass = typename overlying_field<TintGroup>::field_type;
  std::vector<MyMatrix<T>> ListGram;
  // the isometries, acting on the columns: g^T A g = A
  std::vector<MyMatrix<T>> ListIsom;
  std::vector<TintGroup> ListCentralizerOrder;
  Tmass accumulated_mass;
  Tmass target_mass;
  bool complete;
  std::vector<int> ListPrimeUsed;
  std::vector<size_t> ListNewClassPerPrime;
};

inline int64_t ModPosInt64(int64_t a, int64_t m) {
  a %= m;
  return (a < 0) ? a + m : a;
}

inline int64_t PowModInt64(int64_t a, int64_t e, int64_t m) {
  int64_t res = 1 % m;
  a = ModPosInt64(a, m);
  while (e > 0) {
    if (e & 1) {
      res = (res * a) % m;
    }
    a = (a * a) % m;
    e >>= 1;
  }
  return res;
}

inline int64_t InverseModInt64(int64_t u, int64_t m) {
  int64_t a = ModPosInt64(u, m), b = m, x0 = 1, x1 = 0;
  while (b != 0) {
    int64_t q = a / b;
    int64_t t = a - q * b;
    a = b;
    b = t;
    t = x0 - q * x1;
    x0 = x1;
    x1 = t;
  }
  if (a != 1) {
    std::cerr << "GENUS_GROUP: " << u << " is not invertible modulo " << m
              << "\n";
    throw TerminalException{1};
  }
  return ModPosInt64(x0, m);
}

using ModMatrix = std::vector<std::vector<int64_t>>;

inline ModMatrix ModMatrixProduct(ModMatrix const &X, ModMatrix const &Y,
                                  int64_t mod) {
  size_t n = X.size();
  size_t k = Y.size();
  size_t l = Y.empty() ? 0 : Y[0].size();
  ModMatrix Z(n, std::vector<int64_t>(l, 0));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < k; j++) {
      int64_t x = X[i][j];
      if (x == 0) {
        continue;
      }
      for (size_t u = 0; u < l; u++) {
        Z[i][u] = (Z[i][u] + x * Y[j][u]) % mod;
      }
    }
  }
  return Z;
}

/*
  The reduced row echelon form of the row space modulo the prime p: the
  basis rows and their pivot columns. The coordinates of a vector of the
  row space in that basis are its entries at the pivots.
 */
struct RowSpaceModP {
  ModMatrix basis;
  std::vector<int> pivots;
};

inline RowSpaceModP ComputeRowSpaceModP(ModMatrix rows, int64_t p) {
  RowSpaceModP res;
  size_t n_col = rows.empty() ? 0 : rows[0].size();
  size_t r = 0;
  for (auto &row : rows) {
    for (auto &x : row) {
      x = ModPosInt64(x, p);
    }
  }
  for (size_t col = 0; col < n_col && r < rows.size(); col++) {
    size_t piv = r;
    while (piv < rows.size() && rows[piv][col] == 0) {
      piv++;
    }
    if (piv == rows.size()) {
      continue;
    }
    std::swap(rows[r], rows[piv]);
    int64_t inv = InverseModInt64(rows[r][col], p);
    for (auto &x : rows[r]) {
      x = (x * inv) % p;
    }
    for (size_t i = 0; i < rows.size(); i++) {
      if (i != r && rows[i][col] != 0) {
        int64_t c = rows[i][col];
        for (size_t j = 0; j < n_col; j++) {
          rows[i][j] = ModPosInt64(rows[i][j] - c * rows[r][j], p);
        }
      }
    }
    res.pivots.push_back(static_cast<int>(col));
    r++;
  }
  rows.resize(r);
  res.basis = std::move(rows);
  return res;
}

/*
  The data of one eigenvalue at p: lambda modulo p, its lift modulo p^2,
  the projector E_lambda modulo p^2 (on the row vectors) and the eigenspace
  modulo p.
 */
struct EigenspaceData {
  int64_t lambda;
  int64_t lambda2;
  ModMatrix E;
  RowSpaceModP space;
};

template <typename Tint>
ModMatrix ReduceMatrixMod(MyMatrix<Tint> const &M, int64_t mod) {
  int n = M.rows();
  int l = M.cols();
  ModMatrix R(n, std::vector<int64_t>(l));
  Tint mod_T = UniversalScalarConversion<Tint, int64_t>(mod);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < l; j++) {
      Tint r = M(i, j) % mod_T;
      R[i][j] = ModPosInt64(UniversalScalarConversion<int64_t, Tint>(r), mod);
    }
  }
  return R;
}

template <typename Tint>
std::vector<EigenspaceData> ComputeEigenspaces(MyMatrix<Tint> const &P, int m,
                                               int p) {
  std::vector<EigenspaceData> res;
  int n = P.rows();
  int64_t p2 = static_cast<int64_t>(p) * p;
  ModMatrix P2 = ReduceMatrixMod(P, p2);
  for (int64_t lambda = 1; lambda < p; lambda++) {
    if (PowModInt64(lambda, m, p) != 1) {
      continue;
    }
    int64_t lambda2 = PowModInt64(lambda, p, p2);
    // R = (x^m - 1) / (x - lambda2) = sum_k r_k x^k, r_{m-1} = 1,
    // r_{k-1} = lambda2 r_k
    std::vector<int64_t> r(m, 0);
    r[m - 1] = 1;
    for (int k = m - 1; k >= 1; k--) {
      r[k - 1] = (lambda2 * r[k]) % p2;
    }
    // R(P) by Horner
    ModMatrix RP(n, std::vector<int64_t>(n, 0));
    for (int k = m - 1; k >= 0; k--) {
      RP = ModMatrixProduct(RP, P2, p2);
      for (int i = 0; i < n; i++) {
        RP[i][i] = (RP[i][i] + r[k]) % p2;
      }
    }
    int64_t Rlambda = (static_cast<int64_t>(m) % p2) *
                      PowModInt64(lambda2, m - 1, p2) % p2;
    int64_t inv = InverseModInt64(Rlambda, p2);
    for (auto &row : RP) {
      for (auto &x : row) {
        x = (x * inv) % p2;
      }
    }
    RowSpaceModP space = ComputeRowSpaceModP(RP, p);
    if (space.basis.empty()) {
      continue;
    }
    res.push_back({lambda, lambda2, std::move(RP), std::move(space)});
  }
  return res;
}

// The action of h (on the row vectors) on the eigenspace, in its basis.
template <typename Tint>
MyMatrix<Tint> EigenspaceAction(EigenspaceData const &ed, MyMatrix<Tint> const &h,
                                int p) {
  size_t k = ed.space.basis.size();
  int n = h.rows();
  ModMatrix hp = ReduceMatrixMod(h, p);
  MyMatrix<Tint> M(k, k);
  for (size_t i = 0; i < k; i++) {
    std::vector<int64_t> img(n, 0);
    for (int a = 0; a < n; a++) {
      int64_t x = ed.space.basis[i][a];
      if (x == 0) {
        continue;
      }
      for (int b = 0; b < n; b++) {
        img[b] = (img[b] + x * hp[a][b]) % p;
      }
    }
    for (size_t j = 0; j < k; j++) {
      M(i, j) = UniversalScalarConversion<Tint, int64_t>(img[ed.space.pivots[j]]);
    }
  }
  return M;
}

/*
  The eigenvector modulo p^2 for the coefficient vector c of the
  eigenspace, made isotropic modulo p^2. Nothing for an anisotropic line of
  the eigenvalues +-1.
 */
template <typename Tint>
std::optional<std::vector<int>>
EigenvectorNeighborVector(EigenspaceData const &ed, std::vector<int> const &c,
                          MyMatrix<Tint> const &G, int p) {
  int n = G.rows();
  int64_t p2 = static_cast<int64_t>(p) * p;
  size_t k = ed.space.basis.size();
  std::vector<int64_t> x(n, 0);
  for (size_t i = 0; i < k; i++) {
    if (c[i] == 0) {
      continue;
    }
    for (int a = 0; a < n; a++) {
      x[a] = (x[a] + c[i] * ed.space.basis[i][a]) % p;
    }
  }
  auto apply_E = [&](std::vector<int64_t> const &y) -> std::vector<int64_t> {
    std::vector<int64_t> z(n, 0);
    for (int a = 0; a < n; a++) {
      if (y[a] == 0) {
        continue;
      }
      for (int b = 0; b < n; b++) {
        z[b] = (z[b] + y[a] * ed.E[a][b]) % p2;
      }
    }
    return z;
  };
  auto bilinear = [&](std::vector<int64_t> const &u,
                      std::vector<int64_t> const &v) -> Tint {
    Tint sum(0);
    for (int a = 0; a < n; a++) {
      if (u[a] == 0) {
        continue;
      }
      for (int b = 0; b < n; b++) {
        sum += UniversalScalarConversion<Tint, int64_t>(u[a]) * G(a, b) *
               UniversalScalarConversion<Tint, int64_t>(v[b]);
      }
    }
    return sum;
  };
  std::vector<int64_t> v = apply_E(x);
  Tint q = bilinear(v, v);
  // e = v_p(2): the neighbour needs v_p((v, v)) >= e + 2, the isotropy
  // modulo p is v_p >= e + 1, and the valuation e + 1 is corrected.
  int e = (p == 2) ? 1 : 0;
  bool is_pm1 = (ed.lambda == 1 || ed.lambda == p - 1);
  int val = (q == 0) ? e + 2 : ValuationP(q, p);
  if (val < e + 1) {
    if (!is_pm1) {
      std::cerr << "GENUS_GROUP: an eigenvector of an eigenvalue != +-1 is "
                   "not isotropic\n";
      throw TerminalException{1};
    }
    return {};
  }
  if (val == e + 1) {
    if (!is_pm1) {
      std::cerr << "GENUS_GROUP: an eigenvector of an eigenvalue != +-1 is "
                   "not isotropic modulo p^2\n";
      throw TerminalException{1};
    }
    /*
      v <- v - a p u with u an eigenvector modulo p^2 and (v, u) a unit,
      a = (q / p^{e+1}) / (2 (v, u) / p^e) modulo p: then
      (v - a p u, v - a p u) = q - 2 a p (v, u) + a^2 p^2 (u, u) has
      valuation >= e + 2, the last term because (u, u) is even. The
      correction stays in the eigenlattice, so v remains an eigenvector.
     */
    std::optional<std::vector<int64_t>> opt_u;
    Tint bvu;
    for (size_t i = 0; i < k && !opt_u; i++) {
      std::vector<int64_t> u = apply_E(ed.space.basis[i]);
      Tint b = bilinear(v, u);
      if (ResidueModP(b, p) != 0) {
        opt_u = u;
        bvu = b;
      }
    }
    if (!opt_u) {
      // v is in the radical of the eigenlattice modulo p, impossible for a
      // unimodular L_p
      return {};
    }
    Tint p_T(p);
    Tint q_red = q;
    for (int i = 0; i < e + 1; i++) {
      q_red /= p_T;
    }
    Tint twob_red = Tint(2) * bvu;
    for (int i = 0; i < e; i++) {
      twob_red /= p_T;
    }
    int64_t a = (static_cast<int64_t>(ResidueModP(q_red, p)) *
                 InverseModInt64(ResidueModP(twob_red, p), p)) %
                p;
    for (int i = 0; i < n; i++) {
      v[i] = ModPosInt64(v[i] - a * p * (*opt_u)[i], p2);
    }
#ifdef SANITY_CHECK_GENUS_ENUMERATION_GROUP
    Tint q_new = bilinear(v, v);
    if (q_new != 0 && ValuationP(q_new, p) < e + 2) {
      std::cerr << "GENUS_GROUP: the correction to an isotropic vector "
                   "failed\n";
      throw TerminalException{1};
    }
#endif
  }
  std::vector<int> v_int(n);
  for (int i = 0; i < n; i++) {
    v_int[i] = static_cast<int>(v[i]);
  }
  return v_int;
}

/*
  A neighbour pair: the reduced Gram matrix of N and the isometry carried
  over to its reduced basis (on the row vectors).
 */
template <typename T, typename Tint>
std::optional<std::pair<MyMatrix<T>, MyMatrix<T>>>
GenusNeighborPair(MyMatrix<Tint> const &G, MyMatrix<T> const &P,
                  std::vector<int> const &v, int p, std::ostream &os) {
  std::optional<GenusNeighborBasis<T, Tint>> opt =
      GenusNeighborKernel<T, Tint>(G, v, p, os);
  if (!opt) {
    return {};
  }
  LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(opt->GramN, os);
  // the reduced basis of N in the coordinates of L, times p
  MyMatrix<T> Bred = UniversalMatrixConversion<T, Tint>(rec.Pmat * opt->Bp);
  MyMatrix<T> S = Bred * P * Inverse(Bred);
  if (!IsIntegralMatrix(S)) {
    std::cerr << "GENUS_GROUP: the neighbour is not stable under the "
                 "isometry\n";
    throw TerminalException{1};
  }
  return std::make_pair(rec.GramMatRed, S);
}

template <typename T, typename Tint, typename Tgroup> struct PairAutInfo {
  typename Tgroup::Tint order;
  std::vector<MyMatrix<Tint>> ListGenMat;
};

template <typename T> std::vector<MyMatrix<T>> PairForms(MyMatrix<T> const &A,
                                                         MyMatrix<T> const &P) {
  return {A, A * P.transpose()};
}

template <typename T, typename Tint, typename Tgroup>
PairAutInfo<T, Tint, Tgroup> GetPairAutInfo(MyMatrix<T> const &A,
                                            MyMatrix<T> const &P,
                                            std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat = PairForms(A, P);
  IntAutomorphismGroup<Tint, Tgroup> aut =
      ArithmeticAutomorphismGroupMultiple<T, Tint, Tgroup>(ListMat, os);
  PairAutInfo<T, Tint, Tgroup> info;
  info.order = aut.order;
  auto preserves = [&](MyMatrix<Tint> const &h) -> bool {
    MyMatrix<T> h_T = UniversalMatrixConversion<T, Tint>(h);
    for (auto &M : ListMat) {
      if (h_T * M * h_T.transpose() != M) {
        return false;
      }
    }
    return true;
  };
  for (auto &h : aut.ListGen) {
    if (preserves(h)) {
      info.ListGenMat.push_back(h);
    } else if (preserves(h.transpose())) {
      info.ListGenMat.push_back(h.transpose());
    } else {
      std::cerr << "GENUS_GROUP: an automorphism of the pair is wrong\n";
      throw TerminalException{1};
    }
  }
  return info;
}

// The canonical representative of the pair: the reduced forms and the
// isometry in the canonical basis.
template <typename T, typename Tint, typename Tgroup>
std::pair<MyMatrix<T>, MyMatrix<T>> CanonicalPair(MyMatrix<T> const &A,
                                                  MyMatrix<T> const &P,
                                                  std::ostream &os) {
  std::vector<MyMatrix<T>> ListMat = PairForms(A, P);
  MyMatrix<Tint> B = ComputeCanonicalFormMultiple<T, Tint, Tgroup>(ListMat, os);
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  MyMatrix<T> Acan = B_T * A * B_T.transpose();
  MyMatrix<T> Pcan = B_T * P * Inverse(B_T);
  return {Acan, Pcan};
}

template <typename Tint> int FiniteMatrixOrderInt(MyMatrix<Tint> const &P) {
  int n = P.rows();
  MyMatrix<Tint> Id = IdentityMat<Tint>(n);
  MyMatrix<Tint> pw = P;
  for (int k = 1; k <= 100000; k++) {
    if (pw == Id) {
      return k;
    }
    pw = pw * P;
  }
  std::cerr << "GENUS_GROUP: the isometry is not of finite order <= 100000\n";
  throw TerminalException{1};
}

/*
  A good prime: p does not divide m det. Then L_p is unimodular and sigma is
  semisimple modulo p. The prime 2 is good when m and det are odd; there
  only the eigenvalue 1 is in F_2, but the projective space of the fixed
  space over F_2 is the cheapest of all.
 */
template <typename T> bool IsGoodNeighborPrime(int p, int m, T const &det) {
  if (m % p == 0) {
    return false;
  }
  T p_T(p);
  return ResInt(det, p_T) != T(0);
}

inline int NextPrimeAfter(int q) {
  int cand = q + 1;
  while (true) {
    bool is_prime = cand >= 2;
    for (int d = 2; d * d <= cand; d++) {
      if (cand % d == 0) {
        is_prime = false;
        break;
      }
    }
    if (is_prime) {
      return cand;
    }
    cand++;
  }
}

/*
  The enumeration. ListSeed holds pairs (Gram, isometry on the columns).
  With prime != 0 only that prime is used (any good prime). With prime = 0
  the prime 2 when it is good, then the good primes p = 1 mod m, at which
  every eigenvalue of sigma is in F_p, are brought in in increasing order as in the iterate strategy of
  GenusEnumeration, at most max_prime of them. At a prime p != 1 mod m only
  the eigenvalues in F_p give neighbours, and the projective space of the
  fixed space can be large while the useful lines are elsewhere.
 */
template <typename T, typename Tint, typename Tgroup>
GenusEnumerationGroupResult<T, Tgroup> GenusEnumerationGroup(
    std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> const &ListSeed,
    typename GenusEnumerationGroupResult<T, Tgroup>::Tmass const &TotalMass,
    int prime, size_t max_prime, std::ostream &os) {
  using Tresult = GenusEnumerationGroupResult<T, Tgroup>;
  using Tmass = typename Tresult::Tmass;
  using TintGroup = typename Tgroup::Tint;
#ifdef TIMINGS_GENUS_ENUMERATION_GROUP
  MicrosecondTime time_total;
#endif
  Tresult result;
  result.target_mass = TotalMass;
  result.accumulated_mass = Tmass(0);
  result.complete = false;
  // the work list: the canonical pairs (Gram, isometry on the rows)
  std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> ListPairWork;
  std::vector<PairAutInfo<T, Tint, Tgroup>> ListAutInfo;
  std::map<std::pair<MyMatrix<T>, MyMatrix<T>>, size_t> MapCanonic;
  int m = -1;
  T det(0);
  auto f_insert = [&](MyMatrix<T> const &A, MyMatrix<T> const &P) -> bool {
    std::pair<MyMatrix<T>, MyMatrix<T>> can =
        CanonicalPair<T, Tint, Tgroup>(A, P, os);
    std::pair<MyMatrix<T>, MyMatrix<T>> key{can.first,
                                            can.first * can.second.transpose()};
    if (MapCanonic.count(key) == 1) {
      return false;
    }
    PairAutInfo<T, Tint, Tgroup> info =
        GetPairAutInfo<T, Tint, Tgroup>(can.first, can.second, os);
    MapCanonic[key] = ListPairWork.size();
    ListPairWork.push_back(can);
    ListAutInfo.push_back(info);
    result.ListGram.push_back(can.first);
    result.ListIsom.push_back(can.second.transpose());
    result.ListCentralizerOrder.push_back(info.order);
    result.accumulated_mass +=
        Tmass(1) / UniversalScalarConversion<Tmass, TintGroup>(info.order);
#ifdef DEBUG_GENUS_ENUMERATION_GROUP
    os << "GENUS_GROUP: class " << ListPairWork.size()
       << ", |C(sigma)|=" << info.order << ", mass " << result.accumulated_mass
       << " / " << TotalMass << "\n";
#endif
    return true;
  };
  for (auto &eSeed : ListSeed) {
    MyMatrix<T> const &A = eSeed.first;
    MyMatrix<T> P = eSeed.second.transpose();
    if (P * A * P.transpose() != A) {
      std::cerr << "GENUS_GROUP: a seed isometry does not preserve its Gram "
                   "matrix (g^T A g = A is expected)\n";
      throw TerminalException{1};
    }
    MyMatrix<Tint> P_int = UniversalMatrixConversion<Tint, T>(P);
    int m_seed = FiniteMatrixOrderInt(P_int);
    if (m == -1) {
      m = m_seed;
      det = DeterminantMat(A);
    } else if (m != m_seed || DeterminantMat(A) != det) {
      std::cerr << "GENUS_GROUP: the seeds have different orders or "
                   "determinants\n";
      throw TerminalException{1};
    }
    if (result.accumulated_mass == TotalMass) {
      break;
    }
    (void)f_insert(A, P);
  }
  auto expand_at = [&](size_t pos, int p, size_t i_prime) -> void {
    MyMatrix<T> A = ListPairWork[pos].first;
    MyMatrix<T> P = ListPairWork[pos].second;
    MyMatrix<Tint> G_int = UniversalMatrixConversion<Tint, T>(A);
    MyMatrix<Tint> P_int = UniversalMatrixConversion<Tint, T>(P);
    std::vector<EigenspaceData> ListEig = ComputeEigenspaces(P_int, m, p);
    // The small eigenspaces first: their projective spaces are cheap, and
    // the mass may close before a large one is ever cut into orbits.
    std::stable_sort(ListEig.begin(), ListEig.end(),
                     [](EigenspaceData const &x, EigenspaceData const &y) {
                       return x.space.basis.size() < y.space.basis.size();
                     });
    for (auto &ed : ListEig) {
      size_t k = ed.space.basis.size();
      std::vector<MyMatrix<Tint>> ListGenEig;
      for (auto &h : ListAutInfo[pos].ListGenMat) {
        ListGenEig.push_back(EigenspaceAction(ed, h, p));
      }
      ProjectiveLineSelection ListRepr =
          ProjectiveLineOrbitRepresentatives<Tint>(static_cast<int>(k), p,
                                                   ListGenEig);
#ifdef DEBUG_GENUS_ENUMERATION_GROUP
      os << "GENUS_GROUP: class " << (pos + 1) << " p=" << p
         << " lambda=" << ed.lambda << " dim=" << k << " orbits="
         << ListRepr.size() << "\n";
#endif
      for (size_t i_orbit = 0; i_orbit < ListRepr.size(); i_orbit++) {
        std::vector<int> c = ProjectiveLineOfIndex(static_cast<int>(k), p,
                                                   ListRepr.at(i_orbit));
        std::optional<std::vector<int>> opt_v =
            EigenvectorNeighborVector(ed, c, G_int, p);
        if (!opt_v) {
          continue;
        }
        std::optional<std::pair<MyMatrix<T>, MyMatrix<T>>> opt =
            GenusNeighborPair<T, Tint>(G_int, P, *opt_v, p, os);
        if (!opt) {
          continue;
        }
        if (f_insert(opt->first, opt->second)) {
          result.ListNewClassPerPrime[i_prime]++;
        }
        if (result.accumulated_mass == TotalMass) {
          return;
        }
      }
    }
  };
  if (prime != 0) {
    if (!IsGoodNeighborPrime(prime, m, det)) {
      std::cerr << "GENUS_GROUP: the prime " << prime
                << " divides m det\n";
      throw TerminalException{1};
    }
    result.ListPrimeUsed.push_back(prime);
    result.ListNewClassPerPrime.push_back(0);
    for (size_t pos = 0;
         pos < ListPairWork.size() && result.accumulated_mass != TotalMass;
         pos++) {
      expand_at(pos, prime, 0);
    }
  } else {
    std::vector<size_t> next_expand;
    while (result.accumulated_mass != TotalMass) {
      if (result.accumulated_mass > TotalMass) {
        break;
      }
      size_t i_work = result.ListPrimeUsed.size();
      for (size_t i = 0; i < result.ListPrimeUsed.size(); i++) {
        if (next_expand[i] < ListPairWork.size()) {
          i_work = i;
          break;
        }
      }
      if (i_work == result.ListPrimeUsed.size()) {
        if (result.ListPrimeUsed.size() == max_prime) {
          break;
        }
        // 2 first when it is good, then the primes p = 1 mod m, at which
        // all the eigenvalues are in F_p
        int p_new;
        if (result.ListPrimeUsed.empty() && IsGoodNeighborPrime(2, m, det)) {
          p_new = 2;
        } else {
          int q = result.ListPrimeUsed.empty() ? 2 : result.ListPrimeUsed.back();
          p_new = NextPrimeAfter(q);
          while (!IsGoodNeighborPrime(p_new, m, det) || p_new % m != 1) {
            p_new = NextPrimeAfter(p_new);
          }
        }
        result.ListPrimeUsed.push_back(p_new);
        result.ListNewClassPerPrime.push_back(0);
        next_expand.push_back(0);
        continue;
      }
      expand_at(next_expand[i_work], result.ListPrimeUsed[i_work], i_work);
      next_expand[i_work]++;
    }
  }
  result.complete = (result.accumulated_mass == TotalMass);
#ifdef TIMINGS_GENUS_ENUMERATION_GROUP
  os << "|GENUS_GROUP: GenusEnumerationGroup|=" << time_total << "\n";
#endif
  return result;
}

template <typename T, typename Tgroup>
void WriteGenusEnumerationGroupResult(
    std::ostream &os, GenusEnumerationGroupResult<T, Tgroup> const &result,
    std::string const &OutFormat) {
  if (OutFormat == "Summary") {
    os << "primes used (new classes found):";
    for (size_t i = 0; i < result.ListPrimeUsed.size(); i++) {
      os << " " << result.ListPrimeUsed[i] << " ("
         << result.ListNewClassPerPrime[i] << ")";
    }
    os << "\n";
    os << "classes of pairs = " << result.ListGram.size() << "\n";
    os << "accumulated mass = " << result.accumulated_mass << "\n";
    os << "target mass = " << result.target_mass << "\n";
    os << "complete = " << (result.complete ? "true" : "false") << "\n";
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(ListGram:=[";
    for (size_t i = 0; i < result.ListGram.size(); i++) {
      os << (i > 0 ? ",\n" : "");
      WriteMatrixGAP(os, result.ListGram[i]);
    }
    os << "],\n ListIsom:=[";
    for (size_t i = 0; i < result.ListIsom.size(); i++) {
      os << (i > 0 ? ",\n" : "");
      WriteMatrixGAP(os, result.ListIsom[i]);
    }
    os << "],\n ListCentralizerOrder:=[";
    for (size_t i = 0; i < result.ListCentralizerOrder.size(); i++) {
      os << (i > 0 ? ", " : "") << result.ListCentralizerOrder[i];
    }
    os << "],\n accumulated_mass:=" << result.accumulated_mass
       << ",\n target_mass:=" << result.target_mass
       << ",\n complete:=" << (result.complete ? "true" : "false")
       << ",\n ListPrimeUsed:=[";
    for (size_t i = 0; i < result.ListPrimeUsed.size(); i++) {
      os << (i > 0 ? ", " : "") << result.ListPrimeUsed[i];
    }
    os << "],\n ListNewClassPerPrime:=[";
    for (size_t i = 0; i < result.ListNewClassPerPrime.size(); i++) {
      os << (i > 0 ? ", " : "") << result.ListNewClassPerPrime[i];
    }
    os << "]);\n";
    return;
  }
  std::cerr << "GENUS_GROUP: unknown OutFormat=" << OutFormat
            << " (Summary or GAP)\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_GENERA_GENUSENUMERATIONGROUP_H_
// clang-format on
