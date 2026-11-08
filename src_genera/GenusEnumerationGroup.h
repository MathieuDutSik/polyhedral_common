// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_GENUSENUMERATIONGROUP_H_
#define SRC_GENERA_GENUSENUMERATIONGROUP_H_

// clang-format off
#include "GenusEnumeration.h"
#include "LatticeOrbitSplitCanonical.h"
#include "GroupNeighborModP.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>
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
  space of the eigenspace, visited one orbit at a time.

  The classes. A pair is told apart by the canonical form of the pair of
  forms (A, A P^T) (the second one not symmetric), and C_{O(L)}(sigma) is
  the automorphism group of that pair of forms, both by the
  Plesken-Souvignier engine, whose backtrack the second form prunes. For
  sigma = +-I the centraliser is O(L) and the lattice alone is used.

  Finite abelian groups. For G = <P_1, ..., P_k> commuting, of exponent m,
  the G-stable neighbours are those of the common eigenvectors, the
  eigenvectors of a character chi of G with values in F_p (p = 1 mod m,
  p not dividing m det). The projector on the chi-eigenspace modulo p^2 is
  the product of the projectors of the P_i for the eigenvalues chi(P_i),
  and the isotropy correction is needed only for the real characters. The
  classes are told apart by the forms (A, A P_1^T, ..., A P_k^T), whose
  automorphism group is C_{O(L)}(G); the mass is that of GENUS_IsometryMass
  group.
 */

template <typename T, typename Tgroup> struct GenusEnumerationGroupResult {
  using TintGroup = typename Tgroup::Tint;
  using Tmass = typename overlying_field<TintGroup>::field_type;
  std::vector<MyMatrix<T>> ListGram;
  // the isometries, acting on the columns: g^T A g = A, one list per class
  std::vector<std::vector<MyMatrix<T>>> ListListIsom;
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

inline int64_t InverseModPosInt64(int64_t u, int64_t m) {
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
    int64_t inv = InverseModPosInt64(rows[r][col], p);
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
  // the character takes only the values +-1, so the vectors need not be
  // isotropic
  bool is_real;
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

// sigma = +-I: the centraliser is the whole of O(L)
template <typename T> bool IsScalarIsometry(MyMatrix<T> const &P) {
  int n = P.rows();
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j && P(i, j) != 0) {
        return false;
      }
    }
    if (P(i, i) != P(0, 0)) {
      return false;
    }
  }
  return true;
}

template <typename Tint>
std::vector<EigenspaceData> ComputeEigenspaces(MyMatrix<Tint> const &P, int m,
                                               int p) {
  std::vector<EigenspaceData> res;
  int n = P.rows();
  int64_t p2 = static_cast<int64_t>(p) * p;
  if (IsScalarIsometry(P)) {
    // every vector is an eigenvector, for any p (p = 2 included)
    ModMatrix Id(n, std::vector<int64_t>(n, 0));
    for (int i = 0; i < n; i++) {
      Id[i][i] = 1;
    }
    int64_t c = UniversalScalarConversion<int64_t, Tint>(P(0, 0));
    RowSpaceModP space = ComputeRowSpaceModP(Id, p);
    res.push_back({ModPosInt64(c, p), ModPosInt64(c, p2), std::move(Id),
                   std::move(space), true});
    return res;
  }
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
    int64_t inv = InverseModPosInt64(Rlambda, p2);
    for (auto &row : RP) {
      for (auto &x : row) {
        x = (x * inv) % p2;
      }
    }
    RowSpaceModP space = ComputeRowSpaceModP(RP, p);
    if (space.basis.empty()) {
      continue;
    }
    bool is_real = (lambda == 1 || lambda == p - 1);
    res.push_back({lambda, lambda2, std::move(RP), std::move(space), is_real});
  }
  return res;
}

/*
  The joint eigenspaces of commuting P_i, all of order dividing m: the
  products of their projectors modulo p^2, which commute and are
  idempotent, so the product is the projector on the intersection.
 */
template <typename Tint>
std::vector<EigenspaceData>
ComputeJointEigenspaces(std::vector<MyMatrix<Tint>> const &ListP, int m,
                        int p) {
  std::vector<EigenspaceData> res = ComputeEigenspaces(ListP[0], m, p);
  int64_t p2 = static_cast<int64_t>(p) * p;
  for (size_t i = 1; i < ListP.size(); i++) {
    if (IsScalarIsometry(ListP[i])) {
      continue;
    }
    std::vector<EigenspaceData> ListEig = ComputeEigenspaces(ListP[i], m, p);
    std::vector<EigenspaceData> next;
    for (auto const &ed1 : res) {
      for (auto const &ed2 : ListEig) {
        ModMatrix E = ModMatrixProduct(ed1.E, ed2.E, p2);
        RowSpaceModP space = ComputeRowSpaceModP(E, p);
        if (space.basis.empty()) {
          continue;
        }
        next.push_back({ed1.lambda, ed1.lambda2, std::move(E),
                        std::move(space), ed1.is_real && ed2.is_real});
      }
    }
    res = std::move(next);
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
  bool is_pm1 = ed.is_real;
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
                 InverseModPosInt64(ResidueModP(twob_red, p), p)) %
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
  The orbits of the projective space of dimension n - 1 over F_p, visited
  one at a time: f receives the index of the representative of an orbit as
  soon as that orbit is closed and returns true to stop. Unlike
  ProjectiveLineOrbitRepresentatives, the enumeration of pairs usually
  closes its mass after a few orbits, and the rest of the space is never
  traversed.
 */
template <typename Tint, typename F>
void ProjectiveLineOrbitForEach(int n, int p,
                                std::vector<MyMatrix<Tint>> const &ListGen,
                                F f) {
  size_t n_line = NumberProjectiveLines(n, p);
  std::vector<std::vector<int>> ListGenRed;
  for (auto &eGen : ListGen) {
    if (!IsScalarModP(eGen, p)) {
      std::vector<int> gred(n * n);
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          gred[i * n + j] = ResidueModP(eGen(i, j), p);
        }
      }
      ListGenRed.push_back(gred);
    }
  }
  if (ListGenRed.empty()) {
    for (size_t i_line = 0; i_line < n_line; i_line++) {
      if (f(i_line)) {
        return;
      }
    }
    return;
  }
  std::vector<int> w(n);
  auto image = [&](std::vector<int> const &v,
                   std::vector<int> const &gred) -> size_t {
    for (int j = 0; j < n; j++) {
      int64_t sum = 0;
      for (int i = 0; i < n; i++) {
        sum += static_cast<int64_t>(v[i]) * gred[i * n + j];
      }
      w[j] = static_cast<int>(sum % p);
    }
    return IndexOfProjectiveLine(w, p);
  };
  Face status(n_line);
  std::vector<size_t> ListActive;
  for (size_t i_line = 0; i_line < n_line; i_line++) {
    if (status[i_line]) {
      continue;
    }
    status[i_line] = 1;
    ListActive.push_back(i_line);
    while (!ListActive.empty()) {
      size_t pos = ListActive.back();
      ListActive.pop_back();
      std::vector<int> v = ProjectiveLineOfIndex(n, p, pos);
      for (auto &gred : ListGenRed) {
        size_t img = image(v, gred);
        if (!status[img]) {
          status[img] = 1;
          ListActive.push_back(img);
        }
      }
    }
    if (f(i_line)) {
      return;
    }
  }
}

/*
  A neighbour pair: the reduced Gram matrix of N and the isometry carried
  over to its reduced basis (on the row vectors).
 */
template <typename T, typename Tint>
std::optional<std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>>
GenusNeighborPair(MyMatrix<Tint> const &G, std::vector<MyMatrix<T>> const &ListP,
                  std::vector<int> const &v, int p, std::ostream &os) {
  std::optional<GenusNeighborBasis<T, Tint>> opt =
      GenusNeighborKernel<T, Tint>(G, v, p, os);
  if (!opt) {
    return {};
  }
  LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(opt->GramN, os);
  // the reduced basis of N in the coordinates of L, times p
  MyMatrix<T> Bred = UniversalMatrixConversion<T, Tint>(rec.Pmat * opt->Bp);
  MyMatrix<T> BredInv = Inverse(Bred);
  std::vector<MyMatrix<T>> ListS;
  for (auto const &P : ListP) {
    MyMatrix<T> S = Bred * P * BredInv;
    if (!IsIntegralMatrix(S)) {
      std::cerr << "GENUS_GROUP: the neighbour is not stable under the "
                   "isometry\n";
      throw TerminalException{1};
    }
    ListS.push_back(S);
  }
  return std::make_pair(rec.GramMatRed, ListS);
}

template <typename T, typename Tint, typename Tgroup> struct PairAutInfo {
  typename Tgroup::Tint order;
  std::vector<MyMatrix<Tint>> ListGenMat;
};

template <typename T>
std::vector<MyMatrix<T>> PairForms(MyMatrix<T> const &A,
                                   std::vector<MyMatrix<T>> const &ListP) {
  std::vector<MyMatrix<T>> ListMat{A};
  for (auto const &P : ListP) {
    ListMat.push_back(A * P.transpose());
  }
  return ListMat;
}

template <typename T>
bool AllScalarIsometries(std::vector<MyMatrix<T>> const &ListP) {
  for (auto const &P : ListP) {
    if (!IsScalarIsometry(P)) {
      return false;
    }
  }
  return true;
}

template <typename T, typename Tint, typename Tgroup>
PairAutInfo<T, Tint, Tgroup>
GetPairAutInfo(MyMatrix<T> const &A, std::vector<MyMatrix<T>> const &ListP,
               std::ostream &os) {
  if (AllScalarIsometries(ListP)) {
    // The pair forms are A and +-A, whose common automorphisms are those of
    // A alone; the methods of the genus enumeration avoid a
    // Plesken-Souvignier search over the whole of O(L), which is slow on
    // root lattices such as E8^2.
    LatticeAutInfo<T, Tint, Tgroup> lat =
        GetLatticeAutInfo<T, Tint, Tgroup>(A, GenusMethods(), os);
    return {lat.order, lat.ListGenMat};
  }
  std::vector<MyMatrix<T>> ListMat = PairForms(A, ListP);
  // Plesken-Souvignier with the second form: the backtrack is pruned by
  // both forms, which keeps it small even on lattices with huge
  // automorphism groups, where the graph of the characteristic vector
  // family of the first form alone explodes.
  PleskenSouvignierAutomResult<Tint, Tgroup> aut =
      PleskenSouvignierLatticeAutomorphism<T, Tint, Tgroup>(ListMat, os);
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

/*
  The canonical representative of the pair: the reduced forms and the
  isometry in the canonical basis, with the method used (0: orbit splitting
  over the Plesken-Souvignier automorphism group, 1: the graph of the
  characteristic vector family, 2: sigma = +-I, the canonical form of the
  lattice alone). The orbit splitting may decline an input;
  that depends only on the configuration, so isometric pairs take the same
  branch, and the method is part of the key since the two canonical forms
  are not comparable.
 */
template <typename T> struct CanonicalPairData {
  int method;
  MyMatrix<T> Acan;
  std::vector<MyMatrix<T>> ListPcan;
};

template <typename T, typename Tint, typename Tgroup>
CanonicalPairData<T> CanonicalPair(MyMatrix<T> const &A,
                                   std::vector<MyMatrix<T>> const &ListP,
                                   std::ostream &os) {
  if (AllScalarIsometries(ListP)) {
    // the P are fixed by every change of basis; the key is that of A alone
    return {2,
            GenusCanonicalGram<T, Tint, Tgroup>(A, GenusMethods().canonical,
                                                os),
            ListP};
  }
  std::vector<MyMatrix<T>> ListMat = PairForms(A, ListP);
  std::optional<MyMatrix<Tint>> opt =
      OrbitSplitCanonicalizeListMat<T, Tint, Tgroup>(ListMat, os);
  int method = 0;
  MyMatrix<Tint> B;
  if (opt) {
    B = *opt;
  } else {
    method = 1;
    B = ComputeCanonicalFormMultiple<T, Tint, Tgroup>(ListMat, os);
  }
  MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
  MyMatrix<T> Acan = B_T * A * B_T.transpose();
  MyMatrix<T> B_inv = Inverse(B_T);
  std::vector<MyMatrix<T>> ListPcan;
  for (auto const &P : ListP) {
    ListPcan.push_back(B_T * P * B_inv);
  }
  return {method, Acan, ListPcan};
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
  space over F_2 is the cheapest of all. For sigma = +-I every sublattice
  is stable and only det matters.
 */
template <typename T>
bool IsGoodNeighborPrime(int p, int m, T const &det, bool scalar) {
  // for sigma = +-I, the ordinary neighbours at any p not dividing det
  if (!scalar && m % p == 0) {
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
  The good primes in the order of the cost of their neighbours: the number
  of points sum_chi (p^k_chi - 1) / (p - 1) of the projective spaces of the
  joint eigenspaces defined over F_p, of dimensions k_chi. At p != 1 mod m
  only the characters with values in F_p contribute; their neighbours may
  already reach every class, and far more cheaply: for E8+E8 with a group
  of order 6 whose fixed lattice has rank 13, p = 5 visits P^12(F_5) where
  p = 7 visits P^12(F_7), 50 times larger. A prime with no eigenspace over
  F_p gives no neighbour and is left out. The genus fixes the dimensions,
  so the first class gives them.
 */
template <typename T, typename Tint = mpz_class>
std::vector<int> OrderedNeighborPrimes(std::vector<MyMatrix<T>> const &ListP,
                                       int m, T const &det, bool scalar) {
  std::vector<MyMatrix<Tint>> ListP_int;
  for (auto const &P : ListP) {
    ListP_int.push_back(UniversalMatrixConversion<Tint, T>(P));
  }
  // the first good primes with an eigenspace over F_p
  size_t n_cand = 12;
  std::vector<std::pair<double, int>> ListCost;
  int p = 1;
  while (ListCost.size() < n_cand) {
    p = NextPrimeAfter(p);
    if (!IsGoodNeighborPrime(p, m, det, scalar)) {
      continue;
    }
    std::vector<EigenspaceData> ListEig =
        ComputeJointEigenspaces(ListP_int, m, p);
    if (ListEig.empty()) {
      continue;
    }
    double cost = 0;
    for (auto const &ed : ListEig) {
      double k = static_cast<double>(ed.space.basis.size());
      cost += (std::pow(static_cast<double>(p), k) - 1) / (p - 1);
    }
    ListCost.push_back({cost, p});
  }
  std::stable_sort(ListCost.begin(), ListCost.end());
  std::vector<int> res;
  for (auto const &pair : ListCost) {
    res.push_back(pair.second);
  }
  return res;
}

/*
  The elements of the finite group generated by the ListP (on the rows),
  with their inverses, by closure.
 */
template <typename T>
std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>>
GroupElementsAndInverses(std::vector<MyMatrix<T>> const &ListP,
                         size_t max_group) {
  int n = ListP[0].rows();
  auto flatten = [&](MyMatrix<T> const &M) -> std::vector<T> {
    std::vector<T> v;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        v.push_back(M(i, j));
      }
    }
    return v;
  };
  std::vector<MyMatrix<T>> elts{IdentityMat<T>(n)};
  std::set<std::vector<T>> seen{flatten(elts[0])};
  for (size_t pos = 0; pos < elts.size(); pos++) {
    for (auto const &h : ListP) {
      MyMatrix<T> x = elts[pos] * h;
      if (seen.insert(flatten(x)).second) {
        elts.push_back(x);
        if (elts.size() > max_group) {
          std::cerr << "GENUS_GROUP: the group has more than " << max_group
                    << " elements\n";
          throw TerminalException{1};
        }
      }
    }
  }
  std::vector<MyMatrix<T>> inv;
  for (auto const &x : elts) {
    inv.push_back(Inverse(x));
  }
  return {elts, inv};
}

template <typename T>
ModMat ModMatOfMatrix(MyMatrix<T> const &M, int64_t mod) {
  int n = M.rows(), l = M.cols();
  ModMat R(n, std::vector<int64_t>(l));
  T mod_T = UniversalScalarConversion<T, int64_t>(mod);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < l; j++) {
      T r = ResInt(M(i, j), mod_T);
      R[i][j] = ModPosInt64(UniversalScalarConversion<int64_t, T>(r), mod);
    }
  }
  return R;
}

/*
  The G-module data of a class at the prime p: the irreducible types of
  L/pL and the elements of G modulo p and p^2.
 */
struct GroupModuleAtPrime {
  std::vector<IrreducibleType> types;
  std::vector<ModMat> gens;
  std::vector<ModMat> elt2;
  std::vector<ModMat> eltinv2;
  ModMat A;
};

template <typename T>
GroupModuleAtPrime ComputeGroupModuleAtPrime(MyMatrix<T> const &A,
                                             std::vector<MyMatrix<T>> const &ListP,
                                             int64_t p) {
  GroupModuleAtPrime gm;
  std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>> ei =
      GroupElementsAndInverses(ListP, 1000000);
  std::vector<ModMat> elt, eltinv;
  for (size_t i = 0; i < ei.first.size(); i++) {
    elt.push_back(ModMatOfMatrix(ei.first[i], p));
    eltinv.push_back(ModMatOfMatrix(ei.second[i], p));
    gm.elt2.push_back(ModMatOfMatrix(ei.first[i], p * p));
    gm.eltinv2.push_back(ModMatOfMatrix(ei.second[i], p * p));
  }
  for (auto const &P : ListP) {
    gm.gens.push_back(ModMatOfMatrix(P, p));
  }
  gm.A = ModMatOfMatrix(A, p * p);
  ModMat A_p = ModMatOfMatrix(A, p);
  // a fixed seed: the decomposition is deterministic for a given input
  std::mt19937_64 rng(static_cast<uint64_t>(p) * 1000003ULL + A.rows());
  gm.types = IrreducibleTypes(static_cast<size_t>(A.rows()), gm.gens, elt,
                              eltinv, A_p, p, rng);
  return gm;
}

inline double TypeCost(IrreducibleType const &tp, int64_t p) {
  double q = std::pow(static_cast<double>(p), tp.f);
  double m = static_cast<double>(tp.Phi.size());
  return (std::pow(q, m) - 1) / (q - 1);
}

/*
  The odd good primes for a non-abelian group, in the order of the number
  of points sum_tau (q^m_tau - 1) / (q - 1) of their projective spaces.
 */
template <typename T>
std::vector<int> OrderedNeighborPrimesModule(MyMatrix<T> const &A,
                                             std::vector<MyMatrix<T>> const &ListP,
                                             int m, T const &det) {
  std::vector<std::pair<double, int>> ListCost;
  int p = 2;
  while (ListCost.size() < 8) {
    p = NextPrimeAfter(p);
    if (!IsGoodNeighborPrime(p, m, det, false)) {
      continue;
    }
    GroupModuleAtPrime gm = ComputeGroupModuleAtPrime(A, ListP, p);
    double cost = 0;
    for (auto const &tp : gm.types) {
      cost += TypeCost(tp, p);
    }
    ListCost.push_back({cost, p});
  }
  std::stable_sort(ListCost.begin(), ListCost.end());
  std::vector<int> res;
  for (auto const &pair : ListCost) {
    res.push_back(pair.second);
  }
  return res;
}

/*
  The orbits of a finite set of n_point indices under maps image(i_gen,
  idx), visited one at a time as in ProjectiveLineOrbitForEach.
 */
template <typename Fimage, typename F>
void IndexOrbitForEach(size_t n_point, size_t n_gen, Fimage image, F f) {
  Face status(n_point);
  std::vector<size_t> ListActive;
  for (size_t i = 0; i < n_point; i++) {
    if (status[i]) {
      continue;
    }
    status[i] = 1;
    ListActive.push_back(i);
    while (!ListActive.empty()) {
      size_t pos = ListActive.back();
      ListActive.pop_back();
      for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
        size_t img = image(i_gen, pos);
        if (!status[img]) {
          status[img] = 1;
          ListActive.push_back(img);
        }
      }
    }
    if (f(i)) {
      return;
    }
  }
}

/*
  The neighbour of the lattice with Gram matrix A for the lifted submodule
  St: the reduced Gram matrix and the isometries in its basis.
 */
template <typename T, typename Tint>
std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>
GroupNeighborPair(MyMatrix<T> const &A, std::vector<MyMatrix<T>> const &ListP,
                  ModMat const &St, ModMat const &A2, int64_t p,
                  std::ostream &os) {
  int n = A.rows();
  ModMat gens = GroupNeighborGenerators(St, A2, p);
  MyMatrix<Tint> Gen(static_cast<int>(gens.size()), n);
  for (size_t i = 0; i < gens.size(); i++) {
    for (int j = 0; j < n; j++) {
      Gen(static_cast<int>(i), j) =
          UniversalScalarConversion<Tint, int64_t>(gens[i][j]);
    }
  }
  MyMatrix<Tint> Bz = GetZbasis(Gen);
  MyMatrix<T> B = UniversalMatrixConversion<T, Tint>(Bz) /
                  UniversalScalarConversion<T, int64_t>(p);
  MyMatrix<T> GramN = B * A * B.transpose();
  if (!IsIntegralMatrix(GramN)) {
    std::cerr << "GENUS_GROUP: the module neighbour is not integral\n";
    throw TerminalException{1};
  }
  for (int i = 0; i < n; i++) {
    if (!IsInteger(GramN(i, i) / T(2))) {
      std::cerr << "GENUS_GROUP: the module neighbour is not even\n";
      throw TerminalException{1};
    }
  }
  LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(GramN, os);
  MyMatrix<T> Bred = UniversalMatrixConversion<T, Tint>(rec.Pmat) * B;
  MyMatrix<T> BredInv = Inverse(Bred);
  std::vector<MyMatrix<T>> ListS;
  for (auto const &P : ListP) {
    MyMatrix<T> S = Bred * P * BredInv;
    if (!IsIntegralMatrix(S)) {
      std::cerr << "GENUS_GROUP: the module neighbour is not G-stable\n";
      throw TerminalException{1};
    }
    ListS.push_back(S);
  }
  return {rec.GramMatRed, ListS};
}

/*
  The enumeration. ListSeed holds pairs (Gram, isometries on the columns),
  the isometries commuting and generating a finite abelian group G of
  exponent m. With prime != 0 only that prime is used (any good prime).
  With prime = 0 the good primes are brought in one after the other in the
  order of OrderedNeighborPrimes, as in the iterate strategy of
  GenusEnumeration, at most max_prime of them.
 */
template <typename T, typename Tint, typename Tgroup>
GenusEnumerationGroupResult<T, Tgroup> GenusEnumerationAbelianGroup(
    std::vector<std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>> const
        &ListSeed,
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
  // the work list: the canonical pairs (Gram, isometries on the rows)
  std::vector<std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>> ListPairWork;
  std::vector<PairAutInfo<T, Tint, Tgroup>> ListAutInfo;
  std::map<std::pair<int, std::vector<MyMatrix<T>>>, size_t> MapCanonic;
  // m: the exponent of the group
  int m = -1;
  size_t n_gen = 0;
  T det(0);
  bool scalar = false;
  // commuting generators: the eigenvectors; otherwise the G-submodules
  bool abelian = true;
  /*
    The projective spaces are visited in rounds. Building the orbits costs
    the size of the space, as the first orbit covers most of it
    (P^13(F_5) for the fixed lattice of rank 14 of W(A_2) on E8+E8 has 1.5
    10^9 points), while random points find the frequent classes at once
    but may take long to hit a rare one. So: up to small_space points the
    orbits are built; up to large_space points random points come first,
    and the orbits are built before a further prime is brought in, if the
    mass has not closed (the cheapest prime is exhausted first); beyond,
    only random points, four times as many in each final round. The seeds
    are fixed per round. The mass, not the visit, certifies completeness.
   */
  size_t const small_space = 5000000;
  size_t const large_space = 20000000;
  size_t n_sample = 200;
  size_t i_round = 0;
  bool more_rounds = false;
  // the medium spaces whose orbits were deferred, by prime index
  std::vector<uint8_t> deferred;
  bool orbit_phase = false;
  size_t i_prime_cur = 0;
  // 0: nothing to do, 1: build the orbits, 2: random points
  auto visit_mode = [&](size_t n_point) -> int {
    if (orbit_phase) {
      return (n_point > small_space && n_point <= large_space) ? 1 : 0;
    }
    if (n_point <= small_space) {
      return (i_round == 0) ? 1 : 0;
    }
    if (n_point <= large_space) {
      if (i_round == 0) {
        if (deferred.size() <= i_prime_cur) {
          deferred.resize(i_prime_cur + 1, 0);
        }
        deferred[i_prime_cur] = 1;
        return 2;
      }
      return 0;
    }
    more_rounds = true;
    return 2;
  };
  /*
    A class is recognised by isometry tests rather than canonical forms:
    the classes are bucketed by the canonical form of the lattice alone
    (the fast Weyl-chamber one), and a neighbour is tested against those of
    its bucket by Plesken-Souvignier on the pair forms, pruned by the
    known centraliser of the class. The canonical form of a pair is slow
    when the centraliser is large (10^12 for W(A3) on E8+E8, 99% of the
    time), and with random points most neighbours are known classes, for
    which an isometry is found quickly.
   */
  /*
    A known class is recognised in one of two ways, chosen by their
    measured cost: the canonical form of the pair (key in MapCanonic), or
    isometry tests. For the latter the classes with the same lattice are
    put in one bucket, keyed by the canonical form Akey = B A B^T of the
    lattice (the fast Weyl-chamber one); every pair of the bucket is moved
    by B to that canonical basis, then by a fixed LLL reduction R of Akey,
    so that all of them share one Gram matrix and one Z-spanning family,
    computed once per bucket, and a neighbour costs the canonical form of
    its lattice and the Plesken-Souvignier backtrack against the classes of
    its bucket, pruned by their centralisers. Which is cheaper depends on
    the case: the canonical form of the pair takes 99% of the time for
    W(A3) on E8+E8 and less than the isometry tests for two reflections.
    The first neighbours try both in turn, then the cheaper on average is
    used; a class gets the data of a method when that method first meets
    it, so both always see all the classes.
   */
  struct PairBucket {
    MyMatrix<Tint> Rmat;
    MyMatrix<Tint> SHVhalf;
    std::vector<size_t> classes;
  };
  std::map<MyMatrix<T>, PairBucket> MapBucket;
  std::vector<std::vector<MyMatrix<Tint>>> ListFormsRed;
  std::vector<std::vector<MyMatrix<Tint>>> ListAutRed;
  std::vector<uint8_t> ListInBucket;
  std::vector<uint8_t> ListHasKey;
  double t_can = 0, t_iso = 0;
  size_t n_can = 0, n_iso = 0;
  using CanonKey = std::pair<int, std::vector<MyMatrix<T>>>;
  auto canon_key = [&](MyMatrix<T> const &A,
                       std::vector<MyMatrix<T>> const &ListP) -> CanonKey {
    CanonicalPairData<T> can = CanonicalPair<T, Tint, Tgroup>(A, ListP, os);
    return {can.method, PairForms(can.Acan, can.ListPcan)};
  };
  // the bucket data of a pair: its bucket and its forms in the reduced
  // canonical basis, with the basis change M = R B
  struct BucketData {
    MyMatrix<T> Akey;
    MyMatrix<Tint> M;
    std::vector<MyMatrix<Tint>> forms;
  };
  auto bucket_data = [&](MyMatrix<T> const &A,
                         std::vector<MyMatrix<T>> const &ListP) -> BucketData {
    MyMatrix<Tint> B =
        LatticeCanonicalBasis<T, Tint, Tgroup>(A, GenusMethods().canonical, os);
    MyMatrix<T> B_T = UniversalMatrixConversion<T, Tint>(B);
    MyMatrix<T> Akey = B_T * A * B_T.transpose();
    auto iter = MapBucket.find(Akey);
    if (iter == MapBucket.end()) {
      PairBucket bucket;
      LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(Akey, os);
      bucket.Rmat = rec.Pmat;
      bucket.SHVhalf =
          PleskenSouvignierZspanFamilyHalf<T, Tint>(rec.GramMatRed, false, os);
      iter = MapBucket.emplace(Akey, bucket).first;
    }
    MyMatrix<Tint> M = iter->second.Rmat * B;
    MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(M);
    std::vector<MyMatrix<Tint>> forms;
    for (auto const &F : PairForms(A, ListP)) {
      forms.push_back(UniversalMatrixConversion<Tint, T>(
          MyMatrix<T>(M_T * F * M_T.transpose())));
    }
    return {Akey, M, forms};
  };
  auto put_in_bucket = [&](size_t idx, BucketData const &bd) -> void {
    MyMatrix<Tint> Minv = Inverse(bd.M);
    std::vector<MyMatrix<Tint>> gens_red;
    for (auto const &h : ListAutInfo[idx].ListGenMat) {
      gens_red.push_back(bd.M * h * Minv);
    }
    ListFormsRed[idx] = bd.forms;
    ListAutRed[idx] = gens_red;
    ListInBucket[idx] = 1;
    MapBucket[bd.Akey].classes.push_back(idx);
  };
  auto add_class = [&](MyMatrix<T> const &A,
                       std::vector<MyMatrix<T>> const &ListP) -> size_t {
    PairAutInfo<T, Tint, Tgroup> info =
        GetPairAutInfo<T, Tint, Tgroup>(A, ListP, os);
    size_t idx = ListPairWork.size();
    ListPairWork.push_back({A, ListP});
    ListAutInfo.push_back(info);
    ListFormsRed.push_back({});
    ListAutRed.push_back({});
    ListInBucket.push_back(0);
    ListHasKey.push_back(0);
    result.ListGram.push_back(A);
    std::vector<MyMatrix<T>> ListIsom;
    for (auto const &P : ListP) {
      ListIsom.push_back(P.transpose());
    }
    result.ListListIsom.push_back(ListIsom);
    result.ListCentralizerOrder.push_back(info.order);
    result.accumulated_mass +=
        Tmass(1) / UniversalScalarConversion<Tmass, TintGroup>(info.order);
#ifdef DEBUG_GENUS_ENUMERATION_GROUP
    os << "GENUS_GROUP: class " << ListPairWork.size()
       << ", |C(sigma)|=" << info.order << ", mass " << result.accumulated_mass
       << " / " << TotalMass << "\n";
#endif
    return idx;
  };
  auto f_insert_canonical = [&](MyMatrix<T> const &A,
                                std::vector<MyMatrix<T>> const &ListP) -> bool {
    for (size_t idx = 0; idx < ListPairWork.size(); idx++) {
      if (!ListHasKey[idx]) {
        MapCanonic[canon_key(ListPairWork[idx].first,
                             ListPairWork[idx].second)] = idx;
        ListHasKey[idx] = 1;
      }
    }
    CanonKey key = canon_key(A, ListP);
    if (MapCanonic.count(key) == 1) {
      return false;
    }
    size_t idx = add_class(A, ListP);
    MapCanonic[key] = idx;
    ListHasKey[idx] = 1;
    return true;
  };
  auto f_insert_isometry = [&](MyMatrix<T> const &A,
                               std::vector<MyMatrix<T>> const &ListP) -> bool {
    for (size_t idx = 0; idx < ListPairWork.size(); idx++) {
      if (!ListInBucket[idx]) {
        put_in_bucket(idx, bucket_data(ListPairWork[idx].first,
                                       ListPairWork[idx].second));
      }
    }
    BucketData bd = bucket_data(A, ListP);
    PairBucket &bucket = MapBucket[bd.Akey];
    for (size_t i_b = 0; i_b < bucket.classes.size(); i_b++) {
      size_t idx = bucket.classes[i_b];
      std::optional<MyMatrix<Tint>> opt = PleskenSouvignierIsometry<Tint>(
          bd.forms, bucket.SHVhalf, ListFormsRed[idx], bucket.SHVhalf,
          ListAutRed[idx], PSIsometryOptions<Tint>{.depth = 0}, os);
      if (opt) {
        // the class matched last is tried first next time
        std::swap(bucket.classes[i_b], bucket.classes[0]);
        return false;
      }
    }
    size_t idx = add_class(A, ListP);
    put_in_bucket(idx, bd);
    return true;
  };
  auto f_insert = [&](MyMatrix<T> const &A,
                      std::vector<MyMatrix<T>> const &ListP) -> bool {
    bool use_canonical;
    if (AllScalarIsometries(ListP)) {
      use_canonical = true;
    } else if (n_can + n_iso < 6) {
      use_canonical = (n_can <= n_iso);
    } else {
      use_canonical = (t_can / static_cast<double>(n_can) <=
                       t_iso / static_cast<double>(n_iso));
    }
    auto start = std::chrono::steady_clock::now();
    bool res = use_canonical ? f_insert_canonical(A, ListP)
                             : f_insert_isometry(A, ListP);
    double dt = std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - start)
                    .count();
    if (use_canonical) {
      t_can += dt;
      n_can++;
    } else {
      t_iso += dt;
      n_iso++;
    }
    return res;
  };
  for (auto &eSeed : ListSeed) {
    MyMatrix<T> const &A = eSeed.first;
    std::vector<MyMatrix<T>> ListP;
    int m_seed = 1;
    bool commuting = true;
    for (auto const &g : eSeed.second) {
      MyMatrix<T> P = g.transpose();
      if (P * A * P.transpose() != A) {
        std::cerr << "GENUS_GROUP: a seed isometry does not preserve its "
                     "Gram matrix (g^T A g = A is expected)\n";
        throw TerminalException{1};
      }
      for (auto const &Q : ListP) {
        if (P * Q != Q * P) {
          commuting = false;
        }
      }
      MyMatrix<Tint> P_int = UniversalMatrixConversion<Tint, T>(P);
      m_seed = std::lcm(m_seed, FiniteMatrixOrderInt(P_int));
      ListP.push_back(P);
    }
    if (m == -1) {
      if (ListP.empty()) {
        std::cerr << "GENUS_GROUP: a seed has no isometry\n";
        throw TerminalException{1};
      }
      abelian = commuting;
      if (!abelian) {
        // the exponent from all the elements: the orders of the
        // generators do not see the primes of |G|
        std::vector<MyMatrix<T>> elts =
            GroupElementsAndInverses(ListP, 1000000).first;
        for (auto const &x : elts) {
          MyMatrix<Tint> x_int = UniversalMatrixConversion<Tint, T>(x);
          m_seed = std::lcm(m_seed, FiniteMatrixOrderInt(x_int));
        }
      }
      m = m_seed;
      n_gen = ListP.size();
      det = DeterminantMat(A);
      scalar = AllScalarIsometries(ListP);
    } else if ((abelian && m != m_seed) || commuting != abelian ||
               ListP.size() != n_gen || DeterminantMat(A) != det) {
      std::cerr << "GENUS_GROUP: the seeds have different exponents, numbers "
                   "of isometries or determinants\n";
      throw TerminalException{1};
    }
    if (result.accumulated_mass == TotalMass) {
      break;
    }
    (void)f_insert(A, ListP);
  }
  auto expand_at_module = [&](size_t pos, int p, size_t i_prime) -> void {
    MyMatrix<T> A = ListPairWork[pos].first;
    std::vector<MyMatrix<T>> ListP = ListPairWork[pos].second;
    GroupModuleAtPrime gm = ComputeGroupModuleAtPrime(A, ListP, p);
    std::vector<size_t> order(gm.types.size());
    for (size_t i = 0; i < order.size(); i++) {
      order[i] = i;
    }
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
      return TypeCost(gm.types[a], p) < TypeCost(gm.types[b], p);
    });
    std::vector<ModMat> ListCent;
    for (auto const &h : ListAutInfo[pos].ListGenMat) {
      ListCent.push_back(ModMatOfMatrix(UniversalMatrixConversion<T, Tint>(h), p));
    }
    for (size_t i_type : order) {
      IrreducibleType const &tp = gm.types[i_type];
      FiniteFieldModP F{p, tp.f, tp.mu, 1};
      for (int a = 0; a < tp.f; a++) {
        F.q *= p;
      }
      int mdim = static_cast<int>(tp.Phi.size());
      std::vector<ModMat> theta_pow = ThetaPowers(tp, p);
      std::vector<std::vector<std::vector<std::vector<int64_t>>>> ListAct;
      for (auto const &h : ListCent) {
        ListAct.push_back(ActionOnHom(tp, theta_pow, h, p));
      }
      size_t n_point = NumberProjectivePointsFq(mdim, F.q);
      auto image = [&](size_t i_gen, size_t idx) -> size_t {
        std::vector<std::vector<int64_t>> alpha =
            ProjectivePointFqOfIndex(F, mdim, idx);
        return IndexOfProjectivePointFq(
            F, mdim, ApplyActionOnHom(F, ListAct[i_gen], alpha));
      };
#ifdef DEBUG_GENUS_ENUMERATION_GROUP
      size_t n_orbit = 0;
#endif
      bool closed = false;
      auto visit = [&](size_t idx) -> bool {
#ifdef DEBUG_GENUS_ENUMERATION_GROUP
        n_orbit++;
#endif
        std::vector<std::vector<int64_t>> alpha =
            ProjectivePointFqOfIndex(F, mdim, idx);
        ModMat S_rows = HomOfPoint(tp, theta_pow, alpha, p);
        if (!tp.isotypic_isotropic) {
          ModMat A_p = ModMatOfMatrix(A, p);
          ModMat SA = ModMatMul(S_rows, A_p, p);
          for (auto const &x : SA) {
            for (auto const &y : S_rows) {
              int64_t v = 0;
              for (size_t c = 0; c < x.size(); c++) {
                v = (v + MulModQ(x[c], y[c], p)) % p;
              }
              if (v != 0) {
                return false;
              }
            }
          }
        }
        std::optional<ModMat> St =
            LiftIsotropicSubmodule(S_rows, gm.A, gm.elt2, gm.eltinv2, p);
        if (!St) {
          return false;
        }
        std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>> nb =
            GroupNeighborPair<T, Tint>(A, ListP, *St, gm.A, p, os);
        if (f_insert(nb.first, nb.second)) {
          result.ListNewClassPerPrime[i_prime]++;
        }
        closed = (result.accumulated_mass == TotalMass);
        return closed;
      };
      int mode = visit_mode(n_point);
      if (mode == 1) {
        IndexOrbitForEach(n_point, ListAct.size(), image, visit);
      } else if (mode == 2) {
        std::mt19937_64 rng_point(static_cast<uint64_t>(pos) * 7919ULL +
                                  static_cast<uint64_t>(p) * 104729ULL +
                                  i_round);
        std::uniform_int_distribution<size_t> unif(0, n_point - 1);
        for (size_t i_s = 0; i_s < n_sample; i_s++) {
          if (visit(unif(rng_point))) {
            break;
          }
        }
      }
#ifdef DEBUG_GENUS_ENUMERATION_GROUP
      os << "GENUS_GROUP: class " << (pos + 1) << " p=" << p << " type dim="
         << tp.S0.space.basis.size() << " q=" << F.q << " m=" << mdim
         << " orbits visited=" << n_orbit << "\n";
#endif
      if (closed) {
        return;
      }
    }
  };
  auto expand_at = [&](size_t pos, int p, size_t i_prime) -> void {
    i_prime_cur = i_prime;
    if (!abelian) {
      expand_at_module(pos, p, i_prime);
      return;
    }
    MyMatrix<T> A = ListPairWork[pos].first;
    std::vector<MyMatrix<T>> ListP = ListPairWork[pos].second;
    MyMatrix<Tint> G_int = UniversalMatrixConversion<Tint, T>(A);
    std::vector<MyMatrix<Tint>> ListP_int;
    for (auto const &P : ListP) {
      ListP_int.push_back(UniversalMatrixConversion<Tint, T>(P));
    }
    std::vector<EigenspaceData> ListEig =
        ComputeJointEigenspaces(ListP_int, m, p);
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
#ifdef DEBUG_GENUS_ENUMERATION_GROUP
      size_t n_orbit = 0;
#endif
      bool closed = false;
      auto visit_line = [&](size_t i_line) -> bool {
#ifdef DEBUG_GENUS_ENUMERATION_GROUP
            n_orbit++;
#endif
            std::vector<int> c =
                ProjectiveLineOfIndex(static_cast<int>(k), p, i_line);
            std::optional<std::vector<int>> opt_v =
                EigenvectorNeighborVector(ed, c, G_int, p);
            if (!opt_v) {
              return false;
            }
            std::optional<std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>>
                opt = GenusNeighborPair<T, Tint>(G_int, ListP, *opt_v, p, os);
            if (!opt) {
              return false;
            }
            if (f_insert(opt->first, opt->second)) {
              result.ListNewClassPerPrime[i_prime]++;
            }
            closed = (result.accumulated_mass == TotalMass);
            return closed;
          };
      size_t n_line = NumberProjectiveLines(static_cast<int>(k), p);
      int mode = visit_mode(n_line);
      if (mode == 1) {
        ProjectiveLineOrbitForEach<Tint>(static_cast<int>(k), p, ListGenEig,
                                         visit_line);
      } else if (mode == 2) {
        std::mt19937_64 rng_point(static_cast<uint64_t>(pos) * 7919ULL +
                                  static_cast<uint64_t>(p) * 104729ULL +
                                  i_round);
        std::uniform_int_distribution<size_t> unif(0, n_line - 1);
        for (size_t i_s = 0; i_s < n_sample; i_s++) {
          if (visit_line(unif(rng_point))) {
            break;
          }
        }
      }
#ifdef DEBUG_GENUS_ENUMERATION_GROUP
      os << "GENUS_GROUP: class " << (pos + 1) << " p=" << p
         << " lambda=" << ed.lambda << " dim=" << k << " orbits visited="
         << n_orbit << "\n";
#endif
      if (closed) {
        return;
      }
    }
  };
  if (prime != 0) {
    if (!IsGoodNeighborPrime(prime, m, det, scalar) ||
        (!abelian && prime == 2)) {
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
    std::vector<int> ListCandidate =
        abelian ? OrderedNeighborPrimes(ListPairWork[0].second, m, det, scalar)
                : OrderedNeighborPrimesModule(ListPairWork[0].first,
                                              ListPairWork[0].second, m, det);
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
        // the deferred orbits of the medium spaces at the primes in use
        bool did_orbits = false;
        for (size_t i = 0; i < deferred.size() &&
                           result.accumulated_mass < TotalMass;
             i++) {
          if (!deferred[i]) {
            continue;
          }
          deferred[i] = 0;
          did_orbits = true;
          orbit_phase = true;
          size_t n_class = ListPairWork.size();
          for (size_t pos = 0;
               pos < n_class && result.accumulated_mass < TotalMass; pos++) {
            expand_at(pos, result.ListPrimeUsed[i], i);
          }
          orbit_phase = false;
        }
        if (did_orbits) {
          continue;
        }
        if (result.ListPrimeUsed.size() == max_prime ||
            result.ListPrimeUsed.size() == ListCandidate.size()) {
          break;
        }
        result.ListPrimeUsed.push_back(ListCandidate[i_work]);
        result.ListNewClassPerPrime.push_back(0);
        next_expand.push_back(0);
        continue;
      }
      expand_at(next_expand[i_work], result.ListPrimeUsed[i_work], i_work);
      next_expand[i_work]++;
    }
  }
  // more random points while the mass is not closed
  while (more_rounds && result.accumulated_mass < TotalMass && i_round < 6) {
    i_round++;
    n_sample *= 4;
    more_rounds = false;
    for (size_t i_prime = 0;
         i_prime < result.ListPrimeUsed.size() &&
         result.accumulated_mass < TotalMass;
         i_prime++) {
      for (size_t pos = 0;
           pos < ListPairWork.size() && result.accumulated_mass < TotalMass;
           pos++) {
        expand_at(pos, result.ListPrimeUsed[i_prime], i_prime);
      }
    }
  }
  result.complete = (result.accumulated_mass == TotalMass);
#ifdef TIMINGS_GENUS_ENUMERATION_GROUP
  os << "|GENUS_GROUP: GenusEnumerationGroup|=" << time_total << "\n";
#endif
  return result;
}

// The cyclic case: one isometry per seed.
template <typename T, typename Tint, typename Tgroup>
GenusEnumerationGroupResult<T, Tgroup> GenusEnumerationGroup(
    std::vector<std::pair<MyMatrix<T>, MyMatrix<T>>> const &ListSeed,
    typename GenusEnumerationGroupResult<T, Tgroup>::Tmass const &TotalMass,
    int prime, size_t max_prime, std::ostream &os) {
  std::vector<std::pair<MyMatrix<T>, std::vector<MyMatrix<T>>>> ListSeedGroup;
  for (auto const &eSeed : ListSeed) {
    ListSeedGroup.push_back({eSeed.first, {eSeed.second}});
  }
  return GenusEnumerationAbelianGroup<T, Tint, Tgroup>(
      ListSeedGroup, TotalMass, prime, max_prime, os);
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
    if (result.accumulated_mass > result.target_mass) {
      // the classes found weigh more than the target: the target is wrong,
      // for a non-abelian group typically a local factor at p | |G|
      os << "mass exceeded = true\n";
    }
    return;
  }
  if (OutFormat == "GAP") {
    os << "return rec(ListGram:=[";
    for (size_t i = 0; i < result.ListGram.size(); i++) {
      os << (i > 0 ? ",\n" : "");
      WriteMatrixGAP(os, result.ListGram[i]);
    }
    // one isometry: ListIsom as in the cyclic case; otherwise the lists
    bool single = true;
    for (auto const &ListIsom : result.ListListIsom) {
      single = single && (ListIsom.size() == 1);
    }
    os << "],\n " << (single ? "ListIsom" : "ListListIsom") << ":=[";
    for (size_t i = 0; i < result.ListListIsom.size(); i++) {
      os << (i > 0 ? ",\n" : "");
      if (single) {
        WriteMatrixGAP(os, result.ListListIsom[i][0]);
      } else {
        os << "[";
        for (size_t j = 0; j < result.ListListIsom[i].size(); j++) {
          os << (j > 0 ? ", " : "");
          WriteMatrixGAP(os, result.ListListIsom[i][j]);
        }
        os << "]";
      }
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
