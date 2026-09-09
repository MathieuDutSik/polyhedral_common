// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_PLESKENSOUVIGNIER_H_
#define SRC_GROUP_PLESKENSOUVIGNIER_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <limits>
#include <optional>
#include <unordered_map>
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
  The Plesken-Souvignier backtrack algorithm for the automorphism group of
  a lattice and for isometry testing:

  W. Plesken, B. Souvignier, Computing isometries of lattices, Journal of
  Symbolic Computation (1997) 24, 327-334.

  The structure follows Bernd Souvignier's programs AUTO and ISOM through
  the port in Hecke.jl (src/QuadForm/Morphism.jl), with one deviation
  noted in ps_stab.

  The problem solved here is the one of PolytopeEquiStabInt.h --
  stabilizer and equivalence of a configuration of vectors under
  GL_n(Z) -- by a different method. There the full graph of pairwise
  scalar products is built and handed to a graph automorphism program,
  which is quadratic in the number of vectors before the graph code even
  starts. Here nothing quadratic in the family is ever built: images of
  one basis are searched in the family, and the pairwise products are
  computed only along the search. On families of thousands of vectors
  with moderate symmetry this is orders of magnitude cheaper; the graph
  method remains preferable when a canonical form is wanted, which the
  backtrack does not provide.

  Conventions, matching the rest of the repository:
  * Vectors are rows. An automorphism g satisfies g * M * g^T = M for
    every matrix M of the configuration and acts on a row vector v as
    v -> v * g.
  * The isometry returned by PleskenSouvignierIsometry satisfies
    P * M1 * P^T = M2, as TestIntEquivalence_ListMat_Vdiag does.

  The vector family S is given as SHVhalf, one representative per
  antipodal pair; the pair set is what matters and the sign of each
  representative is normalized internally. Requirements on the input:

  (1) ListMat[0] is symmetric positive definite and every other matrix is
      symmetric. Positive definiteness is what makes the pointwise
      stabilizer of a basis trivial, which the order computation uses.
  (2) The set {±v} is invariant under the group being computed. Any set
      of the form {v : v M v^T <= bound for the matrices M} is.
  (3) S contains the standard basis e_1, ..., e_n. The searched basis is
      the standard basis, so its images under any automorphism must be in
      S; this also makes every matrix built from a complete assignment
      integral by construction. The set {v : v ListMat[0] v^T <= max_i
      ListMat[0](i,i)} of the paper satisfies this, and an LLL-reduced
      Gram matrix keeps it small. Building that set requires a short
      vector enumeration, which lives in src_latt; see
      src_latt/LatticePleskenSouvignier.h for the wrappers that produce
      it, LLL-reduce, and transport the answers back.

  A complete assignment of images v_1, ..., v_n in S with
  (v_i, v_j)_M = M(i, j) for every M defines the matrix g with rows
  g_i = v_i; the products against the full-rank standard basis force
  g M g^T = M globally, and |det g| = 1 follows from ListMat[0] being
  definite, so g is in GL_n(Z) with no further test. The search prunes
  with the fingerprint of the paper (Section 4): the basis is reordered
  so that the number of candidate images per level is minimal, and a
  partial assignment is abandoned unless the number of candidates at the
  next level is EXACTLY the fingerprint value, an automorphism mapping
  candidate sets bijectively. The automorphism group is assembled level
  by level along the stabilizer chain of the basis (Section 8), with the
  orbits of the already-found group pruning the candidates, and Schreier
  elements harvested from the orbit computations (ps_stab) supplementing
  the generators found by backtracking.

  The vector sums of Section 5 are implemented with a depth parameter,
  defaulting to n/10 rounded as in Hecke; without them the fingerprint
  alone loses rank-14 lattices of moderate symmetry to combinatorial
  explosion (the paper reports 50 h without result at depth 0 on a
  dimension-16 lattice whose depth-5 run takes 14 s). At the level where
  the images v_1, ..., v_l are chosen, every vector w of the family
  carries the tuple of its scalar products with the last `depth` chosen
  images, taken across all forms and normalized in sign. A partial
  automorphism maps the tuple classes of the identity prefix bijectively
  onto those of the chosen prefix, so the branch dies unless the tuples
  hit the precomputed identity classes with the same multiplicities; and
  it maps the per-class vector sums onto each other, so the Gram
  matrices of the sums have to agree as well. Not implemented: the
  Bacher polynomials (Section 6), an additional invariant for the
  remaining hard cases, pluggable into ps_cand the same way.
 */

// One vector per antipodal pair, sign-normalized so that the first
// non-zero coordinate is positive. A point of the search is a signed
// 1-based index: +k is row k-1 of V, -k its negative.
template <typename Tint> struct PleskenSouvignierVectorSet {
  int n;
  int m;
  MyMatrix<Tint> V;
  std::unordered_map<MyVector<Tint>, int> Map;
  // The signed index of w, or 0 when w is not in the family. A 0 from a
  // vector that an automorphism produced means the family is not
  // invariant, which the callers treat as an input error.
  int find_point(MyVector<Tint> const &w) const {
    int len = w.size();
    int ipos = 0;
    while (ipos < len && w(ipos) == 0) {
      ipos++;
    }
    if (ipos == len) {
      return 0;
    }
    if (w(ipos) > 0) {
      auto iter = Map.find(w);
      if (iter == Map.end()) {
        return 0;
      }
      return iter->second + 1;
    }
    MyVector<Tint> wNeg = -w;
    auto iter = Map.find(wNeg);
    if (iter == Map.end()) {
      return 0;
    }
    return -(iter->second + 1);
  }
};

template <typename Tint>
PleskenSouvignierVectorSet<Tint>
PleskenSouvignierBuildVectorSet(MyMatrix<Tint> const &SHVhalf) {
  int m = SHVhalf.rows();
  int n = SHVhalf.cols();
  PleskenSouvignierVectorSet<Tint> VS;
  VS.n = n;
  VS.m = m;
  VS.V = MyMatrix<Tint>(m, n);
  for (int j = 0; j < m; j++) {
    int ipos = 0;
    while (ipos < n && SHVhalf(j, ipos) == 0) {
      ipos++;
    }
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
    if (ipos == n) {
      std::cerr << "PS: the vector family contains the zero vector\n";
      throw TerminalException{1};
    }
#endif
    Tint sign = (SHVhalf(j, ipos) > 0) ? Tint(1) : Tint(-1);
    for (int i = 0; i < n; i++) {
      VS.V(j, i) = sign * SHVhalf(j, i);
    }
    VS.Map[GetMatrixRow(VS.V, j)] = j;
  }
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  if (static_cast<int>(VS.Map.size()) != m) {
    std::cerr << "PS: the vector family has duplicated antipodal pairs\n";
    throw TerminalException{1};
  }
#endif
  return VS;
}

/*
  The vector-sum data of one level (Section 5 of the paper), computed
  for the identity prefix: the distinct sign-normalized tuples of scalar
  products with the trailing dep basis vectors, their multiplicities,
  and the Gram matrices of the per-tuple vector sums. The search
  compares the same data computed for its chosen prefix.
 */
template <typename Tint> struct PleskenSouvignierVectorSumLevel {
  // Number of trailing positions used at this level, min(depth, level).
  int dep;
  std::unordered_map<MyVector<Tint>, int> TupleMap;
  std::vector<int> counts;
  // Whether the tuple is zero: the sign of a zero-tuple vector is not
  // determined, so it is excluded from the sums (but counted).
  std::vector<uint8_t> is_zero;
  // The sums are compared through their pairwise products when the
  // number of tuple classes is small and through their norms otherwise;
  // the switch depends only on the count, an invariant, so isometric
  // lattices take the same branch.
  bool full_gram;
  // full_gram: SumGram[iMat] is the n_tuple x n_tuple product matrix of
  // the sums (zero-tuple rows and columns zero). Otherwise SumNorm[iMat]
  // holds the diagonal alone.
  std::vector<MyMatrix<Tint>> SumGram;
  std::vector<std::vector<Tint>> SumNorm;
  int n_tuple() const { return counts.size(); }
};

/*
  The search context of one lattice: the family, the product tables, the
  fingerprint, and the group data filled by the automorphism search.

  W[iMat] = V * ListMat[iMat], so that the scalar product of rows j and k
  is the dot product of row j of V with row k of W[iMat], and the product
  of row j with a standard basis vector e_p is the single entry
  W[iMat](j, p). Lengths(j, iMat) is the norm of row j. Nothing quadratic
  in the family size is ever built.
 */
template <typename Tint> struct PleskenSouvignierContext {
  std::vector<MyMatrix<Tint>> ListMat;
  PleskenSouvignierVectorSet<Tint> VS;
  std::vector<MyMatrix<Tint>> W;
  MyMatrix<Tint> Lengths;
  // The basis order per and the fingerprint diagonal (candidate count
  // per level under the identity prefix); std_basis[i] is the signed
  // index of e_{per[i]} in the family.
  std::vector<int> per;
  std::vector<int> fp_diagonal;
  std::vector<int> std_basis;
  // g[i]: generators fixing e_{per[0]}, ..., e_{per[i-1]}; the first
  // nsg[i] of them are Schreier elements from ps_stab, the others come
  // from the backtrack. orders[i] is the orbit length of e_{per[i]}
  // under <g[i], ..., g[n-1]>, so the group order is the product.
  std::vector<std::vector<MyMatrix<Tint>>> g;
  std::vector<int> nsg;
  std::vector<int> orders;
  // The vector-sum invariant: depth 0 disables it, and vs_data[I] holds
  // the identity-prefix data of the level with I images assigned, for
  // I = 1, ..., n - 1.
  int depth;
  std::vector<PleskenSouvignierVectorSumLevel<Tint>> vs_data;
  int n() const { return VS.n; }
};

// The scalar product of the row of signed index s with row k0 (0-based)
// for the iMat-th form.
template <typename Tint>
Tint ps_scal_signed(PleskenSouvignierContext<Tint> const &ctx, int s, int k0,
                    int iMat) {
  int j0 = (s > 0 ? s : -s) - 1;
  int n = ctx.n();
  Tint sum(0);
  for (int i = 0; i < n; i++) {
    sum += ctx.VS.V(j0, i) * ctx.W[iMat](k0, i);
  }
  if (s < 0) {
    return -sum;
  }
  return sum;
}

// The image of the point of signed index s under the matrix A, as a
// signed index. The family is required to be invariant, so a miss is an
// input error, not a data case.
template <typename Tint>
int ps_operate(PleskenSouvignierContext<Tint> const &ctx, int s,
               MyMatrix<Tint> const &A) {
  int j0 = (s > 0 ? s : -s) - 1;
  int n = ctx.n();
  MyVector<Tint> w(n);
  for (int i = 0; i < n; i++) {
    Tint sum(0);
    for (int k = 0; k < n; k++) {
      sum += ctx.VS.V(j0, k) * A(k, i);
    }
    w(i) = (s > 0) ? sum : -sum;
  }
  int im = ctx.VS.find_point(w);
  if (im == 0) {
    std::cerr << "PS: an automorphism maps a family vector outside the "
              << "family, so the family is not invariant\n";
    throw TerminalException{1};
  }
  return im;
}

/*
  The count of the paper's fingerprint: the number of signed vectors
  whose norms match those of e_J and whose products with the first I
  basis vectors e_{per[0]}, ..., e_{per[I-1]} match those of e_J. Both
  products are single matrix entries here, the basis being standard.
 */
template <typename Tint>
int ps_possible(PleskenSouvignierContext<Tint> const &ctx, int I, int J) {
  int m = ctx.VS.m;
  int nbMat = ctx.ListMat.size();
  int count = 0;
  for (int j = 0; j < m; j++) {
    bool good_length = true;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      if (ctx.Lengths(j, iMat) != ctx.ListMat[iMat](J, J)) {
        good_length = false;
        break;
      }
    }
    if (!good_length) {
      continue;
    }
    bool okp = true;
    bool okm = true;
    for (int iMat = 0; iMat < nbMat && (okp || okm); iMat++) {
      for (int k = 0; k < I; k++) {
        Tint const &sc = ctx.W[iMat](j, ctx.per[k]);
        Tint const &target = ctx.ListMat[iMat](J, ctx.per[k]);
        if (sc != target) {
          okp = false;
        }
        if (sc != -target) {
          okm = false;
        }
        if (!okp && !okm) {
          break;
        }
      }
    }
    if (okp) {
      count++;
    }
    if (okm) {
      count++;
    }
  }
  return count;
}

/*
  The fingerprint of Section 4: the basis is reordered so that the
  candidate count of each level is minimal given the previous levels, and
  the diagonal counts are kept for the exact-count pruning of ps_cand.
 */
template <typename Tint>
void ps_fingerprint(PleskenSouvignierContext<Tint> &ctx) {
  int n = ctx.n();
  ctx.per.resize(n);
  for (int i = 0; i < n; i++) {
    ctx.per[i] = i;
  }
  std::vector<std::vector<int>> fp(n, std::vector<int>(n, 0));
  for (int i = 0; i < n; i++) {
    fp[0][i] = ps_possible(ctx, 0, i);
  }
  for (int i = 0; i < n - 1; i++) {
    int mini = i;
    for (int j = i + 1; j < n; j++) {
      if (fp[i][ctx.per[j]] < fp[i][ctx.per[mini]]) {
        mini = j;
      }
    }
    std::swap(ctx.per[mini], ctx.per[i]);
    for (int j = i + 1; j < n; j++) {
      fp[i + 1][ctx.per[j]] = ps_possible(ctx, i + 1, ctx.per[j]);
    }
  }
  ctx.fp_diagonal.resize(n);
  for (int i = 0; i < n; i++) {
    ctx.fp_diagonal[i] = fp[i][ctx.per[i]];
  }
}

// Sign-normalizes the vector in place so that its first non-zero entry
// is positive; returns the sign applied, 0 for the zero vector.
template <typename Tint> int ps_normalize_sign(MyVector<Tint> &v) {
  int len = v.size();
  int ipos = 0;
  while (ipos < len && v(ipos) == 0) {
    ipos++;
  }
  if (ipos == len) {
    return 0;
  }
  if (v(ipos) > 0) {
    return 1;
  }
  for (int i = ipos; i < len; i++) {
    v(i) = -v(i);
  }
  return -1;
}

// The identity-prefix vector-sum data of every level (Section 5): the
// tuples come straight out of the W tables, the products with a basis
// vector being single entries.
template <typename Tint>
void ps_init_vector_sums(PleskenSouvignierContext<Tint> &ctx, int depth) {
  ctx.depth = depth;
  if (depth == 0) {
    return;
  }
  int n = ctx.n();
  int m = ctx.VS.m;
  int nbMat = ctx.ListMat.size();
  ctx.vs_data.resize(n);
  for (int I = 1; I < n; I++) {
    PleskenSouvignierVectorSumLevel<Tint> &lev = ctx.vs_data[I];
    int dep = std::min(depth, I);
    lev.dep = dep;
    int len = nbMat * dep;
    std::vector<MyVector<Tint>> sums;
    for (int j = 0; j < m; j++) {
      MyVector<Tint> tup(len);
      int pos = 0;
      for (int iMat = 0; iMat < nbMat; iMat++) {
        for (int k = I - dep; k < I; k++) {
          tup(pos) = ctx.W[iMat](j, ctx.per[k]);
          pos++;
        }
      }
      int sign = ps_normalize_sign(tup);
      int idx;
      auto iter = lev.TupleMap.find(tup);
      if (iter == lev.TupleMap.end()) {
        idx = lev.counts.size();
        lev.TupleMap[tup] = idx;
        lev.counts.push_back(0);
        sums.push_back(ZeroVector<Tint>(n));
      } else {
        idx = iter->second;
      }
      lev.counts[idx]++;
      if (sign != 0) {
        for (int i = 0; i < n; i++) {
          sums[idx](i) += Tint(sign) * ctx.VS.V(j, i);
        }
      }
    }
    int n_tuple = lev.counts.size();
    lev.full_gram = false;
    lev.SumNorm.resize(nbMat);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      lev.SumNorm[iMat].resize(n_tuple);
      for (int idx = 0; idx < n_tuple; idx++) {
        lev.SumNorm[iMat][idx] =
            EvaluationQuadForm(ctx.ListMat[iMat], sums[idx]);
      }
    }
  }
}

/*
  The runtime side of the vector-sum invariant, at the level with I
  images assigned: the tuple of every family vector must hit the
  identity classes with the same multiplicities, and the per-class sums
  must have the identity norms. Sound because an automorphism extending
  the prefix maps the identity classes bijectively onto the runtime
  ones, tuple value to tuple value, and the class sums accordingly.
 */
template <typename Tint>
bool ps_vector_sum_check(PleskenSouvignierContext<Tint> const &Ci,
                         PleskenSouvignierContext<Tint> const &Co, int I,
                         std::vector<int> const &x) {
  PleskenSouvignierVectorSumLevel<Tint> const &lev = Ci.vs_data[I];
  int n = Ci.n();
  int m = Co.VS.m;
  int nbMat = Ci.ListMat.size();
  int dep = lev.dep;
  int len = nbMat * dep;
  int n_tuple = lev.counts.size();
  std::vector<int> counts_run(n_tuple, 0);
  std::vector<MyVector<Tint>> sums_run(n_tuple, ZeroVector<Tint>(n));
  MyVector<Tint> tup(len);
  for (int j = 0; j < m; j++) {
    int pos = 0;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      for (int k = I - dep; k < I; k++) {
        int xk = x[k];
        int row = (xk > 0 ? xk : -xk) - 1;
        Tint sc(0);
        for (int i = 0; i < n; i++) {
          sc += Co.VS.V(j, i) * Co.W[iMat](row, i);
        }
        if (xk < 0) {
          sc = -sc;
        }
        tup(pos) = sc;
        pos++;
      }
    }
    int sign = ps_normalize_sign(tup);
    auto iter = lev.TupleMap.find(tup);
    if (iter == lev.TupleMap.end()) {
      return false;
    }
    int idx = iter->second;
    counts_run[idx]++;
    if (counts_run[idx] > lev.counts[idx]) {
      return false;
    }
    if (sign != 0) {
      for (int i = 0; i < n; i++) {
        sums_run[idx](i) += Tint(sign) * Co.VS.V(j, i);
      }
    }
  }
  // Every count is bounded by its identity value and the totals agree,
  // so the multisets agree.
  for (int iMat = 0; iMat < nbMat; iMat++) {
    for (int idx = 0; idx < n_tuple; idx++) {
      Tint nrm = EvaluationQuadForm(Co.ListMat[iMat], sums_run[idx]);
      if (nrm != lev.SumNorm[iMat][idx]) {
        return false;
      }
    }
  }
  return true;
}

/*
  The context. with_fingerprint is set for the lattice whose basis is
  searched (the only one for automorphisms, the first one for an
  isometry); the fingerprint, the location of the standard basis and the
  vector-sum data are only needed there. depth = -1 selects the default
  of n/10 rounded, as Hecke does; depth = 0 disables the vector sums.
 */
template <typename Tint>
PleskenSouvignierContext<Tint>
PleskenSouvignierBuildContext(std::vector<MyMatrix<Tint>> const &ListMat,
                              MyMatrix<Tint> const &SHVhalf,
                              bool with_fingerprint, int depth,
                              [[maybe_unused]] std::ostream &os) {
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  PleskenSouvignierContext<Tint> ctx;
  ctx.ListMat = ListMat;
  ctx.VS = PleskenSouvignierBuildVectorSet(SHVhalf);
  int n = ctx.VS.n;
  int m = ctx.VS.m;
  int nbMat = ListMat.size();
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  for (auto &eMat : ListMat) {
    if (!IsSymmetricMatrix(eMat)) {
      std::cerr << "PS: the matrices of the configuration have to be "
                << "symmetric\n";
      throw TerminalException{1};
    }
  }
#endif
  ctx.W.reserve(nbMat);
  for (int iMat = 0; iMat < nbMat; iMat++) {
    ctx.W.push_back(ctx.VS.V * ListMat[iMat]);
  }
  ctx.Lengths = MyMatrix<Tint>(m, nbMat);
  for (int j = 0; j < m; j++) {
    for (int iMat = 0; iMat < nbMat; iMat++) {
      Tint sum(0);
      for (int i = 0; i < n; i++) {
        sum += ctx.VS.V(j, i) * ctx.W[iMat](j, i);
      }
      ctx.Lengths(j, iMat) = sum;
    }
  }
  ctx.g.assign(n, {});
  ctx.nsg.assign(n, 0);
  ctx.orders.assign(n, 1);
  ctx.depth = 0;
  // -Id preserves every symmetric form and every antipodal family.
  ctx.g[0].push_back(-IdentityMat<Tint>(n));
  if (with_fingerprint) {
    ps_fingerprint(ctx);
    ctx.std_basis.resize(n);
    for (int i = 0; i < n; i++) {
      MyVector<Tint> e = ZeroVector<Tint>(n);
      e(ctx.per[i]) = 1;
      int k = ctx.VS.find_point(e);
      if (k == 0) {
        std::cerr << "PS: the vector family does not contain the standard "
                  << "basis vector e_" << ctx.per[i] << ". The family must "
                  << "contain e_1, ..., e_n; the set of vectors of norm at "
                  << "most max_i M(i,i) does\n";
        throw TerminalException{1};
      }
      ctx.std_basis[i] = k;
    }
    if (depth == -1) {
      // Hecke's default round(n / 10): 1 from dimension 5 up, 2 from 15.
      depth = (n + 5) / 10;
    }
    ps_init_vector_sums(ctx, depth);
  }
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  os << "|PS: BuildContext|=" << time << "\n";
#endif
  return ctx;
}

/*
  The candidates for the image of the level-I basis vector, given the
  images x[0], ..., x[I-1] already chosen: the signed vectors whose norms
  match e_{per[I]} and whose products with the chosen images match the
  Gram entries. Returns false unless the count is EXACTLY the fingerprint
  value: an automorphism maps the candidate set of the identity prefix
  bijectively to this one, so any other count kills the branch (Section 4
  of the paper).

  Two contexts so that the same function serves the isometry test: the
  basis, fingerprint and target products come from Ci, the vectors from
  Co. For automorphisms both are the same context.
 */
template <typename Tint>
bool ps_cand(PleskenSouvignierContext<Tint> const &Ci,
             PleskenSouvignierContext<Tint> const &Co, int I,
             std::vector<int> const &x, std::vector<int> &candidates) {
  int m = Co.VS.m;
  int nbMat = Ci.ListMat.size();
  int target = Ci.fp_diagonal[I];
  candidates.clear();
  // The vector-sum invariant kills most of the prefixes that the
  // count of candidates alone cannot.
  if (I >= 1 && Ci.depth > 0) {
    if (!ps_vector_sum_check(Ci, Co, I, x)) {
      return false;
    }
  }
  for (int j = 0; j < m; j++) {
    bool okp = true;
    bool okm = true;
    for (int iMat = 0; iMat < nbMat; iMat++) {
      if (Co.Lengths(j, iMat) != Ci.ListMat[iMat](Ci.per[I], Ci.per[I])) {
        okp = false;
        okm = false;
        break;
      }
      for (int k = 0; k < I; k++) {
        int xk = x[k];
        int row = (xk > 0 ? xk : -xk) - 1;
        int n = Ci.n();
        Tint sc(0);
        for (int i = 0; i < n; i++) {
          sc += Co.VS.V(j, i) * Co.W[iMat](row, i);
        }
        if (xk < 0) {
          sc = -sc;
        }
        Tint const &tgt = Ci.ListMat[iMat](Ci.per[I], Ci.per[k]);
        if (sc != tgt) {
          okp = false;
        }
        if (sc != -tgt) {
          okm = false;
        }
        if (!okp && !okm) {
          break;
        }
      }
      if (!okp && !okm) {
        break;
      }
    }
    if (okp) {
      if (static_cast<int>(candidates.size()) >= target) {
        return false;
      }
      candidates.push_back(j + 1);
    }
    if (okm) {
      if (static_cast<int>(candidates.size()) >= target) {
        return false;
      }
      candidates.push_back(-(j + 1));
    }
  }
  return static_cast<int>(candidates.size()) == target;
}

// The matrix of a complete assignment: row per[i] is the vector of
// signed index x[i]. Integral by construction, and in GL_n(Z) whenever
// the products matched (see the header comment).
template <typename Tint>
MyMatrix<Tint> ps_matgen(PleskenSouvignierContext<Tint> const &Ci,
                         PleskenSouvignierContext<Tint> const &Co,
                         std::vector<int> const &x) {
  int n = Ci.n();
  MyMatrix<Tint> X(n, n);
  for (int i = 0; i < n; i++) {
    int xi = x[i];
    int row = (xi > 0 ? xi : -xi) - 1;
    for (int k = 0; k < n; k++) {
      X(Ci.per[i], k) = (xi > 0) ? Co.VS.V(row, k) : -Co.VS.V(row, k);
    }
  }
  return X;
}

// The orbit of a set of signed points under a list of matrices.
template <typename Tint>
std::vector<int> ps_orbit(PleskenSouvignierContext<Tint> const &ctx,
                          std::vector<int> const &pts,
                          std::vector<MyMatrix<Tint>> const &ListGen) {
  int m = ctx.VS.m;
  std::vector<uint8_t> flag(2 * m + 1, 0);
  std::vector<int> orb;
  for (auto &pt : pts) {
    if (flag[pt + m] == 0) {
      flag[pt + m] = 1;
      orb.push_back(pt);
    }
  }
  size_t cnd = 0;
  while (cnd < orb.size()) {
    for (auto &eGen : ListGen) {
      int im = ps_operate(ctx, orb[cnd], eGen);
      if (flag[im + m] == 0) {
        flag[im + m] = 1;
        orb.push_back(im);
      }
    }
    cnd++;
  }
  return orb;
}

// The orbit length of one point, computed no further than cap: the
// caller only compares it with cap.
template <typename Tint>
int ps_orbitlen(PleskenSouvignierContext<Tint> const &ctx, int pt, int cap,
                std::vector<MyMatrix<Tint>> const &ListGen) {
  int m = ctx.VS.m;
  std::vector<uint8_t> flag(2 * m + 1, 0);
  std::vector<int> orb{pt};
  flag[pt + m] = 1;
  size_t cnd = 0;
  while (cnd < orb.size() && static_cast<int>(orb.size()) < cap) {
    for (auto &eGen : ListGen) {
      int im = ps_operate(ctx, orb[cnd], eGen);
      if (flag[im + m] == 0) {
        flag[im + m] = 1;
        orb.push_back(im);
        if (static_cast<int>(orb.size()) >= cap) {
          break;
        }
      }
    }
    cnd++;
  }
  return orb.size();
}

// Removes the members of orb from the ordered candidate list.
inline void ps_setdiff(std::vector<int> &candidates,
                       std::vector<int> const &orb, int m) {
  std::vector<uint8_t> flag(2 * m + 1, 0);
  for (auto &pt : orb) {
    flag[pt + m] = 1;
  }
  std::vector<int> reduced;
  reduced.reserve(candidates.size());
  for (auto &c : candidates) {
    if (flag[c + m] == 0) {
      reduced.push_back(c);
    }
  }
  candidates = std::move(reduced);
}

/*
  The element X1 * G * X2^{-1}: x1 is the assignment of an element X1
  mapping the basis prefix somewhere, x2 of an element X2 with the same
  image of e_{per[I]}, and G a generator mapping the image of X1 to the
  image of X2, so the product stabilizes every basis vector on which the
  two words agree. X2 is in GL_n(Z), so the inverse is integral.
 */
template <typename Tint>
MyMatrix<Tint> ps_stabil(PleskenSouvignierContext<Tint> const &ctx,
                         std::vector<int> const &x1,
                         std::vector<int> const &x2,
                         MyMatrix<Tint> const &G) {
  int n = ctx.n();
  std::vector<int> x(n);
  for (int i = 0; i < n; i++) {
    x[i] = ps_operate(ctx, x1[i], G);
  }
  MyMatrix<Tint> XG = ps_matgen(ctx, ctx, x);
  MyMatrix<Tint> X2 = ps_matgen(ctx, ctx, x2);
  MyMatrix<Tint> S = XG * Inverse(X2);
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  if (S * X2 != XG) {
    std::cerr << "PS: the stabilizer element does not solve S * X2 = XG\n";
    throw TerminalException{1};
  }
#endif
  return S;
}

// A deterministic xorshift, standing in for the random prodding of the
// original: the repository is deterministic throughout, and canonical
// runs must not depend on a seed.
struct PSDeterministicRng {
  uint64_t state = UINT64_C(0x9E3779B97F4A7C15);
  uint64_t next() {
    state ^= state << 13;
    state ^= state >> 7;
    state ^= state << 17;
    return state;
  }
  int next_index(int siz) {
    return static_cast<int>(next() % static_cast<uint64_t>(siz));
  }
};

/*
  Schreier harvesting along the orbit of e_{per[I]} (Section 8 and the
  stab routine of AUTO): while computing the orbit under the generators
  fixing the earlier basis vectors, every coincidence of two words
  yields an element fixing e_{per[0]}, ..., e_{per[I]}, which is kept
  when it enlarges the orbit of a deeper basis vector. These elements
  are much cheaper than ones found by backtracking, and the heuristic
  budget (Maxfail, Rest) bounds the redundancy; the budget only affects
  how many generators are harvested, never correctness.

  Deviation from the Hecke port: the level j at which a harvested
  element is tested is the first index where the two words differ, per
  Souvignier's original code; the port breaks on the first index where
  they agree, which is always I itself.
 */
template <typename Tint>
void ps_stab(PleskenSouvignierContext<Tint> &ctx, int I,
             [[maybe_unused]] std::ostream &os) {
  int n = ctx.n();
  int m = ctx.VS.m;
  int Rest = 0;
  for (int i = I; i < n; i++) {
    if (ctx.fp_diagonal[i] > 1 && ctx.orders[i] < ctx.fp_diagonal[i]) {
      Rest++;
    }
  }
  int Maxfail = Rest;
  for (int i = 0; i < n; i++) {
    if (ctx.fp_diagonal[i] > 1) {
      Maxfail++;
    }
  }
  std::vector<MyMatrix<Tint>> H;
  for (int i = I; i < n; i++) {
    for (auto &eGen : ctx.g[i]) {
      H.push_back(eGen);
    }
  }
  if (H.empty()) {
    return;
  }
  std::vector<std::vector<int>> w(2 * m + 1);
  std::vector<uint8_t> flag(2 * m + 1, 0);
  std::vector<int> orb;
  orb.push_back(ctx.std_basis[I]);
  flag[orb[0] + m] = 1;
  w[orb[0] + m] = ctx.std_basis;
  size_t cnd = 0;
  int fail = 0;
  PSDeterministicRng rng;
  while (cnd < orb.size() && fail < Maxfail + Rest) {
    size_t nH = H.size();
    for (size_t i = 0; i < nH; i++) {
      if (fail >= Maxfail + Rest) {
        break;
      }
      size_t use_cnd = cnd;
      size_t use_i = i;
      if (fail >= Maxfail) {
        // The regular harvest failed Maxfail times: Rest more elements
        // are taken by prodding arbitrary orbit points, which often
        // helps the deeper levels.
        use_cnd = rng.next_index(orb.size());
        use_i = rng.next_index(nH);
      }
      int im = ps_operate(ctx, orb[use_cnd], H[use_i]);
      if (flag[im + m] == 0) {
        flag[im + m] = 1;
        std::vector<int> wim(n);
        for (int j = 0; j < n; j++) {
          wim[j] = ps_operate(ctx, w[orb[use_cnd] + m][j], H[use_i]);
        }
        w[im + m] = std::move(wim);
        orb.push_back(im);
      } else {
        // A coincidence: the first level where the two words differ is
        // where the Schreier element can act.
        int j = I;
        while (j < n) {
          if (ps_operate(ctx, w[orb[use_cnd] + m][j], H[use_i]) !=
              w[im + m][j]) {
            break;
          }
          j++;
        }
        if (j < n && (ctx.orders[j] < ctx.fp_diagonal[j] || fail >= Maxfail)) {
          std::vector<int> x1(n);
          for (int k = 0; k < n; k++) {
            x1[k] = w[orb[use_cnd] + m][k];
          }
          MyMatrix<Tint> S = ps_stabil(ctx, x1, w[im + m], H[use_i]);
          std::vector<MyMatrix<Tint>> Hj{S};
          for (int k = j; k < n; k++) {
            for (auto &eGen : ctx.g[k]) {
              Hj.push_back(eGen);
            }
          }
          int tmplen =
              ps_orbitlen(ctx, ctx.std_basis[j], ctx.fp_diagonal[j], Hj);
          if (tmplen > ctx.orders[j] || fail >= Maxfail) {
            ctx.orders[j] = tmplen;
            ctx.nsg[j]++;
            ctx.g[j].insert(ctx.g[j].begin() + (ctx.nsg[j] - 1), S);
            H.push_back(S);
            if (fail < Maxfail) {
              fail = 0;
            } else {
              fail++;
            }
          } else {
            fail++;
          }
        } else {
          fail++;
        }
      }
    }
    if (fail < Maxfail) {
      cnd++;
    }
  }
}

// The backtrack extension: whether the partial assignment x[0..step-1]
// with the candidate list for level step extends to a complete
// automorphism; on success x holds it.
template <typename Tint>
bool ps_aut_extend(PleskenSouvignierContext<Tint> const &ctx, int step,
                   std::vector<int> &x, std::vector<int> const &cand_step) {
  int n = ctx.n();
  for (auto &c : cand_step) {
    x[step] = c;
    if (step == n - 1) {
      return true;
    }
    std::vector<int> cand_next;
    if (ps_cand(ctx, ctx, step + 1, x, cand_next)) {
      if (ps_aut_extend(ctx, step + 1, x, cand_next)) {
        return true;
      }
    }
  }
  return false;
}

/*
  The automorphism search along the stabilizer chain of the basis
  (Section 8): at level step the group G_step fixing the earlier basis
  vectors is assembled by trying, for each candidate image of e_{per[step]}
  outside the orbit of the group found so far, to extend it to a full
  automorphism. Success gives a new generator and the orbit swallows a
  chunk of candidates; failure kills the whole orbit of the candidate,
  a group element mapping one dead end to another dead end.
 */
template <typename Tint>
void ps_auto(PleskenSouvignierContext<Tint> &ctx, std::ostream &os) {
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  int n = ctx.n();
  int m = ctx.VS.m;
  for (int step = 0; step < n; step++) {
    std::vector<MyMatrix<Tint>> H;
    for (int i = step; i < n; i++) {
      for (auto &eGen : ctx.g[i]) {
        H.push_back(eGen);
      }
    }
    std::vector<int> x(n, 0);
    for (int i = 0; i < step; i++) {
      x[i] = ctx.std_basis[i];
    }
    std::vector<int> cand_step;
    if (ctx.fp_diagonal[step] > 1) {
      [[maybe_unused]] bool test = ps_cand(ctx, ctx, step, x, cand_step);
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
      if (!test) {
        std::cerr << "PS: the identity prefix has a wrong candidate count, "
                  << "the fingerprint is inconsistent\n";
        throw TerminalException{1};
      }
#endif
    } else {
      cand_step.push_back(ctx.std_basis[step]);
    }
    std::vector<int> orb = ps_orbit(ctx, {ctx.std_basis[step]}, H);
    ctx.orders[step] = orb.size();
    ps_setdiff(cand_step, orb, m);
    std::vector<int> bad;
#ifdef DEBUG_PLESKEN_SOUVIGNIER
    os << "PS: step " << step << ", fingerprint " << ctx.fp_diagonal[step]
       << ", candidates outside the orbit " << cand_step.size() << "\n";
#endif
    while (!cand_step.empty()) {
      int im = cand_step[0];
      x[step] = im;
      for (int i = step + 1; i < n; i++) {
        x[i] = 0;
      }
      bool found = false;
      if (step < n - 1) {
        std::vector<int> cand_next;
        if (ps_cand(ctx, ctx, step + 1, x, cand_next)) {
          found = ps_aut_extend(ctx, step + 1, x, cand_next);
        }
      } else {
        found = true;
      }
      if (!found) {
        std::vector<int> oc = ps_orbit(ctx, {im}, H);
        ps_setdiff(cand_step, oc, m);
        bad.push_back(im);
      } else {
        ctx.g[step].push_back(ps_matgen(ctx, ctx, x));
        H.push_back(ctx.g[step].back());
        orb = ps_orbit(ctx, {ctx.std_basis[step]}, H);
        ctx.orders[step] = orb.size();
        ps_setdiff(cand_step, orb, m);
        if (!bad.empty()) {
          std::vector<int> oc = ps_orbit(ctx, bad, H);
          ps_setdiff(cand_step, oc, m);
        }
#ifdef DEBUG_PLESKEN_SOUVIGNIER
        os << "PS: step " << step << ", new generator, orbit "
           << ctx.orders[step] << ", candidates left " << cand_step.size()
           << "\n";
#endif
      }
    }
    if (step < n - 1 && ctx.orders[step] > 1) {
      ps_stab(ctx, step, os);
    }
  }
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  os << "|PS: auto|=" << time << "\n";
#endif
}

/*
  The result of the automorphism computation. The group order is the
  product of the orbit lengths along the stabilizer chain; it is
  returned as the list of factors so that the caller multiplies in
  whatever integer type its orders live in.
 */
template <typename Tint> struct PleskenSouvignierAutomResult {
  std::vector<MyMatrix<Tint>> ListGen;
  std::vector<int> ListOrbitSize;
};

template <typename Tord>
Tord PleskenSouvignierGroupOrder(std::vector<int> const &ListOrbitSize) {
  Tord order(1);
  for (auto &siz : ListOrbitSize) {
    order *= Tord(siz);
  }
  return order;
}

template <typename Tint>
PleskenSouvignierAutomResult<Tint>
PleskenSouvignierAutomorphism(std::vector<MyMatrix<Tint>> const &ListMat,
                              MyMatrix<Tint> const &SHVhalf, std::ostream &os,
                              int depth = -1) {
  PleskenSouvignierContext<Tint> ctx =
      PleskenSouvignierBuildContext(ListMat, SHVhalf, true, depth, os);
  ps_auto(ctx, os);
  PleskenSouvignierAutomResult<Tint> result;
  for (int i = 0; i < ctx.n(); i++) {
    for (int j = ctx.nsg[i]; j < static_cast<int>(ctx.g[i].size()); j++) {
      result.ListGen.push_back(ctx.g[i][j]);
    }
  }
  result.ListOrbitSize = ctx.orders;
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  for (auto &eGen : result.ListGen) {
    for (auto &eMat : ListMat) {
      if (eGen * eMat * eGen.transpose() != eMat) {
        std::cerr << "PS: a generator does not preserve the configuration\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return result;
}

/*
  A bounded harvest of elements of <ListGen> stabilizing the point pt,
  from the coincidences of the orbit transversal. Used only for pruning
  the isometry search, where any subgroup of the true stabilizer gives
  valid orbit deletions, so the bounds trade pruning power for time and
  nothing else.
 */
template <typename Tint>
std::vector<MyMatrix<Tint>>
ps_isostab(PleskenSouvignierContext<Tint> const &ctx, int pt,
           std::vector<MyMatrix<Tint>> const &ListGen) {
  int m = ctx.VS.m;
  int n = ctx.n();
  const size_t max_gens = 16;
  const int max_fail = 32;
  std::vector<MyMatrix<Tint>> H;
  if (ListGen.empty()) {
    return H;
  }
  std::vector<int> orb{pt};
  std::vector<uint8_t> flag(2 * m + 1, 0);
  std::vector<MyMatrix<Tint>> w(1);
  std::vector<int> wpos(2 * m + 1, -1);
  flag[pt + m] = 1;
  wpos[pt + m] = 0;
  w[0] = IdentityMat<Tint>(n);
  size_t cnd = 0;
  int fail = 0;
  while (cnd < orb.size() && fail < max_fail && H.size() < max_gens) {
    for (auto &eGen : ListGen) {
      int im = ps_operate(ctx, orb[cnd], eGen);
      if (flag[im + m] == 0) {
        flag[im + m] = 1;
        wpos[im + m] = w.size();
        w.push_back(w[wpos[orb[cnd] + m]] * eGen);
        orb.push_back(im);
      } else {
        MyMatrix<Tint> B = w[wpos[orb[cnd] + m]] * eGen;
        MyMatrix<Tint> const &Bknown = w[wpos[im + m]];
        if (B != Bknown) {
          MyMatrix<Tint> S = B * Inverse(Bknown);
          bool is_known = false;
          for (auto &eH : H) {
            if (eH == S) {
              is_known = true;
              break;
            }
          }
          if (is_known) {
            fail++;
          } else {
            H.push_back(S);
            if (H.size() >= max_gens) {
              break;
            }
          }
        } else {
          fail++;
        }
      }
    }
    cnd++;
  }
  return H;
}

/*
  The isometry search: candidates for the images of the basis of the
  first lattice among the vectors of the second, with the orbit of a
  failed candidate under the stabilizer of the earlier choices deleted
  wholesale. ListGenAut2 is any set of automorphisms of the second
  configuration; the group they generate powers the orbit deletions, so
  a caller that knows Aut of the second lattice should pass it, and an
  empty list still leaves -Id.
 */
template <typename Tint>
bool ps_iso_extend(PleskenSouvignierContext<Tint> const &Ci,
                   PleskenSouvignierContext<Tint> const &Co, int step,
                   std::vector<int> &x, std::vector<int> &cand_step,
                   std::vector<MyMatrix<Tint>> const &H) {
  int n = Ci.n();
  int m = Co.VS.m;
  while (!cand_step.empty()) {
    int im = cand_step[0];
    x[step] = im;
    if (step == n - 1) {
      return true;
    }
    std::vector<int> cand_next;
    if (ps_cand(Ci, Co, step + 1, x, cand_next)) {
      std::vector<MyMatrix<Tint>> Hnext = ps_isostab(Co, im, H);
      if (ps_iso_extend(Ci, Co, step + 1, x, cand_next, Hnext)) {
        return true;
      }
    }
    std::vector<int> oc = ps_orbit(Co, {im}, H);
    if (oc.size() == 1) {
      // ps_setdiff would only remove im itself; keep the common case
      // cheap.
      cand_step.erase(cand_step.begin());
    } else {
      ps_setdiff(cand_step, oc, m);
    }
  }
  return false;
}

// The isometry, as P with P * ListMat1[i] * P^T = ListMat2[i] for all i,
// or nothing. The two families have to be built by the same invariant
// rule: vectors of norm at most B for the SAME bound B on both sides
// (see the src_latt wrappers).
template <typename Tint>
std::optional<MyMatrix<Tint>>
PleskenSouvignierIsometry(std::vector<MyMatrix<Tint>> const &ListMat1,
                          MyMatrix<Tint> const &SHVhalf1,
                          std::vector<MyMatrix<Tint>> const &ListMat2,
                          MyMatrix<Tint> const &SHVhalf2,
                          std::vector<MyMatrix<Tint>> const &ListGenAut2,
                          std::ostream &os, int depth = -1) {
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  MicrosecondTime time;
#endif
  if (ListMat1.size() != ListMat2.size()) {
    return {};
  }
  if (SHVhalf1.rows() != SHVhalf2.rows()) {
    return {};
  }
  PleskenSouvignierContext<Tint> Ci =
      PleskenSouvignierBuildContext(ListMat1, SHVhalf1, true, depth, os);
  PleskenSouvignierContext<Tint> Co =
      PleskenSouvignierBuildContext(ListMat2, SHVhalf2, false, 0, os);
  int n = Ci.n();
  std::vector<int> x(n, 0);
  std::vector<int> cand0;
  if (!ps_cand(Ci, Co, 0, x, cand0)) {
    return {};
  }
  std::vector<MyMatrix<Tint>> H = ListGenAut2;
  H.push_back(-IdentityMat<Tint>(n));
  bool found = ps_iso_extend(Ci, Co, 0, x, cand0, H);
#ifdef TIMINGS_PLESKEN_SOUVIGNIER
  os << "|PS: isometry found=" << found << "|=" << time << "\n";
#endif
  if (!found) {
    return {};
  }
  // The assignment gives Q with Q * M2 * Q^T = M1; the repository
  // convention is the inverse direction.
  MyMatrix<Tint> Q = ps_matgen(Ci, Co, x);
  MyMatrix<Tint> P = Inverse(Q);
#ifdef SANITY_CHECK_PLESKEN_SOUVIGNIER
  for (size_t iMat = 0; iMat < ListMat1.size(); iMat++) {
    if (P * ListMat1[iMat] * P.transpose() != ListMat2[iMat]) {
      std::cerr << "PS: the isometry does not map ListMat1 to ListMat2\n";
      throw TerminalException{1};
    }
  }
#endif
  return P;
}

// clang-format off
#endif  // SRC_GROUP_PLESKENSOUVIGNIER_H_
// clang-format on
