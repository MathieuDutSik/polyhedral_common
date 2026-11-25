// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_MASSGROUPCOMPONENTS_H_
#define SRC_GENERA_MASSGROUPCOMPONENTS_H_

// clang-format off
#include "MassMultiplicityLocalFactor.h"
#include "MassCyclotomicHermitian.h"
#include "MassAbelianHermitian.h"
#include "MassAbelianQuadratic.h"
#include "MassQuaternionic.h"
#include "MassAlbertTypeTwo.h"
#include "MassUnitaryInnerForm.h"
// clang-format on

#ifdef SANITY_CHECK
# define SANITY_CHECK_ISOMETRY_MASS
#endif

/*
  The components of the mass of a pair (L, G) for a non-abelian finite
  group G (MassIsometry.h): the isotypic components V_rho, End_G(V_rho)
  and its center, and the mass of each according to End_G (Q, a CM or a
  totally real center, a quaternion algebra, Albert type II, a unitary
  inner form), with its local factors at the primes dividing |G|
  (FiniteGroupComponents).
 */

/*
  The rational isotypic components of a finite group G (elements elts, on
  the rows). The rational class sums R_K = sum of the classes of the g^k,
  k prime to the order of g, act on the isotypic part of a rational
  irreducible rho by an integer, and these integers tell the rho apart, so
  the components are the common eigenspaces of the R_K; their eigenvalues
  are found numerically and checked exactly.

  On a component where the generators commute, G acts through a cyclic
  quotient and the component is that of the cyclic case for the element
  acting with the largest order. Otherwise, End_G(V_rho) = Q is recognized
  as a center of dimension 1, End_G of dimension m^2, and a common
  eigenspace F of some elements h for eigenvalues +-1 of dimension m: F is
  then the image l (x) U of an idempotent of Q[H] of rank one on W, and X
  = L cap F, with the restriction of the form, is the multiplicity lattice
  (Section on non-abelian groups of the paper). Other components get no
  mass.
 */
/*
  An End_G-submodule of dimension dim of a block (basis E, rows, in the
  coordinates of L) as End_G x: End_G y for y a pure tensor of W (x) U is
  l (x) U, l the line of y, of dimension (multiplicity) x [K : Q]. The y
  tried are the basis vectors of the block, then their pairwise sums; with
  multiplicity one every nonzero y is a pure tensor, so the first one
  succeeds. EndBlk is End_G in the block coordinates, acting on the rows as
  the elements do. Returned: the rows of an integral basis of L cap F, or
  nothing. F is a valid choice for the multiplicity lattice: its orthogonal
  projection commutes with End_G (closed under the adjoint), so it lies in
  the span of G, and at p not dividing |G| a line through a primitive
  vector gives L_p = Z_p[G] X_p over the maximal order Z_p[G] e_rho.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
EndomorphismSpanSubmodule(MyMatrix<T> const &E,
                          std::vector<MyMatrix<T>> const &EndBlk, int dim) {
  int k = E.rows();
  int n = E.cols();
  auto try_vector = [&](MyVector<T> const &y) -> std::optional<MyMatrix<T>> {
    MyMatrix<T> span(static_cast<int>(EndBlk.size()), k);
    for (size_t i = 0; i < EndBlk.size(); i++) {
      MyVector<T> z = EndBlk[i].transpose() * y;
      for (int j = 0; j < k; j++) {
        span(static_cast<int>(i), j) = z(j);
      }
    }
    if (RankMat(span) != dim) {
      return {};
    }
    // L cap F, F the rows of span E: the integral x killed by the
    // annihilator of F
    MyMatrix<T> F = span * E;
    MyMatrix<T> Ann = NullspaceMat(MyMatrix<T>(F.transpose()));
    MyMatrix<Tint> Ann_i(Ann.rows(), n);
    // the annihilator rows with their denominators cleared
    for (int i = 0; i < Ann.rows(); i++) {
      Tint den(1);
      for (int j = 0; j < n; j++) {
        den = LCMpair(den, GetDenominator_z(Ann(i, j)));
      }
      for (int j = 0; j < n; j++) {
        Ann_i(i, j) = UniversalScalarConversion<Tint, T>(
            Ann(i, j) * UniversalScalarConversion<T, Tint>(den));
      }
    }
    MyMatrix<Tint> K_i = NullspaceIntTrMat(Ann_i);
    if (K_i.rows() != dim) {
      return {};
    }
    return UniversalMatrixConversion<T, Tint>(K_i);
  };
  for (int i = 0; i < k; i++) {
    MyVector<T> y = ZeroVector<T>(k);
    y(i) = T(1);
    if (std::optional<MyMatrix<T>> opt = try_vector(y)) {
      return opt;
    }
  }
  for (int i = 0; i < k; i++) {
    for (int j = i + 1; j < k; j++) {
      MyVector<T> y = ZeroVector<T>(k);
      y(i) = T(1);
      y(j) = T(1);
      if (std::optional<MyMatrix<T>> opt = try_vector(y)) {
        return opt;
      }
    }
  }
  return {};
}

/*
  The search of a submodule of a block by kernels. The block is {x : x C = 0
  for C in conditions} (rows, as the elements); the candidates are the
  subspaces S_c = {x in the block : x M_c = 0}. Returned: the first S_c of
  dimension dim, otherwise the first S_a cap S_b of dimension dim with
  owner[a] != owner[b] (the pairs in lexicographic order), as the rows of an
  integral basis; or nothing.

  The dimensions are read modulo a prime p first, by linear algebra over
  F_p in the coordinates of the block: the basis B of the block once, each
  S_c as the left kernel of B M_c, and dim (S_a cap S_b) = dim S_a + dim S_b
  - rank [S_a; S_b]. A candidate of the right dimension modulo p is then
  computed exactly (an integral kernel) and kept only if its dimension over
  Q is dim. The matrices are integral, and the dimension over Q is at most
  the one modulo p, so a candidate can be missed only when p divides a
  minor; nothing wrong is ever returned, a miss only leaves the component
  without a submodule.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
FirstKernelOfDimension(std::vector<MyMatrix<T>> const &conditions,
                       std::vector<MyMatrix<T>> const &cand,
                       std::vector<size_t> const &owner, int n, int dim) {
  int64_t const p = 2147483647;
  auto to_mod = [&](MyMatrix<T> const &M) -> ModRowsI64 {
    ModRowsI64 R(M.rows(), std::vector<int64_t>(M.cols()));
    for (int i = 0; i < M.rows(); i++) {
      for (int j = 0; j < M.cols(); j++) {
#ifdef SANITY_CHECK_ISOMETRY_MASS
        if (!IsInteger(M(i, j))) {
          std::cerr << "ISOM_MASS: a condition is not integral\n";
          throw TerminalException{1};
        }
#endif
        R[i][j] = ResModOfI64(GetNumerator_z(M(i, j)), p);
      }
    }
    return R;
  };
  auto hconcat = [&](std::vector<MyMatrix<T>> const &l) -> MyMatrix<T> {
    MyMatrix<T> Big(n, n * static_cast<int>(l.size()));
    for (size_t a = 0; a < l.size(); a++) {
      for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
          Big(i, static_cast<int>(a) * n + j) = l[a](i, j);
        }
      }
    }
    return Big;
  };
  // the exact integral kernel of the conditions and the candidates cs
  auto exact = [&](std::vector<size_t> const &cs) -> std::optional<MyMatrix<T>> {
    std::vector<MyMatrix<T>> all = conditions;
    for (size_t c : cs) {
      all.push_back(cand[c]);
    }
    MyMatrix<Tint> Big_i =
        UniversalMatrixConversion<Tint, T>(hconcat(all));
    MyMatrix<Tint> K_i = NullspaceIntTrMat(MyMatrix<Tint>(Big_i.transpose()));
    if (K_i.rows() != dim) {
      return {};
    }
    return UniversalMatrixConversion<T, Tint>(K_i);
  };
  // the block modulo p, then each S_c in its coordinates
  ModRowsI64 B;
  if (conditions.empty()) {
    B.assign(n, std::vector<int64_t>(n, 0));
    for (int i = 0; i < n; i++) {
      B[i][i] = 1;
    }
  } else {
    B = LeftNullspaceModPI64(to_mod(hconcat(conditions)), p);
  }
  std::vector<ModRowsI64> S(cand.size());
  for (size_t c = 0; c < cand.size(); c++) {
    S[c] = B.empty() ? ModRowsI64{} : LeftNullspaceModPI64(ModRowsProductI64(B, to_mod(cand[c]), p), p);
  }
  for (size_t a = 0; a < cand.size(); a++) {
    if (static_cast<int>(S[a].size()) == dim) {
      if (std::optional<MyMatrix<T>> opt = exact({a})) {
        return opt;
      }
    }
  }
  for (size_t a = 0; a < cand.size(); a++) {
    if (static_cast<int>(S[a].size()) < dim) {
      continue;
    }
    for (size_t b = a + 1; b < cand.size(); b++) {
      if (owner[a] == owner[b] || static_cast<int>(S[b].size()) < dim) {
        continue;
      }
      ModRowsI64 stacked = S[a];
      stacked.insert(stacked.end(), S[b].begin(), S[b].end());
      int d = static_cast<int>(S[a].size() + S[b].size()) -
              static_cast<int>(RankModPrimeI64(stacked, p));
      if (d == dim) {
        if (std::optional<MyMatrix<T>> opt = exact({a, b})) {
          return opt;
        }
      }
    }
  }
  return {};
}

/*
  An End_G-submodule of the block of dimension dim, as the integral kernel
  of Phi_d(g), for one element or two, d dividing the order of g, together
  with the conditions of the block (row convention, as the elements): the
  rows of a basis of L cap F, or nothing.
 */
template <typename T, typename Tint>
std::optional<MyMatrix<T>>
CyclotomicKernelSubmodule(std::vector<MyMatrix<T>> const &elts,
                          std::vector<MyMatrix<T>> const &conditions, int n,
                          int dim) {
  std::vector<std::pair<size_t, MyMatrix<T>>> cand;
  for (size_t i = 0; i < elts.size(); i++) {
    int o = IsometryOrder(elts[i]);
    for (int dd = 1; dd <= o; dd++) {
      if (o % dd == 0) {
        cand.push_back({i, EvaluatePolynomialMatrix<T>(CyclotomicPolynomial(dd),
                                                       elts[i])});
      }
    }
  }
  std::vector<MyMatrix<T>> mats;
  std::vector<size_t> owner;
  for (auto const &c : cand) {
    owner.push_back(c.first);
    mats.push_back(c.second);
  }
  return FirstKernelOfDimension<T, Tint>(conditions, mats, owner, n, dim);
}

template <typename T, typename Tint>
std::vector<IsometryComponent<T>>
FiniteGroupComponents(MyMatrix<T> const &A, std::vector<MyMatrix<T>> const &ListP,
                      std::vector<MyMatrix<T>> const &elts,
                      [[maybe_unused]] std::ostream &os) {
  int n = A.rows();
  /*
    The group computations (classes, inverses, orders, powers) on the
    elements as integer matrices, flattened, with products in 64 bits
    (checked) and a hash index: exact rational products were most of the
    time of the components for groups of a few hundred elements.
   */
  size_t n_elt = elts.size();
  using IMat = std::vector<int64_t>;
  struct IMatHash {
    size_t operator()(IMat const &v) const {
      size_t h = 1469598103934665603ULL;
      for (int64_t x : v) {
        h ^= static_cast<size_t>(x) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
      }
      return h;
    }
  };
  std::vector<IMat> ielts(n_elt, IMat(n * n));
  for (size_t i = 0; i < n_elt; i++) {
    for (int u = 0; u < n; u++) {
      for (int v = 0; v < n; v++) {
        T x = elts[i](u, v);
#ifdef SANITY_CHECK_ISOMETRY_MASS
        if (!IsInteger(x)) {
          std::cerr << "ISOM_MASS: a group element is not integral\n";
          throw TerminalException{1};
        }
#endif
        ielts[i][u * n + v] = UniversalScalarConversion<int64_t, T>(x);
      }
    }
  }
  auto imul = [&](IMat const &X, IMat const &Y) -> IMat {
    IMat Z(n * n, 0);
    for (int u = 0; u < n; u++) {
      for (int w = 0; w < n; w++) {
        int64_t xv = X[u * n + w];
        if (xv != 0) {
          for (int v = 0; v < n; v++) {
            __int128 z = static_cast<__int128>(Z[u * n + v]) +
                         static_cast<__int128>(xv) * Y[w * n + v];
            if (z > INT64_MAX || z < INT64_MIN) {
              std::cerr << "ISOM_MASS: a product of group elements overflows\n";
              throw TerminalException{1};
            }
            Z[u * n + v] = static_cast<int64_t>(z);
          }
        }
      }
    }
    return Z;
  };
  std::unordered_map<IMat, size_t, IMatHash> iindex;
  for (size_t i = 0; i < n_elt; i++) {
    iindex[ielts[i]] = i;
  }
  auto ifind = [&](IMat const &X) -> size_t {
    auto iter = iindex.find(X);
#ifdef SANITY_CHECK_ISOMETRY_MASS
    if (iter == iindex.end()) {
      std::cerr << "ISOM_MASS: a product is not in the group\n";
      throw TerminalException{1};
    }
#endif
    return iter->second;
  };
  size_t i_id = 0;
  {
    IMat Id(n * n, 0);
    for (int u = 0; u < n; u++) {
      Id[u * n + u] = 1;
    }
    i_id = ifind(Id);
  }
  // the orders and the inverses, by the powers
  std::vector<int> iorder(n_elt);
  std::vector<size_t> iinv(n_elt);
  for (size_t i = 0; i < n_elt; i++) {
    IMat pw = ielts[i];
    size_t prev = i;
    int o = 1;
    while (true) {
      size_t ip = ifind(pw);
      if (ip == i_id) {
        break;
      }
      prev = ip;
      pw = imul(pw, ielts[i]);
      o++;
    }
    iorder[i] = o;
    iinv[i] = (o == 1) ? i : prev;
  }
  auto conj = [&](size_t j, size_t u) -> size_t {
    return ifind(imul(imul(ielts[iinv[j]], ielts[u]), ielts[j]));
  };
  // the index of x^a
  auto ipow = [&](size_t i, int64_t a) -> size_t {
    IMat pw = ielts[i_id];
    int64_t e2 = a % iorder[i];
    for (int64_t t = 0; t < e2; t++) {
      pw = imul(pw, ielts[i]);
    }
    return ifind(pw);
  };
  // the rational classes
  std::vector<int> cls(n_elt, -1);
  int n_cls = 0;
  for (size_t i = 0; i < n_elt; i++) {
    if (cls[i] != -1) {
      continue;
    }
    std::vector<size_t> queue{i};
    cls[i] = n_cls;
    while (!queue.empty()) {
      size_t u = queue.back();
      queue.pop_back();
      std::vector<size_t> nb;
      for (size_t j = 0; j < n_elt; j++) {
        nb.push_back(conj(j, u));
      }
      // the powers prime to the order
      int o = iorder[u];
      IMat pw = ielts[u];
      for (int k = 2; k < o; k++) {
        pw = imul(pw, ielts[u]);
        if (std::gcd(k, o) == 1) {
          nb.push_back(ifind(pw));
        }
      }
      for (size_t v : nb) {
        if (cls[v] == -1) {
          cls[v] = n_cls;
          queue.push_back(v);
        }
      }
    }
    n_cls++;
  }
  std::vector<MyMatrix<T>> ListR(n_cls, ZeroMatrix<T>(n, n));
  for (size_t i = 0; i < n_elt; i++) {
    ListR[cls[i]] += elts[i];
  }
  // the common eigenspaces, refined class sum by class sum
  struct Block {
    MyMatrix<T> E;
    std::vector<MyMatrix<T>> conditions;
  };
  std::vector<Block> blocks{{IdentityMat<T>(n), {}}};
  for (auto const &Rk : ListR) {
    std::vector<Block> next;
    for (auto const &blk : blocks) {
      int k = blk.E.rows();
      // the action of R_K on the block, numerically
      MyMatrix<T> act =
          blk.E * Rk * blk.E.transpose() *
          Inverse(MyMatrix<T>(blk.E * blk.E.transpose()));
      Eigen::MatrixXd act_d(k, k);
      for (int i = 0; i < k; i++) {
        for (int j = 0; j < k; j++) {
          act_d(i, j) = UniversalScalarConversion<double, T>(act(i, j));
        }
      }
      Eigen::EigenSolver<Eigen::MatrixXd> es(act_d, false);
      std::set<int64_t> ListLambda;
      for (int i = 0; i < k; i++) {
        ListLambda.insert(static_cast<int64_t>(std::llround(es.eigenvalues()(i).real())));
      }
      int tot = 0;
      for (int64_t lambda : ListLambda) {
        MyMatrix<T> M = Rk - UniversalScalarConversion<T, int64_t>(lambda) * IdentityMat<T>(n);
        MyMatrix<T> EM = blk.E * M;
        MyMatrix<Tint> EM_i = UniversalMatrixConversion<Tint, T>(EM);
        MyMatrix<Tint> K_i = NullspaceIntTrMat(MyMatrix<Tint>(EM_i.transpose()));
        if (K_i.rows() == 0) {
          continue;
        }
        MyMatrix<T> K = UniversalMatrixConversion<T, Tint>(K_i) * blk.E;
        std::vector<MyMatrix<T>> cond = blk.conditions;
        cond.push_back(M);
        next.push_back({K, cond});
        tot += K.rows();
      }
#ifdef SANITY_CHECK_ISOMETRY_MASS
      if (tot != k) {
        std::cerr << "ISOM_MASS: the class sums are not diagonalizable with "
                     "integer eigenvalues on a block\n";
        throw TerminalException{1};
      }
#endif
    }
    blocks = std::move(next);
  }
  std::vector<IsometryComponent<T>> res;
  for (auto const &blk : blocks) {
    int k = blk.E.rows();
    MyMatrix<T> Gram = blk.E * A * blk.E.transpose();
    MyMatrix<T> GramInv = Inverse(Gram);
    auto action = [&](MyMatrix<T> const &P) -> MyMatrix<T> {
      return blk.E * P * A * blk.E.transpose() * GramInv;
    };
    std::vector<MyMatrix<T>> ListAct;
    for (auto const &P : ListP) {
      ListAct.push_back(action(P));
    }
    bool commuting = true;
    for (size_t i = 0; i < ListAct.size(); i++) {
      for (size_t j = i + 1; j < ListAct.size(); j++) {
        if (ListAct[i] * ListAct[j] != ListAct[j] * ListAct[i]) {
          commuting = false;
        }
      }
    }
    if (commuting) {
      // a cyclic quotient: the element of largest order on the block
      size_t i_best = 0;
      int o_best = 0;
      for (size_t i = 0; i < n_elt; i++) {
        int o = IsometryOrder(action(elts[i]));
        if (o > o_best) {
          o_best = o;
          i_best = i;
        }
      }
      res.push_back(
          BuildIsometryComponent<T, Tint>(A, blk.E, elts[i_best], o_best));
      continue;
    }
    IsometryComponent<T> comp =
        BuildIsometryComponent<T, Tint>(A, blk.E, ListP[0], 1);
    comp.d = 0;
    // set by the method that gives the mass; a component left without one
    // stops the program at the end of the loop
    comp.kind = ComponentKind::Unsupported;
    for (size_t i = 1; i < ListP.size(); i++) {
      comp.extra_sigma.push_back(action(ListP[i]));
    }
    // End_G of the block and its center
    auto commutant = [&](std::vector<MyMatrix<T>> const &ListX)
        -> std::vector<MyMatrix<T>> {
      // the Y with X Y = Y X for all X, by the nullspace
      int kk = k * k;
      MyMatrix<T> Eq(kk, kk * static_cast<int>(ListX.size()));
      for (size_t ix = 0; ix < ListX.size(); ix++) {
        MyMatrix<T> const &X = ListX[ix];
        for (int u = 0; u < k; u++) {
          for (int v = 0; v < k; v++) {
            int col = static_cast<int>(ix) * kk + u * k + v;
            for (int w = 0; w < kk; w++) {
              Eq(w, col) = 0;
            }
            // (X Y - Y X)(u, v) = sum_w X(u,w) Y(w,v) - Y(u,w) X(w,v)
            for (int w = 0; w < k; w++) {
              Eq(w * k + v, col) += X(u, w);
              Eq(u * k + w, col) -= X(w, v);
            }
          }
        }
      }
      MyMatrix<T> ker = NullspaceMat(Eq);
      std::vector<MyMatrix<T>> res_c;
      for (int r = 0; r < ker.rows(); r++) {
        MyMatrix<T> Y(k, k);
        for (int u = 0; u < k; u++) {
          for (int v = 0; v < k; v++) {
            Y(u, v) = ker(r, u * k + v);
          }
        }
        res_c.push_back(Y);
      }
      return res_c;
    };
    // the class sums, the exponent N of G and the a with sigma_a fixing
    // the center K of End_G (the subfield of Q(zeta_N) it is)
    auto center_field = [&]()
        -> std::tuple<std::vector<MyMatrix<T>>, int64_t, std::vector<int64_t>> {
      // the conjugacy classes and their sums on the block
      std::vector<int> ccls(n_elt, -1);
      int n_ccls = 0;
      for (size_t i = 0; i < n_elt; i++) {
        if (ccls[i] != -1) {
          continue;
        }
        for (size_t j = 0; j < n_elt; j++) {
          ccls[conj(j, i)] = n_ccls;
        }
        n_ccls++;
      }
      std::vector<MyMatrix<T>> ClassSum(n_ccls, ZeroMatrix<T>(n, n));
      for (size_t i = 0; i < n_elt; i++) {
        ClassSum[ccls[i]] += elts[i];
      }
      int64_t Nexp = 1;
      for (size_t i = 0; i < n_elt; i++) {
        Nexp = std::lcm(Nexp, static_cast<int64_t>(iorder[i]));
      }
      // the class sums on the block, once
      std::vector<MyMatrix<T>> ActClassSum;
      for (auto const &Cs : ClassSum) {
        ActClassSum.push_back(action(Cs));
      }
      std::vector<int64_t> Hgens;
      for (int64_t a = 1; a < Nexp; a++) {
        if (std::gcd(a, Nexp) != 1) {
          continue;
        }
        bool fixes = true;
        for (size_t i = 0; i < n_elt && fixes; i++) {
          int ci = ccls[i];
          int ca = ccls[ipow(i, a)];
          if (ci != ca && ActClassSum[ci] != ActClassSum[ca]) {
            fixes = false;
          }
        }
        if (fixes) {
          Hgens.push_back(a);
        }
      }
      return {ClassSum, Nexp, Hgens};
    };
    std::vector<MyMatrix<T>> End = commutant(ListAct);
    /*
      The center of End_G. End_G is the commutant of the image A of Q[G] on
      the block, A is semisimple, so the center of End_G is that of A'' = A:
      the span of the conjugacy class sums on the block. Computed as the
      commutant of A and End_G together, it was a system of k^2 unknowns
      and (1 + dim End_G) k^2 equations, beyond 13 GB for 2I diagonal on
      I^4 (k = 32, dim End_G = 128).
     */
    std::vector<MyMatrix<T>> CenterBasis = [&]() -> std::vector<MyMatrix<T>> {
      std::vector<int> ccls(n_elt, -1);
      int n_ccls = 0;
      for (size_t i = 0; i < n_elt; i++) {
        if (ccls[i] != -1) {
          continue;
        }
        for (size_t j = 0; j < n_elt; j++) {
          ccls[conj(j, i)] = n_ccls;
        }
        n_ccls++;
      }
      std::vector<MyMatrix<T>> Sums(n_ccls, ZeroMatrix<T>(n, n));
      for (size_t i = 0; i < n_elt; i++) {
        Sums[ccls[i]] += elts[i];
      }
      std::vector<MyMatrix<T>> Acts;
      MyMatrix<T> Fl(n_ccls, k * k);
      for (int ic = 0; ic < n_ccls; ic++) {
        Acts.push_back(action(Sums[ic]));
        for (int u = 0; u < k; u++) {
          for (int v = 0; v < k; v++) {
            Fl(ic, u * k + v) = Acts[ic](u, v);
          }
        }
      }
      std::vector<MyMatrix<T>> basis;
      for (int i_row : TMat_SelectRowCol(Fl).ListRowSelect) {
        basis.push_back(Acts[i_row]);
      }
      return basis;
    }();
#ifdef SANITY_CHECK_ISOMETRY_MASS
    {
      std::vector<MyMatrix<T>> ListBoth = ListAct;
      ListBoth.insert(ListBoth.end(), End.begin(), End.end());
      if (commutant(ListBoth).size() != CenterBasis.size()) {
        std::cerr << "ISOM_MASS: the class sums do not span the center\n";
        throw TerminalException{1};
      }
    }
#endif
    int c = static_cast<int>(CenterBasis.size());
    int e = static_cast<int>(End.size());
    int m = static_cast<int>(std::llround(std::sqrt(static_cast<double>(e))));
    bool central_q = (c == 1);
    {
      // End_G = M_r(D), D a quaternion algebra (Albert type III), in the
      // columns convention
      std::vector<MyMatrix<T>> EndC, CenterC, GensC, EltsC;
      for (auto const &Y : End) {
        EndC.push_back(Y.transpose());
      }
      for (auto const &Z : CenterBasis) {
        CenterC.push_back(Z.transpose());
      }
      int rq = QuaternionicRank(Gram, EndC, CenterC);
      // Albert type II over D: the same reduction to X, the mass of
      // AlbertTypeTwoRankOneMass (r = 1) or AlbertTypeTwoOrthogonalMass.
      // For r >= 2 the split M_{2r}(Q) of type I has the same dimensions:
      // D is only seen on X, and when it is split (or G does not span an
      // M_s(D)) the block goes on to the type I path.
      bool type2 = false;
      bool type2_split = false;
      std::optional<AbelianField> FKq_holder;
      T dKq_holder(1);
      if (rq == 0) {
        int r2 = AlbertTypeTwoRankCandidate(Gram, EndC, CenterC);
        if (r2 == 1) {
          // the center K, its discriminant, and D a division algebra
          if (c > 1) {
            auto [ClassSum2, Nexp2, Hgens2] = center_field();
            (void)ClassSum2;
            FKq_holder = MakeAbelianRealField(Nexp2, Hgens2);
            DirichletGroup Gd = ComputeDirichletGroup(FKq_holder->n);
            for (auto const &t : FieldCharacters(Gd, FKq_holder->H, false)) {
              dKq_holder *= UniversalScalarConversion<T, int64_t>(PrimitiveCharacter(Gd, t).first);
            }
          }
          if (EndIsQuaternionDivision<T, Tint>(EndC, k, c, dKq_holder,
                                               QuaternionicMassBounds(), os)) {
            rq = 1;
            type2 = true;
          }
        } else if (r2 >= 2 && c == 1) {
          rq = r2;
          type2 = true;
        }
      }
      auto quat_mass = [&](MyMatrix<T> const &Gx,
                           std::vector<MyMatrix<T>> const &Ds,
                           QuaternionicMassBounds const &qb) -> std::optional<T> {
        if (!type2) {
          std::optional<T> qm = QuaternionicComponentMass<T, Tint>(
              Gx, Ds, FKq_holder, dKq_holder, rq, qb, os);
#ifdef SANITY_CHECK_ALBERT_TYPE_TWO_MASS_EXTENSIVE
          // the measures of the Tamagawa route on type III, rank one over Q:
          // a second mass computation, expensive
          std::cerr << "ISOM_MASS: SANITY_CHECK_ALBERT_TYPE_TWO_MASS_EXTENSIVE\n";
          if (qm && rq == 1 && c == 1) {
            std::optional<T> tm = TypeThreeRankOneTamagawaMass<T, Tint>(
                Gx, Ds, AlbertTypeTwoMassBounds{qb}, os);
            os << "ISOM_MASS: type III rank one, mass " << *qm
               << ", by the Tamagawa number " << (tm ? *tm : T(-1)) << "\n";
            if (tm && *tm != *qm) {
              std::cerr << "ISOM_MASS: the Tamagawa route disagrees on type III\n";
              throw TerminalException{1};
            }
          }
          std::cerr << "ISOM_MASS: SANITY_CHECK_ALBERT_TYPE_TWO_MASS_EXTENSIVE, done\n";
#endif
          return qm;
        }
        if (rq == 1) {
          AlbertTypeTwoMassBounds b2;
          b2.quat = qb;
          return AlbertTypeTwoRankOneMass<T, Tint>(Gx, Ds, FKq_holder,
                                                   dKq_holder, b2, os);
        }
        if (!QuaternionSpanIsDivision<T, Tint>(Ds, Gx.rows(), rq, qb, os)) {
          type2_split = true;
          return {};
        }
        return AlbertTypeTwoOrthogonalMass<T, Tint>(
            Gx, Ds, rq, AlbertTypeTwoMassBounds{qb, QuadLValueMethod::Mpfr}, os);
      };
      if (rq > 0) {
        for (auto const &M : ListAct) {
          GensC.push_back(M.transpose());
        }
        for (auto const &x : elts) {
          EltsC.push_back(action(x).transpose());
        }
#ifdef SANITY_CHECK_ISOMETRY_MASS
        for (auto const &C : GensC) {
          if (MyMatrix<T>(C.transpose() * Gram * C) != Gram) {
            std::cerr << "ISOM_MASS: a generator is not an isometry of the "
                         "block\n";
            throw TerminalException{1};
          }
        }
#endif
        std::optional<AbelianField> FKq;
        T dKq(1);
        if (c > 1) {
          auto [ClassSumq, Nexpq, Hgensq] = center_field();
          (void)ClassSumq;
          FKq = MakeAbelianRealField(Nexpq, Hgensq);
#ifdef SANITY_CHECK_ISOMETRY_MASS
          if (AbelianFieldDegree(*FKq) != c) {
            std::cerr << "ISOM_MASS: the center has degree " << c
                      << " but the field found has degree "
                      << AbelianFieldDegree(*FKq) << "\n";
            throw TerminalException{1};
          }
#endif
          DirichletGroup Gd = ComputeDirichletGroup(FKq->n);
          for (auto const &t : FieldCharacters(Gd, FKq->H, false)) {
            dKq *= UniversalScalarConversion<T, int64_t>(PrimitiveCharacter(Gd, t).first);
          }
        }
        comp.kind = ComponentKind::Precomputed;
        QuaternionicMassBounds qbounds;
        // G spans M_s(D) on the block, of dimension 4 c s^2
        int dimA = 0;
        {
          MyMatrix<T> Fl(static_cast<int>(EltsC.size()), k * k);
          for (size_t i = 0; i < EltsC.size(); i++) {
            for (int u = 0; u < k; u++) {
              for (int v = 0; v < k; v++) {
                Fl(static_cast<int>(i), u * k + v) = EltsC[i](u, v);
              }
            }
          }
          dimA = RankMat(Fl);
        }
        int sq = static_cast<int>(std::llround(
            std::sqrt(static_cast<double>(dimA) / static_cast<double>(4 * c))));
        if (4 * c * sq * sq != dimA || k != 4 * c * sq * rq) {
#ifdef SANITY_CHECK_ISOMETRY_MASS
          if (!(type2 && rq >= 2)) {
            std::cerr << "ISOM_MASS: G spans an algebra of dimension " << dimA
                      << " on a quaternionic block of rank " << k << "\n";
            throw TerminalException{1};
          }
#endif
          type2_split = true;
        }
        if (!type2) {
          FKq_holder = FKq;
          dKq_holder = dKq;
        }
        if (type2_split) {
          // the split M_{2r}(Q): the type I path below
        } else if (sq == 1) {
          // X = L_rho, on which G spans D
          comp.precomputed_mass = quat_mass(Gram, EltsC, qbounds);
          comp.multiplicity_factor = T(1);
        } else {
          /*
            X = L cap F, F = e V_rho for a primitive self-adjoint idempotent
            e of the span of G: any End_G-submodule of V_rho of dimension 4
            c r (its orthogonal projection commutes with End_G, so it is in
            the span of G, and it is self-adjoint). F is found as the kernel
            of Phi_d(g) for one or two elements, d | ord(g). D acts on F by
            the e g e restricted to F, and the mass of the component is
            that of X times the local factors [U(X_p) : Stab(L_p)].
           */
          int dimF = 4 * c * rq;
          std::vector<MyMatrix<T>> cand;
          for (size_t i = 0; i < n_elt; i++) {
            int o = iorder[i];
            for (int dd = 1; dd <= o; dd++) {
              if (o % dd == 0) {
                cand.push_back(EvaluatePolynomialMatrix<T>(
                    CyclotomicPolynomial(dd), elts[i]));
              }
            }
          }
          // every pair allowed, two kernels of one element included
          std::vector<size_t> owner(cand.size());
          for (size_t a = 0; a < cand.size(); a++) {
            owner[a] = a;
          }
          std::optional<MyMatrix<T>> Xq = FirstKernelOfDimension<T, Tint>(
              blk.conditions, cand, owner, n, dimF);
          if (Xq) {
            MyMatrix<T> const &Xb = *Xq;
            MyMatrix<T> GramX = Xb * A * Xb.transpose();
            MyMatrix<T> GramXinv = Inverse(GramX);
            // e g e on X (rows), the orthogonal projection onto F
            std::vector<MyMatrix<T>> DspanX;
            for (auto const &x : elts) {
              DspanX.push_back(MyMatrix<T>(Xb * x * A * Xb.transpose() * GramXinv)
                                   .transpose());
            }
            comp.precomputed_mass = quat_mass(GramX, DspanX, qbounds);
            std::vector<MyMatrix<T>> OmegaQ;
            for (auto const &W2 : QuatIntegralPoints<T, Tint>(DspanX, dimF)) {
              OmegaQ.push_back(W2.transpose());
            }
            auto [ClassSumq, Nexpq, Hgensq] = center_field();
            (void)Nexpq;
            (void)Hgensq;
            MultiplicityFactorBounds mbounds;
            if (!type2_split) {
              comp.multiplicity_factor = MultiplicityLocalFactor<T, Tint>(
                  A, Xb, blk.E, elts, ListP, OmegaQ, ClassSumq, true, mbounds,
                  os);
            }
          }
#ifdef DEBUG_ISOMETRY_MASS
          os << "ISOM_MASS: quaternionic with s=" << sq << ", X "
             << (Xq ? "found" : "not found") << ", factor "
             << (comp.multiplicity_factor ? "found" : "not found") << "\n";
#endif
        }
#ifdef DEBUG_ISOMETRY_MASS
        os << "ISOM_MASS: quaternionic component of rank " << k << ", r=" << rq
           << ", center of degree " << c << ", mass "
           << (comp.precomputed_mass ? "found" : "not found") << "\n";
#endif
        if (!type2_split) {
          res.push_back(comp);
          continue;
        }
        comp.kind = ComponentKind::Unsupported;
        comp.precomputed_mass.reset();
        comp.multiplicity_factor.reset();
      }
    }
    if (central_q && m * m == e) {
      // a common eigenspace of dimension m, for one or two elements
      std::optional<MyMatrix<T>> found;
      MyMatrix<T> found_basis;
      // the common +-1-eigenspaces of one or two elements
      std::vector<MyMatrix<T>> cand;
      std::vector<size_t> owner;
      for (size_t i = 0; i < n_elt; i++) {
        for (int sgn : {1, -1}) {
          cand.push_back(elts[i] - UniversalScalarConversion<T, int64_t>(sgn) * IdentityMat<T>(n));
          owner.push_back(i);
        }
      }
      if (std::optional<MyMatrix<T>> X =
              FirstKernelOfDimension<T, Tint>(blk.conditions, cand, owner, n, m)) {
        found = MyMatrix<T>(*X * A * X->transpose());
        found_basis = *X;
      }
      if (!found) {
        // the kernels of Phi_d(g) for the d dividing the order of g, alone
        // or by pairs: any End_G-submodule of dimension m does (its
        // orthogonal projection is a primitive idempotent of the span of G)
        std::optional<MyMatrix<T>> Xq = CyclotomicKernelSubmodule<T, Tint>(
            elts, blk.conditions, n, m);
        if (Xq) {
          found = MyMatrix<T>(*Xq * A * Xq->transpose());
          found_basis = *Xq;
        }
      }
      if (!found) {
        // End_G x for a pure tensor x
        if (std::optional<MyMatrix<T>> X =
                EndomorphismSpanSubmodule<T, Tint>(blk.E, End, m)) {
          found = MyMatrix<T>(*X * A * X->transpose());
          found_basis = *X;
        }
      }
      if (!found && m == 1) {
        // multiplicity one: the centralizer on the component is {+-1},
        // which preserves every lattice; the mass is that of any lattice of
        // rank one (1/2), without the line F and without a local factor
        MyMatrix<T> one(1, 1);
        one(0, 0) = T(2);
        comp.kind = ComponentKind::Multiplicity;
        comp.multiplicity_gram = one;
        res.push_back(comp);
        continue;
      }
      if (found) {
        comp.kind = ComponentKind::Multiplicity;
        comp.multiplicity_gram = found;
        MultiplicityFactorBounds mbounds;
        comp.multiplicity_factor = MultiplicityLocalFactor<T, Tint>(
            A, found_basis, blk.E, elts, ListP, {}, {}, false, mbounds, os);
      }
    }
    // with a center of degree c, End_G = M_m(K) has dimension m^2 c
    if (c > 1 && e % c == 0) {
      m = static_cast<int>(std::llround(std::sqrt(static_cast<double>(e / c))));
    }
    if (c > 1 && m * m * c == e) {
      /*
        End_G = M_m(K) with K = Q(chi) the center, of degree c. The
        multiplicity space is X = L cap F, F the image of an idempotent of
        Q[H] of K-rank one on W, of dimension m c, with the action of the
        order generated by the class sums: Hermitian over K when K is CM
        (the adjoint acts nontrivially on the center), quadratic over K
        when K is totally real. K is the subfield of Q(zeta_N), N the
        exponent of G, fixed by the a with C_g = C_{g^a} on V_rho for every
        class (sigma_a sends the central character of C_g to that of
        C_{g^a}). F is stable under End_G, so it has no dimension m c when
        End_G is not split (M_{m/2}(D), D a quaternion algebra over K, has
        no module of that dimension), and the component is then left
        without a mass.
       */
      std::vector<MyMatrix<T>> const &Center = CenterBasis;
      MyMatrix<T> Gblk = blk.E * A * blk.E.transpose();
      MyMatrix<T> GblkInv = Inverse(Gblk);
      bool is_cm = false;
      for (auto const &Z : Center) {
        MyMatrix<T> Zs = (GblkInv * Z * Gblk).transpose();
        if (Zs != Z) {
          is_cm = true;
        }
      }
      auto [ClassSum, Nexp, Hgens] = center_field();
      std::optional<MyMatrix<T>> found_basis;
      if (Nexp > 2) {
        std::vector<MyMatrix<T>> cand;
        std::vector<size_t> owner;
        for (size_t i = 0; i < n_elt; i++) {
          for (int sgn : {1, -1}) {
            cand.push_back(elts[i] -
                           UniversalScalarConversion<T, int64_t>(sgn) * IdentityMat<T>(n));
            owner.push_back(i);
          }
        }
        found_basis = FirstKernelOfDimension<T, Tint>(blk.conditions, cand,
                                                      owner, n, m * c);
        if (!found_basis) {
          found_basis = CyclotomicKernelSubmodule<T, Tint>(
              elts, blk.conditions, n, m * c);
        }
      }
      if (!found_basis) {
        // End_G x for a pure tensor x
        found_basis =
            EndomorphismSpanSubmodule<T, Tint>(blk.E, End, m * c);
      }
#ifdef DEBUG_ISOMETRY_MASS
      os << "ISOM_MASS: center of degree " << c << " (" << (is_cm ? "CM" : "real")
         << "), m=" << m << ", Nexp=" << Nexp << ", F "
         << (found_basis ? "found" : "not found") << "\n";
#endif
      // the unitary group on the whole block, by the Tamagawa number
      auto unitary_block_mass = [&]() -> std::optional<T> {
        std::vector<MyMatrix<T>> EltsBlk;
        for (auto const &x : elts) {
          EltsBlk.push_back(action(x).transpose());
        }
        return UnitaryInnerFormMass<T, Tint>(
            Gblk, EltsBlk, MakeAbelianCMField(Nexp, Hgens), m,
            QuaternionicMassBounds(), os);
      };
      // the mass of the block by the Tamagawa number, printed for a
      // comparison: a second mass computation, expensive
#ifdef SANITY_CHECK_UNITARY_INNER_FORM_MASS_EXTENSIVE
      std::cerr << "ISOM_MASS: SANITY_CHECK_UNITARY_INNER_FORM_MASS_EXTENSIVE\n";
      if (found_basis && is_cm && m <= 2) {
        std::optional<T> um = unitary_block_mass();
        os << "ISOM_MASS: unitary mass of the block by the Tamagawa number "
           << (um ? *um : T(-1)) << "\n";
      }
      std::cerr << "ISOM_MASS: SANITY_CHECK_UNITARY_INNER_FORM_MASS_EXTENSIVE, done\n";
#endif
      if (!found_basis && is_cm && m == 2) {
        /*
          End_G = D, a quaternion division algebra over the CM center K (a
          complex character of Schur index 2): U(D, tau) is an inner form
          of U_2 over K^+. When G spans the commutant of D on the block
          (of dimension 4 c), X is the whole block.
         */
        int dimA = 0;
        {
          MyMatrix<T> Fl(static_cast<int>(n_elt), k * k);
          for (size_t i = 0; i < n_elt; i++) {
            MyMatrix<T> Ai = action(elts[i]);
            for (int u2 = 0; u2 < k; u2++) {
              for (int v2 = 0; v2 < k; v2++) {
                Fl(static_cast<int>(i), u2 * k + v2) = Ai(u2, v2);
              }
            }
          }
          dimA = RankMat(Fl);
        }
        if (dimA == 4 * c && k == 4 * c) {
          comp.kind = ComponentKind::Precomputed;
          comp.precomputed_mass = unitary_block_mass();
          comp.multiplicity_factor = T(1);
        }
#ifdef DEBUG_ISOMETRY_MASS
        os << "ISOM_MASS: End_G a quaternion division algebra over the CM "
           << "center, G spans an algebra of dimension " << dimA << ", mass "
           << (comp.precomputed_mass ? "found" : "not found") << "\n";
#endif
      }
      if (found_basis) {
        MyMatrix<T> const &Xb = *found_basis;
        int mX = Xb.rows();
        MyMatrix<T> XX = Xb * Xb.transpose();
        MyMatrix<T> XXinv = Inverse(XX);
        // the class sums on X, then the ring they generate
        auto on_X = [&](MyMatrix<T> const &Mc) -> MyMatrix<T> {
          return Xb * Mc * Xb.transpose() * XXinv;
        };
        std::vector<MyMatrix<T>> gens{IdentityMat<T>(mX)};
        for (auto const &Cs : ClassSum) {
          gens.push_back(on_X(Cs));
        }
        auto zbasis = [&](std::vector<MyMatrix<T>> const &L)
            -> std::vector<MyMatrix<T>> {
          MyMatrix<Tint> F(static_cast<int>(L.size()), mX * mX);
          for (size_t i = 0; i < L.size(); i++) {
            for (int u = 0; u < mX; u++) {
              for (int v = 0; v < mX; v++) {
                F(static_cast<int>(i), u * mX + v) =
                    UniversalScalarConversion<Tint, T>(L[i](u, v));
              }
            }
          }
          MyMatrix<Tint> Bz = GetZbasis(F);
          std::vector<MyMatrix<T>> res_b;
          for (int i = 0; i < Bz.rows(); i++) {
            MyMatrix<T> Y(mX, mX);
            for (int u = 0; u < mX; u++) {
              for (int v = 0; v < mX; v++) {
                Y(u, v) = UniversalScalarConversion<T, Tint>(Bz(i, u * mX + v));
              }
            }
            res_b.push_back(Y);
          }
          return res_b;
        };
        std::vector<MyMatrix<T>> Omega = zbasis(gens);
        for (int iter = 0; iter < 8; iter++) {
          std::vector<MyMatrix<T>> prods = Omega;
          for (auto const &X1 : Omega) {
            for (auto const &X2 : Omega) {
              prods.push_back(X1 * X2);
            }
          }
          std::vector<MyMatrix<T>> next = zbasis(prods);
          if (next.size() == Omega.size()) {
            Omega = next;
            break;
          }
          Omega = next;
        }
        // the multiplier ring of X in K (the integral matrices of the span
        // of the class sums on X), which can be larger than the order the
        // class sums generate: X is an O_K-lattice when it is O_K
        if (static_cast<int>(Omega.size()) == c) {
          std::vector<MyMatrix<T>> Mult = QuatIntegralPoints<T, Tint>(Omega, mX);
          if (static_cast<int>(Mult.size()) == c) {
            Omega = Mult;
          }
        }
        AbelianField Fk = is_cm ? MakeAbelianCMField(Nexp, Hgens)
                                  : MakeAbelianRealField(Nexp, Hgens);
        // the order is O_K when its discriminant is that of K (the product
        // of the conductors of the characters)
        bool maximal = (static_cast<int>(Omega.size()) == c);
        if (maximal) {
          MyMatrix<T> Tr(c, c);
          for (int i = 0; i < c; i++) {
            for (int j = 0; j < c; j++) {
              MyMatrix<T> P = Omega[i] * Omega[j];
              T t(0);
              for (int u = 0; u < mX; u++) {
                t += P(u, u);
              }
              Tr(i, j) = t / UniversalScalarConversion<T, int64_t>(m);
            }
          }
          T dO = T_abs(DeterminantMat(Tr));
          DirichletGroup Gd = ComputeDirichletGroup(Fk.n);
          T dK(1);
          for (auto const &t : FieldCharacters(Gd, Fk.H, false)) {
            dK *= UniversalScalarConversion<T, int64_t>(PrimitiveCharacter(Gd, t).first);
          }
          maximal = (dO == dK);
        }
#ifdef DEBUG_ISOMETRY_MASS
        os << "ISOM_MASS: the order of the class sums on X has rank "
           << Omega.size() << ", maximal " << maximal << "\n";
#endif
        if (maximal) {
          comp.kind = is_cm ? ComponentKind::AbelianHermitian
                            : ComponentKind::AbelianQuadratic;
          comp.center_field = Fk;
          comp.center_gram = Xb * A * Xb.transpose();
          comp.center_omega = Omega;
          MultiplicityFactorBounds mbounds;
          comp.multiplicity_factor = MultiplicityLocalFactor<T, Tint>(
              A, Xb, blk.E, elts, ListP, Omega, ClassSum, false, mbounds, os);
        } else if (static_cast<int>(Omega.size()) == c) {
          /*
            The order O of X in K (its multiplier ring) is not O_K, so X is
            not an O_K-lattice. O_K is found by adding x = y / p, y in O, x
            integral (K commutative: the characteristic polynomial of x on
            F, by the traces of its powers, is integral), at the primes of
            [O_K : O], until the discriminant is d_K. Then S = O_K X is an
            O_K-lattice in the same space, of mass given by the Hermitian
            (or quadratic) mass, and mass(X) = mass(S) prod_p vol U(S_p) /
            vol U(X_p) over the primes of [S : X], the volumes from the
            orders of the images modulo p^c and the Lie lattices
            (MassQuaternionic.h, which holds for any order acting).
           */
          comp.kind = ComponentKind::Precomputed;
          auto zbasis_rat = [&](std::vector<MyMatrix<T>> const &L)
              -> std::vector<MyMatrix<T>> {
            MyMatrix<T> F(static_cast<int>(L.size()), mX * mX);
            for (size_t i = 0; i < L.size(); i++) {
              for (int u = 0; u < mX; u++) {
                for (int v = 0; v < mX; v++) {
                  F(static_cast<int>(i), u * mX + v) = L[i](u, v);
                }
              }
            }
            MyMatrix<T> Bz = GetZbasis(F);
            std::vector<MyMatrix<T>> res_b;
            for (int i = 0; i < Bz.rows(); i++) {
              MyMatrix<T> Y(mX, mX);
              for (int u = 0; u < mX; u++) {
                for (int v = 0; v < mX; v++) {
                  Y(u, v) = Bz(i, u * mX + v);
                }
              }
              res_b.push_back(Y);
            }
            return res_b;
          };
          T dKval(1);
          {
            DirichletGroup Gd = ComputeDirichletGroup(Fk.n);
            for (auto const &t : FieldCharacters(Gd, Fk.H, false)) {
              dKval *= UniversalScalarConversion<T, int64_t>(PrimitiveCharacter(Gd, t).first);
            }
          }
          T tr_scale = UniversalScalarConversion<T, int64_t>(mX / c);
          auto disc_of = [&](std::vector<MyMatrix<T>> const &L) -> T {
            int d0 = static_cast<int>(L.size());
            MyMatrix<T> G(d0, d0);
            for (int i = 0; i < d0; i++) {
              for (int j = 0; j < d0; j++) {
                G(i, j) = QuatTrace(MyMatrix<T>(L[i] * L[j])) / tr_scale;
              }
            }
            return T_abs(DeterminantMat(G));
          };
          // x integral over Z: the elementary symmetric functions of its
          // eigenvalues, from the traces of its powers, are integers
          auto integral = [&](MyMatrix<T> const &x) -> bool {
            std::vector<T> pk(mX + 1), ek(mX + 1);
            MyMatrix<T> pw = IdentityMat<T>(mX);
            for (int kk = 1; kk <= mX; kk++) {
              pw = pw * x;
              pk[kk] = QuatTrace(pw);
            }
            ek[0] = T(1);
            for (int kk = 1; kk <= mX; kk++) {
              T sacc(0);
              for (int i = 1; i <= kk; i++) {
                T term = ek[kk - i] * pk[i];
                sacc += (i % 2 == 1) ? term : T(-term);
              }
              ek[kk] = sacc / UniversalScalarConversion<T, int64_t>(kk);
              if (!IsInteger(ek[kk])) {
                return false;
              }
            }
            return true;
          };
          auto ring_of = [&](std::vector<MyMatrix<T>> L) -> std::vector<MyMatrix<T>> {
            for (int iter = 0; iter < 20; iter++) {
              std::vector<MyMatrix<T>> prods = L;
              for (auto const &X1 : L) {
                for (auto const &X2 : L) {
                  prods.push_back(X1 * X2);
                }
              }
              std::vector<MyMatrix<T>> next = zbasis_rat(prods);
              if (next.size() == L.size() && disc_of(next) == disc_of(L)) {
                return next;
              }
              L = next;
            }
            return L;
          };
          std::vector<MyMatrix<T>> OK = Omega;
          bool ok_found = false;
          for (int iter = 0; iter < 64; iter++) {
            T q = disc_of(OK) / dKval;
            if (q == T(1)) {
              ok_found = true;
              break;
            }
            if (!IsInteger(q)) {
              break;
            }
            std::vector<Tint> pr = FactorsInt(UniversalScalarConversion<Tint, T>(q));
            std::sort(pr.begin(), pr.end());
            pr.erase(std::unique(pr.begin(), pr.end()), pr.end());
            bool enlarged = false;
            for (Tint const &p_i : pr) {
              int64_t pp = UniversalScalarConversion<int64_t, Tint>(p_i);
              size_t n_el = MyPow<size_t>(static_cast<size_t>(pp), c);
              for (size_t sidx = 1; sidx < n_el && !enlarged; sidx++) {
                size_t tt = sidx;
                MyMatrix<T> x = ZeroMatrix<T>(mX, mX);
                for (int i = 0; i < c; i++) {
                  int64_t a = static_cast<int64_t>(tt % static_cast<size_t>(pp));
                  tt /= static_cast<size_t>(pp);
                  if (a != 0) {
                    x += UniversalScalarConversion<T, int64_t>(a) * OK[i];
                  }
                }
                x /= UniversalScalarConversion<T, int64_t>(pp);
                if (integral(x)) {
                  std::vector<MyMatrix<T>> L = OK;
                  L.push_back(x);
                  OK = ring_of(zbasis_rat(L));
                  enlarged = true;
                }
              }
              if (enlarged) {
                break;
              }
            }
            if (!enlarged) {
              break;
            }
          }
#ifdef DEBUG_ISOMETRY_MASS
          os << "ISOM_MASS: O_K on X " << (ok_found ? "found" : "not found")
             << ", [O_K : O]^2 = " << disc_of(Omega) / dKval << "\n";
#endif
          if (ok_found && static_cast<int>(OK.size()) == c) {
            MyMatrix<T> GramX = Xb * A * Xb.transpose();
            // S = O_K X, rows in the coordinates of X
            std::vector<MyVector<T>> rows;
            for (auto const &w : OK) {
              for (int i = 0; i < mX; i++) {
                rows.push_back(GetMatrixRow(w, i));
              }
            }
            MyMatrix<T> Rall(static_cast<int>(rows.size()), mX);
            for (size_t i = 0; i < rows.size(); i++) {
              for (int j = 0; j < mX; j++) {
                Rall(static_cast<int>(i), j) = rows[i](j);
              }
            }
            MyMatrix<T> Sb = GetZbasis(Rall);
            MyMatrix<T> SbInv = Inverse(Sb);
            MyMatrix<T> GramS = Sb * GramX * Sb.transpose();
            // integral by a scaling, which changes neither U nor the mass
            T den = GetDenominatorMatrix(GramS);
            GramS *= den;
            std::vector<MyMatrix<T>> OmS;
            for (auto const &w : OK) {
              OmS.push_back(MyMatrix<T>(Sb * w * SbInv));
            }
            std::optional<T> mass_S =
                is_cm ? HermitianMassAbelianLattice(GramS, OmS, Fk)
                      : QuadraticMassAbelianLattice(GramS, OmS, Fk,
                                                    QuadLValueMethod::Mpfr);
            T index_SX = T_abs(DeterminantMat(SbInv));
            std::optional<T> mass_X;
            if (mass_S) {
              mass_X = *mass_S;
              std::vector<MyMatrix<T>> gX, gS;
              for (auto const &w : Omega) {
                gX.push_back(MyMatrix<T>(w.transpose()));
              }
              for (auto const &w : OmS) {
                gS.push_back(MyMatrix<T>(w.transpose()));
              }
              std::vector<MyMatrix<T>> uX =
                  QuatIntegralCommutant<T, Tint>(GramX, gX, true);
              std::vector<MyMatrix<T>> uS =
                  QuatIntegralCommutant<T, Tint>(GramS, gS, true);
              std::vector<MyMatrix<T>> LamX =
                  QuatIntegralCommutant<T, Tint>(GramX, gX, false);
              std::vector<MyMatrix<T>> LamS =
                  QuatIntegralCommutant<T, Tint>(GramS, gS, false);
              T discX = QuatTraceDiscriminant(uX);
              T discS = QuatTraceDiscriminant(uS);
              // the integral forms B(w x, y), w in the order
              auto forms = [&](MyMatrix<T> const &Bz,
                               std::vector<MyMatrix<T>> const &W)
                  -> std::vector<MyMatrix<T>> {
                std::vector<MyMatrix<T>> L;
                for (auto const &w : W) {
                  L.push_back(MyMatrix<T>(w.transpose() * Bz));
                }
                return QuatIntegralPoints<T, Tint>(L, Bz.rows());
              };
              std::vector<MyMatrix<T>> FormsX = forms(GramX, gX);
              std::vector<MyMatrix<T>> FormsS = forms(GramS, gS);
              MultiplicityFactorBounds ib;
              std::vector<Tint> pr_SX =
                  FactorsInt(UniversalScalarConversion<Tint, T>(index_SX));
              std::sort(pr_SX.begin(), pr_SX.end());
              pr_SX.erase(std::unique(pr_SX.begin(), pr_SX.end()), pr_SX.end());
              for (Tint const &p_i : pr_SX) {
                int64_t pp = UniversalScalarConversion<int64_t, Tint>(p_i);
                int cc = (pp == 2) ? 2 : 1;
                auto nX = QuatLocalImageOrders<T, Tint>(GramX, gX, LamX, FormsX,
                                                        pp, cc, ib, os);
                auto nS = QuatLocalImageOrders<T, Tint>(GramS, gS, LamS, FormsS,
                                                        pp, cc, ib, os);
                if (!nX || !nS) {
                  mass_X = {};
                  break;
                }
                int vX = RationalValuation(discX, pp);
                int vS = RationalValuation(discS, pp);
#ifdef SANITY_CHECK_ISOMETRY_MASS
                if ((vX - vS) % 2 != 0) {
                  std::cerr << "ISOM_MASS: the Lie discriminants differ by an "
                               "odd power of " << pp << "\n";
                  throw TerminalException{1};
                }
#endif
                T ratio = UniversalScalarConversion<T, Tint>(nS->first) /
                          UniversalScalarConversion<T, Tint>(nX->first);
                int ex = (vX - vS) / 2;
                for (int i = 0; i < (ex >= 0 ? ex : -ex); i++) {
                  if (ex >= 0) {
                    ratio *= UniversalScalarConversion<T, int64_t>(pp);
                  } else {
                    ratio /= UniversalScalarConversion<T, int64_t>(pp);
                  }
                }
#ifdef DEBUG_ISOMETRY_MASS
                os << "ISOM_MASS: p=" << pp << " |U(X_p)|=" << nX->first
                   << " |U(S_p)|=" << nS->first << " ratio=" << ratio << "\n";
#endif
                *mass_X *= ratio;
              }
            }
            comp.precomputed_mass = mass_X;
            if (mass_X) {
              MultiplicityFactorBounds mbounds;
              comp.multiplicity_factor = MultiplicityLocalFactor<T, Tint>(
                  A, Xb, blk.E, elts, ListP, Omega, ClassSum, true, mbounds, os);
            }
          }
        }
      }
    }
#ifdef DEBUG_ISOMETRY_MASS
    os << "ISOM_MASS: non-cyclic component of rank " << k << ", dim End_G="
       << e << ", multiplicity lattice "
       << (comp.multiplicity_gram ? "found" : "not found") << "\n";
#endif
    if (comp.kind == ComponentKind::Unsupported) {
      std::cerr << "ISOM_MASS: a component of rank " << k << " with End_G of "
                << "dimension " << e << " and a center of degree " << c
                << " has no mass method: ";
      if (c == 1 && m * m == e) {
        std::cerr << "End_G = M_" << m << "(Q), no End_G-submodule F of "
                  << "dimension " << m << " found by the searches";
      } else if (c > 1 && m * m * c == e) {
        std::cerr << "End_G of dimension m^2 c, m=" << m << ", no "
                  << "End_G-submodule F of dimension m c found by the "
                  << "searches, or End_G not split";
      } else {
        std::cerr << "End_G not recognized";
      }
      std::cerr << ", please insert the code\n";
      throw TerminalException{1};
    }
    res.push_back(comp);
  }
  return res;
}

// clang-format off
#endif  // SRC_GENERA_MASSGROUPCOMPONENTS_H_
// clang-format on
