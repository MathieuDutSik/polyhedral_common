// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_DELAUNAY_IMINIMUM_H_
#define SRC_DELAUNAY_IMINIMUM_H_

// clang-format off
#include "LatticeDelaunay.h"
#include "Enumeration_k_space.h"
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_IMINIMUM
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_IMINIMUM
#endif

/*
  Computation of the i-minimum mu_i(L)^2 of a lattice L given by its
  Gram matrix: the supremum over the i-dimensional affine subspaces X
  of the squared distance d(X, L)^2. For a fixed rational direction
  V the supremum over the affine subspaces of direction V is the
  covering radius of the projection of L onto the orthogonal
  complement of V, and the supremum over all directions is attained
  by a rational direction.
  ---
  The value is obtained with a certificate by the procedure
  BOUND(L', i, T) which certifies mu(pi_{V^perp}(L'))^2 <= T for
  every rational direction V of dimension i of a section L', up to a
  finite list of explicitly returned exceptional sublattices L' cap V:
  * (base) if dim L' = i + 1 the projection along V is one
    dimensional of squared covolume det(L')/det(L' cap V), so the
    value is det(L')/(4 det(L' cap V)) exactly and the exceptions are
    the primitive rank i sublattices of determinant below
    det(L')/(4T).
  * (free directions) by the transference bound of Banaszczyk the
    directions not orthogonal to any primitive dual vector of norm
    below (dim L' - i)^2/(4T) are certified.
  * (trapped directions) the directions orthogonal to a subspace W
    spanned by short dual vectors satisfy the projection inequality
    mu(pi_{V^perp}(L'))^2 <= mu(pi_W(L'))^2 + mu(pi'(L' cap W^perp))^2
    where the label loss mu(pi_W(L'))^2 is computed exactly, so it
    suffices to run BOUND(L' cap W^perp, i, T - loss). A worklist of
    subspaces W is processed; when a sub-certification fails the
    failing constraint is joined to W and the enlarged subspace is
    put back on the worklist.
  The driver enumerates the orbits of the rank i sublattices of
  determinant at most M, sets R^2 to the maximum of the covering
  radii of the projections, runs BOUND(L, i, R^2) and disposes of
  the exceptions: an exceptional determinant at most M is covered by
  the enumeration, and an exception beyond M has its exact value
  computed directly; if a value exceeds R^2 the bound M is raised to
  reach that determinant and the process restarts with the improved
  lower bound. On success mu_i(L)^2 = R^2 with certificate.
  ---
  Compared with the prototype implementation certify_mu1.py
  distributed with the paper, the label losses are computed exactly
  in every dimension (through the Delaunay tessellation of the label
  lattice) and at every level of the recursion the trap subspaces
  are reduced modulo the automorphisms of the full lattice that
  preserve the chain of traps leading to the section, which shrinks
  the worklists by orders of magnitude on the very symmetric
  lattices. The reduction is tested on chain-colored invariant
  vector configurations through the ListMat_Vdiag machinery, and it
  is sound because the reducing automorphisms extend to the full
  lattice: they transport the certification statements, the failure
  joins and the values of the exceptional directions alike.
 */

// The Gram matrix of the projection of the lattice onto the orthogonal
// complement of the span of SubBasis. SubBasis must be saturated so that
// its completion to a basis of Z^n projects to a basis of the projected
// lattice.
template <typename T, typename Tint>
MyMatrix<T> GetProjectedGramMatrix(MyMatrix<T> const &GramMat,
                                   MyMatrix<Tint> const &SubBasis) {
  int n = GramMat.rows();
  MyMatrix<Tint> TheCompl = SubspaceCompletionInt(SubBasis, n);
  MyMatrix<T> TheCompl_T = UniversalMatrixConversion<T, Tint>(TheCompl);
  MyMatrix<T> TheProj = GetOrthogonalProjector(GramMat, SubBasis);
  // For a row vector x the projection onto span(SubBasis) is x P^T,
  // so the projection onto the orthogonal complement is x (I - P^T).
  MyMatrix<T> ProjCompl = TheCompl_T - TheCompl_T * TheProj.transpose();
  MyMatrix<T> RedGram = ProjCompl * GramMat * ProjCompl.transpose();
#ifdef SANITY_CHECK_IMINIMUM
  MyMatrix<T> SubBasis_T = UniversalMatrixConversion<T, Tint>(SubBasis);
  MyMatrix<T> eGram = SubBasis_T * GramMat * SubBasis_T.transpose();
  T det_prod = DeterminantMat(RedGram) * DeterminantMat(eGram);
  T det_gram = DeterminantMat(GramMat);
  if (det_prod != det_gram) {
    std::cerr << "GetProjectedGramMatrix: det(RedGram) * det(eGram) = "
              << det_prod << " but det(GramMat) = " << det_gram << "\n";
    throw TerminalException{1};
  }
#endif
  return RedGram;
}

// The exact squared covering radius of a two dimensional lattice:
// Gauss reduction of the Gram matrix to 0 <= 2b <= a <= c makes the
// fundamental triangle of squared sides (a, c, a + c - 2b) non
// obtuse, and the covering radius is its circumradius.
template <typename T, typename Tint>
T CoveringRadiusSquaredDim2(MyMatrix<T> const &G) {
  T a = G(0, 0);
  T b = G(0, 1);
  T c = G(1, 1);
  T two(2);
  while (true) {
    if (a > c) {
      std::swap(a, c);
    }
    Tint k = UniversalNearestScalarInteger<Tint, T>(b / a);
    T k_T = UniversalScalarConversion<T, Tint>(k);
    c = c - two * k_T * b + k_T * k_T * a;
    b = b - k_T * a;
    if (a <= c && two * T_abs(b) <= a) {
      break;
    }
  }
  b = T_abs(b);
  T third = a + c - two * b;
  T den = T(4) * (a * c - b * b);
  return a * c * third / den;
}

// The exact squared covering radius of the label lattice pi_W(L) for
// a subspace W spanned by dual vectors. The input is the Gram matrix
// of L* cap W, whose dual is the label lattice. Closed forms in
// dimensions one and two, the Delaunay tessellation beyond.
template <typename T, typename Tint, typename Tgroup>
T LabelLatticeCovSqr(MyMatrix<T> const &GramDualSection, std::ostream &os) {
  int k = GramDualSection.rows();
  if (k == 1) {
    return T(1) / (T(4) * GramDualSection(0, 0));
  }
  MyMatrix<T> LabelGram = Inverse(GramDualSection);
  if (k == 2) {
    return CoveringRadiusSquaredDim2<T, Tint>(LabelGram);
  }
  return ComputeCoveringRadiusSquared<T, Tint, Tgroup>(LabelGram, os);
}

// The saturated sublattice basis of the subspace spanned by the rows
// of Y (the rows may be dependent): the integral kernel of the
// integral kernel.
template <typename Tint>
MyMatrix<Tint> SaturationSpanOfRows(MyMatrix<Tint> const &Y) {
  MyMatrix<Tint> N = NullspaceIntTrMat(Y);
  if (N.rows() == 0) {
    return IdentityMat<Tint>(Y.cols());
  }
  return NullspaceIntTrMat(N);
}

// Scales a nonzero rational row to the primitive integral row with
// the same span.
template <typename T, typename Tint>
MyVector<Tint> PrimitiveIntegralRow(MyVector<T> const &v) {
  MyVector<T> vRed = RemoveFractionVector(v);
  return UniversalVectorConversion<Tint, T>(vRed);
}

template <typename Tint> struct IminimumFailTuple {
  // The constraints in the dual coordinates of the section the
  // failure is reported from; an empty list means that every
  // direction of the section is uncertified.
  std::vector<MyVector<Tint>> constraints;
};

template <typename Tint> struct IminimumBoundResult {
  std::vector<IminimumFailTuple<Tint>> fails;
};

// The chain-colored configuration used to reduce the traps. The rows
// are the invariant family of the full lattice (color 0) followed by
// the invariant families of the nested section sublattices of the
// chain of traps (colors 1, 2, ...). Two candidate continuations of
// the same chain are interchangeable exactly when an automorphism of
// the full lattice preserves every chain member and maps one to the
// other. Such an automorphism transports every certification
// statement, the values of the exceptional directions in the full
// lattice included, which makes the reduction sound.
template <typename T> struct IminimumChainConf {
  MyMatrix<T> conf;
  std::vector<T> Vdiag;
};

// The state of one certification run of BOUND(L, i, R^2).
template <typename T, typename Tint, typename Tgroup> struct IminimumBound {
  MyMatrix<T> GramMat;
  int i;
  std::ostream &os;
  std::vector<MyMatrix<T>> ListMat;
  IminimumChainConf<T> top_conf;
  // The exceptions, as canonicalized sublattices of the full lattice.
  std::unordered_set<MyMatrix<Tint>> exceptions;
  // Safety bound on the total number of traps processed.
  size_t n_trap = 0;
  size_t max_trap;
  bool aborted = false;
#ifdef DEBUG_IMINIMUM
  size_t n_aut_skip = 0;
#endif
  IminimumBound(MyMatrix<T> const &_GramMat, int const &_i,
                size_t const &_max_trap, std::ostream &_os)
      : GramMat(_GramMat), i(_i), os(_os), ListMat{_GramMat},
        max_trap(_max_trap) {
    MyMatrix<Tint> W =
        ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat, os);
    MyMatrix<T> W_T = UniversalMatrixConversion<T, Tint>(W);
    std::vector<T> Vdiag(W_T.rows(), T(0));
    top_conf = {std::move(W_T), std::move(Vdiag)};
  }

  // The chain configuration extended by the invariant family of the
  // sublattice K of the full lattice, with the color eColor.
  IminimumChainConf<T> extend_chain_conf(IminimumChainConf<T> const &prefix,
                                         MyMatrix<Tint> const &K,
                                         T const &eColor) {
    MyMatrix<T> K_T = UniversalMatrixConversion<T, Tint>(K);
    MyMatrix<T> eGram = K_T * GramMat * K_T.transpose();
    MyMatrix<Tint> SHV =
        ExtractInvariantVectorFamilyZbasis<T, Tint>(eGram, os);
    MyMatrix<Tint> invar = SHV * K;
    MyMatrix<T> invar_T = UniversalMatrixConversion<T, Tint>(invar);
    int n_prev = prefix.conf.rows();
    int n_new = invar_T.rows();
    MyMatrix<T> conf(n_prev + n_new, prefix.conf.cols());
    std::vector<T> Vdiag(n_prev + n_new);
    for (int u = 0; u < n_prev; u++) {
      conf.row(u) = prefix.conf.row(u);
      Vdiag[u] = prefix.Vdiag[u];
    }
    for (int u = 0; u < n_new; u++) {
      conf.row(n_prev + u) = invar_T.row(u);
      Vdiag[n_prev + u] = eColor;
    }
    return {std::move(conf), std::move(Vdiag)};
  }

  // The generators of the stabilizer of the chain, acting on the
  // dual coordinates of the section S by right multiplication: for a
  // generator g of the automorphisms of the chain-colored
  // configuration, the induced action on the section coordinates is
  // the integral matrix A with S g = A S, and the action on the dual
  // coordinates is (A^{-1})^T.
  std::vector<MyMatrix<Tint>>
  chain_stab_dual_gens(IminimumChainConf<T> const &prefix,
                       [[maybe_unused]] MyMatrix<Tint> const &S,
                       MyMatrix<T> const &S_T) {
    std::vector<MyMatrix<T>> LGen = GetIntAutomorphism_ListMat_Vdiag<T, Tgroup>(
        prefix.conf, ListMat, prefix.Vdiag, os);
    MyMatrix<T> Fac = S_T.transpose() * Inverse(MyMatrix<T>(
        S_T * S_T.transpose()));
    std::vector<MyMatrix<Tint>> l_dgen;
    for (auto &g_T : LGen) {
      MyMatrix<T> A_T = S_T * g_T * Fac;
      MyMatrix<Tint> A = UniversalMatrixConversion<Tint, T>(A_T);
#ifdef SANITY_CHECK_IMINIMUM
      MyMatrix<Tint> g = UniversalMatrixConversion<Tint, T>(g_T);
      if (A * S != S * g) {
        std::cerr << "chain_stab_dual_gens: the generator does not "
                     "stabilize the section\n";
        throw TerminalException{1};
      }
#endif
      MyMatrix<T> B_T = Inverse(A_T).transpose();
      MyMatrix<Tint> B = UniversalMatrixConversion<Tint, T>(B_T);
      l_dgen.push_back(B);
    }
    return l_dgen;
  }

  // Certifies mu(pi_{V^perp})^2 <= target for the rational directions
  // V of dimension i of the section with basis rows S (a saturated
  // sublattice of the full lattice, in canonical form), up to the
  // exceptions collected in this->exceptions and the returned failure
  // tuples. The prefix is the chain configuration of the traps
  // leading to S and depth is the length of that chain.
  IminimumBoundResult<Tint> bound_rec(MyMatrix<Tint> const &S, T const &target,
                                      IminimumChainConf<T> const &prefix,
                                      int const &depth) {
    int r = S.rows();
    IminimumBoundResult<Tint> result;
    if (r <= i || aborted) {
      return result;
    }
    if (target <= 0) {
      result.fails.push_back(IminimumFailTuple<Tint>{});
      return result;
    }
    MyMatrix<T> S_T = UniversalMatrixConversion<T, Tint>(S);
    MyMatrix<T> G_S = S_T * GramMat * S_T.transpose();
    if (r == i + 1) {
      // Base case: the value along V is det(S)/(4 det(S cap V))
      // exactly, so the exceptions are the primitive rank i
      // sublattices of determinant strictly below det(S)/(4 target).
      T detS = DeterminantMat(G_S);
      T MaxDet = detS / (T(4) * target);
      std::vector<MyMatrix<Tint>> l_sub =
          Rankin_k_level<T, Tint>(G_S, i, MaxDet, os);
      for (auto &X : l_sub) {
        MyMatrix<T> X_T = UniversalMatrixConversion<T, Tint>(X);
        T det_X = DeterminantMat(MyMatrix<T>(X_T * G_S * X_T.transpose()));
        if (T(4) * target * det_X < detS) {
          MyMatrix<Tint> Xfull = CanonicalizeSublatticeBasis(
              MyMatrix<Tint>(X * S));
          exceptions.insert(Xfull);
        }
      }
      return result;
    }
    int m = r - i;
    T s2 = T(m * m) / (T(4) * target);
    MyMatrix<T> G_S_inv = Inverse(G_S);
    std::vector<MyVector<Tint>> l_short =
        computeLevel_GramMat<T, Tint>(G_S_inv, s2, os);
    std::vector<std::vector<MyVector<Tint>>> work;
    for (auto &y : l_short) {
      if (!IsVectorPrimitive(y)) {
        continue;
      }
      MyVector<T> y_T = UniversalVectorConversion<T, Tint>(y);
      T norm = y_T.dot(G_S_inv * y_T);
      if (norm < s2) {
        work.push_back({y});
      }
    }
#ifdef DEBUG_IMINIMUM
    os << "IMINIMUM: bound_rec r=" << r << " target=" << target
       << " |short duals|=" << work.size() << "\n";
#endif
    std::unordered_set<MyMatrix<Tint>> seen;
    // The generators of the stabilizer of the chain, computed once:
    // each processed trap has its whole orbit marked as seen, so the
    // equivalent traps are skipped by the plain lookup below.
    std::vector<MyMatrix<Tint>> l_dgen;
    if (!work.empty()) {
      l_dgen = chain_stab_dual_gens(prefix, S, S_T);
    }
    T eColor(depth + 1);
    while (!work.empty()) {
      std::vector<MyVector<Tint>> tup = work.back();
      work.pop_back();
      MyMatrix<Tint> Y = MatrixFromVectorFamily(tup);
      MyMatrix<Tint> Ys = CanonicalizeSublatticeBasis(
          SaturationSpanOfRows(Y));
      if (seen.count(Ys) == 1) {
        continue;
      }
      std::vector<MyMatrix<Tint>> orbit = OrbitSublattice(l_dgen, Ys);
      for (auto &oYs : orbit) {
        seen.insert(oYs);
      }
#ifdef DEBUG_IMINIMUM
      n_aut_skip += orbit.size() - 1;
#endif
      int kW = Ys.rows();
      if (kW >= r) {
        result.fails.push_back({tup});
        continue;
      }
      n_trap += 1;
      if (n_trap > max_trap) {
        aborted = true;
        break;
      }
      MyMatrix<T> Ys_T = UniversalMatrixConversion<T, Tint>(Ys);
      MyMatrix<T> GramDualSection = Ys_T * G_S_inv * Ys_T.transpose();
      T loss = LabelLatticeCovSqr<T, Tint, Tgroup>(GramDualSection, os);
      T sub_target = target - loss;
#ifdef DEBUG_IMINIMUM
      os << "IMINIMUM: trap dim=" << kW << " loss=" << loss
         << " sub_target=" << sub_target << "\n";
#endif
      if (sub_target <= 0) {
        result.fails.push_back({tup});
        continue;
      }
      MyMatrix<Tint> K = CanonicalizeSublatticeBasis(NullspaceIntTrMat(Ys));
      MyMatrix<Tint> Ssub = CanonicalizeSublatticeBasis(
          MyMatrix<Tint>(K * S));
      IminimumBoundResult<Tint> sub;
      if (r - kW > i + 1) {
        // The sub call will process its own traps, which needs the
        // extended chain configuration.
        IminimumChainConf<T> child_conf =
            extend_chain_conf(prefix, Ssub, eColor);
        sub = bound_rec(Ssub, sub_target, child_conf, depth + 1);
      } else {
        // Base case or trivial call: the chain configuration is not
        // used.
        sub = bound_rec(Ssub, sub_target, prefix, depth + 1);
      }
      if (sub.fails.empty()) {
        continue;
      }
      // Conversion of the failing constraints from the dual
      // coordinates of the sub section to the dual coordinates of
      // the present section: through the common geometric
      // coordinates, y'' G_sub^{-1} Ssub for the vector and
      // c G S^T for the dual coordinates.
      MyMatrix<T> Ssub_T = UniversalMatrixConversion<T, Tint>(Ssub);
      MyMatrix<T> G_sub = Ssub_T * GramMat * Ssub_T.transpose();
      MyMatrix<T> Conv = Inverse(G_sub) * Ssub_T * GramMat * S_T.transpose();
      for (auto &ftup : sub.fails) {
        if (ftup.constraints.size() == 0) {
          result.fails.push_back({tup});
        } else {
          std::vector<MyVector<Tint>> newtup = tup;
          for (auto &ycon : ftup.constraints) {
            MyVector<T> ycon_T = UniversalVectorConversion<T, Tint>(ycon);
            MyVector<T> z = Conv.transpose() * ycon_T;
            newtup.push_back(PrimitiveIntegralRow<T, Tint>(z));
          }
          work.push_back(newtup);
        }
      }
    }
    return result;
  }
};

template <typename T, typename Tint> struct IminimumOrbitEntry {
  MyMatrix<Tint> representative;
  T det;
  T covsqr;
  size_t orbit_size;
};

template <typename T, typename Tint> struct IminimumResult {
  int n;
  int i;
  bool certified;
  // The certified value of mu_i(L)^2, or the best proven lower bound
  // when the certification did not succeed.
  T IminimumSqr;
  T MaxDetEnumerated;
  // The determinants of the exceptional directions of the run of
  // BOUND, all disposed of by the enumeration or by direct
  // computation of their value.
  std::vector<T> exceptional_dets;
  std::vector<IminimumOrbitEntry<T, Tint>> l_orbit;
  std::string failure_reason;
};

// The full computation of mu_i(L)^2 with certificate.
template <typename T, typename Tint, typename Tgroup>
IminimumResult<T, Tint>
ComputeIminimumCertified(MyMatrix<T> const &GramMat, int const &i,
                         std::ostream &os) {
  int n = GramMat.rows();
  size_t max_attempt = 20;
  size_t max_bound_failure = 3;
  size_t max_trap = 100000;
  // The starting determinant bound: the minimal determinant of a
  // rank i sublattice, so that the enumeration is nonempty.
  ResultKRankinMinOrbits<T, Tint> res_min =
      Rankin_k_minimum_orbits<T, Tint, Tgroup>(GramMat, i, os);
  T MaxDet = res_min.min;
  std::unordered_map<MyMatrix<Tint>, T> map_covsqr;
  auto get_covsqr = [&](MyMatrix<Tint> const &X) -> T {
    MyMatrix<Tint> Xcan = CanonicalizeSublatticeBasis(X);
    auto iter = map_covsqr.find(Xcan);
    if (iter != map_covsqr.end()) {
      return iter->second;
    }
    MyMatrix<T> RedGram = GetProjectedGramMatrix(GramMat, Xcan);
    T covsqr = ComputeCoveringRadiusSquared<T, Tint, Tgroup>(RedGram, os);
    map_covsqr[Xcan] = covsqr;
    return covsqr;
  };
  auto get_det = [&](MyMatrix<Tint> const &X) -> T {
    MyMatrix<T> X_T = UniversalMatrixConversion<T, Tint>(X);
    return DeterminantMat(MyMatrix<T>(X_T * GramMat * X_T.transpose()));
  };
  auto enumerate_at = [&](T const &M)
      -> std::pair<std::vector<IminimumOrbitEntry<T, Tint>>, T> {
    std::vector<SublatticeOrbit<Tint>> l_orbit =
        Rankin_k_level_orbits<T, Tint, Tgroup>(GramMat, i, M, os);
    std::vector<IminimumOrbitEntry<T, Tint>> l_ent;
    T rmax(0);
    for (auto &eOrbit : l_orbit) {
      T det = get_det(eOrbit.representative);
      T covsqr = get_covsqr(eOrbit.representative);
      if (covsqr > rmax) {
        rmax = covsqr;
      }
      l_ent.push_back(
          {eOrbit.representative, det, covsqr, eOrbit.orbit_size});
    }
    return {std::move(l_ent), std::move(rmax)};
  };
  size_t n_bound_failure = 0;
  T prev_fail_R2(0);
  std::string failure_reason;
  std::vector<IminimumOrbitEntry<T, Tint>> l_entry;
  T R2(0);
  for (size_t attempt = 0; attempt < max_attempt; attempt++) {
#ifdef DEBUG_IMINIMUM
    os << "IMINIMUM: attempt=" << attempt << " MaxDet=" << MaxDet << "\n";
#endif
    std::pair<std::vector<IminimumOrbitEntry<T, Tint>>, T> pairEnum =
        enumerate_at(MaxDet);
    l_entry = std::move(pairEnum.first);
    R2 = pairEnum.second;
    // Stabilization of the candidate value: the determinant bound is
    // doubled until the maximum stops increasing, so that BOUND does
    // not run on a value that a slightly deeper enumeration already
    // improves.
    while (true) {
      T MaxDet2 = T(2) * MaxDet;
      std::pair<std::vector<IminimumOrbitEntry<T, Tint>>, T> pair2 =
          enumerate_at(MaxDet2);
      MaxDet = MaxDet2;
      bool is_stable = (pair2.second == R2);
      l_entry = std::move(pair2.first);
      R2 = pair2.second;
      if (is_stable) {
        break;
      }
    }
#ifdef DEBUG_IMINIMUM
    os << "IMINIMUM: |l_entry|=" << l_entry.size() << " MaxDet=" << MaxDet
       << " R2=" << R2 << "\n";
#endif
    IminimumBound<T, Tint, Tgroup> bnd(GramMat, i, max_trap, os);
    MyMatrix<Tint> S_top = IdentityMat<Tint>(n);
    IminimumBoundResult<Tint> res =
        bnd.bound_rec(S_top, R2, bnd.top_conf, 0);
#ifdef DEBUG_IMINIMUM
    os << "IMINIMUM: BOUND |fails|=" << res.fails.size()
       << " |exceptions|=" << bnd.exceptions.size()
       << " n_trap=" << bnd.n_trap << " aut_skip=" << bnd.n_aut_skip
       << "\n";
#endif
    if (bnd.aborted || res.fails.size() > 0) {
      // The certification failed. Raising the determinant bound can
      // only help through a larger R2, so give up when R2 stalls.
      if (attempt > 0 && R2 == prev_fail_R2) {
        n_bound_failure += 1;
      } else {
        n_bound_failure = 1;
      }
      prev_fail_R2 = R2;
      if (n_bound_failure >= max_bound_failure) {
        if (bnd.aborted) {
          failure_reason = "the number of trap subspaces exceeded " +
                           std::to_string(max_trap);
        } else {
          failure_reason = std::to_string(res.fails.size()) +
                           " unresolved trap constraints";
        }
        break;
      }
      MaxDet = T(2) * MaxDet;
      continue;
    }
    // Disposal of the exceptions. The trap reduction transports the
    // exceptional directions by automorphisms of the full lattice, so
    // the values of the listed representatives settle their whole
    // classes. R2 is the maximum over the full enumeration up to
    // MaxDet, so an exceptional determinant at most MaxDet is
    // automatically covered; an exception beyond it has its exact
    // value computed directly.
    std::set<T> set_exc_det;
    T MaxBadDet(0);
    size_t n_bad = 0;
    for (auto &X : bnd.exceptions) {
      T det = get_det(X);
      set_exc_det.insert(det);
      if (det <= MaxDet) {
        continue;
      }
      T val = get_covsqr(X);
      if (val > R2) {
        n_bad += 1;
        if (det > MaxBadDet) {
          MaxBadDet = det;
        }
      }
    }
    if (n_bad == 0) {
      std::vector<T> exceptional_dets(set_exc_det.begin(), set_exc_det.end());
#ifdef SANITY_CHECK_IMINIMUM
      bool is_attained = false;
      for (auto &entry : l_entry) {
        if (entry.covsqr == R2) {
          is_attained = true;
        }
      }
      if (!is_attained && R2 > 0) {
        std::cerr << "IMINIMUM: the certified value is not attained in the "
                     "enumeration\n";
        throw TerminalException{1};
      }
#endif
      return {n,  i,  true, R2, MaxDet, std::move(exceptional_dets),
              std::move(l_entry), ""};
    }
    // An exceptional direction beyond the enumeration bound has a
    // value above the candidate: raise the bound to reach it and
    // restart with the improved lower bound.
#ifdef DEBUG_IMINIMUM
    os << "IMINIMUM: n_bad=" << n_bad << " MaxBadDet=" << MaxBadDet << "\n";
#endif
    MaxDet = MaxBadDet;
    n_bound_failure = 0;
  }
  if (failure_reason.size() == 0) {
    failure_reason = "the maximal number of attempts was reached";
  }
  return {n,  i,  false, R2, MaxDet, {}, std::move(l_entry),
          std::move(failure_reason)};
}

// clang-format off
#endif  // SRC_DELAUNAY_IMINIMUM_H_
// clang-format on
