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
  lattice) and the recursion is shared through the orbits of the
  sections under the automorphism group of the full lattice: every
  section met by the recursion is registered together with the
  automorphism mapping the orbit representative onto it, the subtree
  is computed once on the representative, and the result is
  transported everywhere else --- the failing constraints map
  through the automorphism, and the exceptional directions need no
  mapping since their disposal only uses their values and
  determinants in the full lattice, which are invariant. On the very
  symmetric lattices this shrinks the computation by orders of
  magnitude.
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

// The certification results are shared through the orbits of the
// sections under Aut(L). If sigma in Aut(L) maps the section S1 onto
// the section S2, it transports the whole certification statement of
// BOUND(S1, T) to BOUND(S2, T): the failing constraints map through
// sigma, and the exceptional directions need no mapping at all,
// since their disposal only uses their values and determinants in
// the full lattice, which are Aut(L)-invariant. A subtree is
// therefore computed once per orbit of sections and replayed, with
// the mapping, everywhere else.
template <typename T, typename Tint, typename Tgroup> struct IminimumBound {
  MyMatrix<T> GramMat;
  int i;
  std::ostream &os;
  // The generators of the automorphism group of the full lattice.
  std::vector<MyMatrix<Tint>> l_autgen;
  // The exceptions, as canonicalized sublattices of the full lattice.
  std::unordered_set<MyMatrix<Tint>> exceptions;
  // The orbit classes of the sections. Every section met by the
  // recursion is registered with its class index and the
  // automorphism mapping the class representative onto it; the
  // certification results are attached to the class, in the
  // coordinates of the representative.
  struct SectionEntry {
    T target;
    IminimumBoundResult<Tint> result;
    std::vector<MyMatrix<Tint>> excs;
  };
  struct SectionClass {
    MyMatrix<Tint> rep;
    std::vector<SectionEntry> entries;
  };
  std::vector<SectionClass> classes;
  std::unordered_map<MyMatrix<Tint>, std::pair<size_t, MyMatrix<Tint>>>
      map_section;
  // Cache of the label losses in dimension at least 3, keyed by the
  // LLL reduced Gram matrix of the dual section: equal keys are
  // isometric label lattices, so the reuse is sound (unequal keys
  // may recompute, which is only a cost).
  std::unordered_map<MyMatrix<T>, T> map_loss;
  // Safety bound on the total number of traps processed.
  size_t n_trap = 0;
  size_t max_trap;
  bool aborted = false;
#ifdef DEBUG_IMINIMUM
  size_t n_compute = 0;
  size_t n_replay = 0;
#endif
  IminimumBound(MyMatrix<T> const &_GramMat, int const &_i,
                size_t const &_max_trap, std::ostream &_os)
      : GramMat(_GramMat), i(_i), os(_os), max_trap(_max_trap) {
    SublatticeStabEquiData<T, Tint> data =
        GetSublatticeStabEquiData<T, Tint>(GramMat, os);
    l_autgen = ArithmeticAutomorphismGroupMultiple_inner<T, Tint, Tgroup>(
        data.ListMat, data.W, os);
  }

  // Registers the orbit of the section S under Aut(L): every member
  // is stored with the automorphism mapping S onto it, S becomes the
  // representative of the new class and its index is returned.
  size_t register_section_orbit(MyMatrix<Tint> const &S) {
    size_t i_class = classes.size();
    classes.push_back({S, {}});
    int n = GramMat.rows();
    std::vector<std::pair<MyMatrix<Tint>, MyMatrix<Tint>>> l_active;
    MyMatrix<Tint> eId = IdentityMat<Tint>(n);
    map_section[S] = {i_class, eId};
    l_active.push_back({S, eId});
    size_t pos = 0;
    // The enumeration of a very large orbit costs more than the
    // sharing it enables, so it is capped: the members beyond the cap
    // are simply registered later as their own classes, which only
    // duplicates work and does not affect the soundness.
    size_t max_orbit = 20000;
    while (pos < l_active.size() && l_active.size() <= max_orbit) {
      std::pair<MyMatrix<Tint>, MyMatrix<Tint>> pairSM = l_active[pos];
      pos++;
      for (auto &eGen : l_autgen) {
        MyMatrix<Tint> Simg = CanonicalizeSublatticeBasis(
            MyMatrix<Tint>(pairSM.first * eGen));
        if (map_section.count(Simg) == 0) {
          MyMatrix<Tint> sigma = pairSM.second * eGen;
          map_section[Simg] = {i_class, sigma};
          l_active.push_back({std::move(Simg), std::move(sigma)});
        }
      }
    }
    return i_class;
  }

  // Transports a result from the coordinates of the representative
  // Srep to the coordinates of the section S = Srep sigma: with the
  // integral matrix U defined by Srep sigma = U S, a failing
  // constraint of dual coordinates y becomes y (U^{-1})^T.
  IminimumBoundResult<Tint> map_result(IminimumBoundResult<Tint> const &res,
                                       MyMatrix<Tint> const &Srep,
                                       MyMatrix<Tint> const &sigma,
                                       MyMatrix<Tint> const &S) {
    if (res.fails.empty()) {
      return res;
    }
    MyMatrix<T> Srep_T = UniversalMatrixConversion<T, Tint>(Srep);
    MyMatrix<T> sigma_T = UniversalMatrixConversion<T, Tint>(sigma);
    MyMatrix<T> S_T = UniversalMatrixConversion<T, Tint>(S);
    MyMatrix<T> U_T = Srep_T * sigma_T * S_T.transpose() *
                      Inverse(MyMatrix<T>(S_T * S_T.transpose()));
#ifdef SANITY_CHECK_IMINIMUM
    MyMatrix<Tint> U = UniversalMatrixConversion<Tint, T>(U_T);
    if (UniversalMatrixConversion<T, Tint>(U) != U_T ||
        MyMatrix<Tint>(U * S) != MyMatrix<Tint>(Srep * sigma)) {
      std::cerr << "map_result: the transport matrix is incorrect\n";
      throw TerminalException{1};
    }
#endif
    MyMatrix<T> V_T = Inverse(U_T).transpose();
    MyMatrix<Tint> Vmap = UniversalMatrixConversion<Tint, T>(V_T);
    MyMatrix<Tint> Vmap_tr = Vmap.transpose();
    IminimumBoundResult<Tint> res_map;
    for (auto &ftup : res.fails) {
      std::vector<MyVector<Tint>> con;
      for (auto &y : ftup.constraints) {
        MyVector<Tint> y2 = Vmap_tr * y;
        con.push_back(y2);
      }
      res_map.fails.push_back({std::move(con)});
    }
    return res_map;
  }

  // Certifies mu(pi_{V^perp})^2 <= target for the rational directions
  // V of dimension i of the section with basis rows S (a saturated
  // sublattice of the full lattice, in canonical form), up to the
  // exceptions collected in this->exceptions (also appended to p_exc
  // when given) and the returned failure tuples. The computation
  // happens on the representative of the orbit class of S and the
  // result is transported.
  IminimumBoundResult<Tint> bound_rec(MyMatrix<Tint> const &S, T const &target,
                                      std::vector<MyMatrix<Tint>> *p_exc) {
    int r = S.rows();
    IminimumBoundResult<Tint> result;
    if (r <= i || aborted) {
      return result;
    }
    if (target <= 0) {
      result.fails.push_back(IminimumFailTuple<Tint>{});
      return result;
    }
    size_t i_class;
    auto iter = map_section.find(S);
    if (iter == map_section.end()) {
      i_class = register_section_orbit(S);
    } else {
      i_class = iter->second.first;
    }
    // The iterator can be invalidated by the registrations of the
    // recursion, so the mapping is copied out.
    MyMatrix<Tint> sigma = map_section.at(S).second;
    // A result for a smaller target is a stronger statement and can
    // be replayed conservatively; an entry with failures is only
    // replayed on the exact target, so that its failing constraints
    // do not spam the caller with unneeded joins.
    auto find_best = [&]() -> SectionEntry * {
      SectionEntry *best = nullptr;
      for (auto &cc : classes[i_class].entries) {
        bool usable = cc.result.fails.empty() ? (cc.target <= target)
                                              : (cc.target == target);
        if (usable) {
          if (!best || cc.target > best->target) {
            best = &cc;
          }
        }
      }
      return best;
    };
    SectionEntry *best = find_best();
    if (best == nullptr) {
      // The recursion registers new classes, which can reallocate the
      // vector, so the representative is copied out.
      MyMatrix<Tint> rep = classes[i_class].rep;
      SectionEntry eEnt = compute_entry(rep, target);
      if (aborted) {
        // A partial computation is not cached; the driver discards
        // the whole run anyway.
        return eEnt.result;
      }
      classes[i_class].entries.push_back(std::move(eEnt));
      best = &classes[i_class].entries.back();
    }
#ifdef DEBUG_IMINIMUM
    else {
      n_replay += 1;
    }
#endif
    for (auto &eX : best->excs) {
      exceptions.insert(eX);
      if (p_exc) {
        p_exc->push_back(eX);
      }
    }
    if (classes[i_class].rep == S) {
      return best->result;
    }
    return map_result(best->result, classes[i_class].rep, sigma, S);
  }

  // The actual computation of BOUND on the class representative S.
  SectionEntry compute_entry(MyMatrix<Tint> const &S, T const &target) {
#ifdef DEBUG_IMINIMUM
    n_compute += 1;
#endif
    int r = S.rows();
    IminimumBoundResult<Tint> result;
    std::vector<MyMatrix<Tint>> local_excs;
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
          local_excs.push_back(Xfull);
        }
      }
      return {target, std::move(result), std::move(local_excs)};
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
    os << "IMINIMUM: compute_entry r=" << r << " target=" << target
       << " |short duals|=" << work.size() << "\n";
#endif
    std::unordered_set<MyMatrix<Tint>> seen;
    while (!work.empty()) {
      std::vector<MyVector<Tint>> tup = work.back();
      work.pop_back();
      MyMatrix<Tint> Y = MatrixFromVectorFamily(tup);
      MyMatrix<Tint> Ys = CanonicalizeSublatticeBasis(
          SaturationSpanOfRows(Y));
      if (seen.count(Ys) == 1) {
        continue;
      }
      seen.insert(Ys);
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
      T loss;
      if (kW <= 2) {
        loss = LabelLatticeCovSqr<T, Tint, Tgroup>(GramDualSection, os);
      } else {
        LLLreduction<T, Tint> recLLL =
            LLLreducedBasis<T, Tint>(GramDualSection, os);
        auto iterL = map_loss.find(recLLL.GramMatRed);
        if (iterL != map_loss.end()) {
          loss = iterL->second;
        } else {
          loss = LabelLatticeCovSqr<T, Tint, Tgroup>(GramDualSection, os);
          map_loss[recLLL.GramMatRed] = loss;
        }
      }
      T sub_target = target - loss;
#ifdef DEBUG_IMINIMUM
      os << "IMINIMUM: trap dim=" << kW << " loss=" << loss
         << " sub_target=" << sub_target << "\n";
#endif
      if (sub_target <= 0) {
        result.fails.push_back({tup});
        continue;
      }
      // The sub target is rounded down to a fixed grid: certifying a
      // smaller target is a stronger statement, so the rounding is
      // conservative, and it collapses the near-identical targets of
      // the equivalent sections onto common cache entries.
      T quant(1024);
      Tint num_q = UniversalFloorScalarInteger<Tint, T>(sub_target * quant);
      sub_target = UniversalScalarConversion<T, Tint>(num_q) / quant;
      if (sub_target <= 0) {
        result.fails.push_back({tup});
        continue;
      }
      MyMatrix<Tint> K = CanonicalizeSublatticeBasis(NullspaceIntTrMat(Ys));
      MyMatrix<Tint> Ssub = CanonicalizeSublatticeBasis(
          MyMatrix<Tint>(K * S));
      IminimumBoundResult<Tint> sub =
          bound_rec(Ssub, sub_target, &local_excs);
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
    return {target, std::move(result), std::move(local_excs)};
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
  size_t max_trap = 1000000;
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
    IminimumBoundResult<Tint> res = bnd.bound_rec(S_top, R2, nullptr);
#ifdef DEBUG_IMINIMUM
    os << "IMINIMUM: BOUND |fails|=" << res.fails.size()
       << " |exceptions|=" << bnd.exceptions.size()
       << " n_trap=" << bnd.n_trap << " n_compute=" << bnd.n_compute
       << " n_replay=" << bnd.n_replay << "\n";
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
