// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_IGUSA_IGUSA_FACE_ENUM_H_
#define SRC_IGUSA_IGUSA_FACE_ENUM_H_

// clang-format off
#include "igusa_initial_facet.h"
#include <map>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_IGUSA_FACE_ENUM
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_IGUSA_FACE_ENUM
#endif

#ifdef TIMINGS
#define TIMINGS_IGUSA_FACE_ENUM
#endif

/*
  The enumeration of the facets of P = conv(I_n) by the recursive adjacency
  decomposition on the faces (paper, sections "The structure of the faces"
  to "Connectivity").

  A face is stored with a valid pair (F, b), valid on all of I_n, whose zero
  set on P is exactly the face. From it: the rank r, a basis B of the
  lattice W, the vertices of the bounded face Phi0 of P_r (in the basis B,
  from the representatives of the orbits of vertices of P_r), their
  barycenter sigma and the canonical functional M = B^T sigma^{-1} B.

  * Equivalence: the faces Psi contained in a chain of faces Phi_1, ... are
    equivalent under the common stabilizer iff the lists
    (sigma(Psi)^{-1}, M_{Phi_1} restricted to W(Psi), ...) are equivalent
    under GL_r(Z). Canonical forms give a complete key.
  * Minimization: only positive semidefinite functionals over I_n, reduced
    to the rank of their image, where they are positive definite.
  * A functional h valid on a face Phi is tested, and extended to a pair
    valid on I_n, by minimizing h + mu f_Phi over I_n for increasing mu.
  * The flip of Theta subset Psi subset Phi: the smallest s such that
    g + s f >= 0 on Phi, g the functional of Theta, f the one of Psi.
  * The facets of a face Phi of rank r: the vertical ones from the polytope
    Phi0, the non vertical ones by the walk of Proposition
    PROP_nonvertical_connected, starting from Psi_w or from a given seed.
 */

// tr(F X) >= b
template <typename T> struct IgusaPair {
  MyMatrix<T> F;
  T b;
};

template <typename T> T igusa_pair_eval(IgusaPair<T> const &p, MyMatrix<T> const &X) {
  return frobenius_inner(p.F, X) - p.b;
}

template <typename T> struct IgusaFace {
  IgusaPair<T> pair;
  int rank;
  // r x n basis of W (rows)
  MyMatrix<T> B;
  MyMatrix<T> F0;
  // vertices of Phi0 in the basis B, and their type (representative index)
  std::vector<MyMatrix<T>> l_vert;
  std::vector<size_t> l_type;
  MyMatrix<T> sigma;
  // B^T sigma^{-1} B, zero for P itself
  MyMatrix<T> M;
  // dimension of Phi0
  int dim0;
};

template <typename T> int igusa_face_dim(IgusaFace<T> const &face, int n) {
  int r = face.rank;
  return face.dim0 + n * (n + 1) / 2 - r * (r + 1) / 2;
}

template <typename T, typename Tint, typename Tgroup> struct IgusaFaceContext {
  IgusaParameters<T> params;
  int n;
  // the vertex representatives of P_r by rank
  std::map<int, std::vector<MyMatrix<T>>> map_repr;
  std::map<int, IgusaSpace<T, Tint>> map_space;
  std::ostream &os;
  size_t n_face_built = 0;
  size_t n_flip = 0;
  size_t n_minimization = 0;

  IgusaSpace<T, Tint> &get_space(int r) {
    auto iter = map_space.find(r);
    if (iter != map_space.end()) {
      return iter->second;
    }
    LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(r);
    IgusaSpace<T, Tint> space = build_igusa_space<T, Tint>(LinSpa, params, os);
    igusa_insert_initial_cuts(space);
    map_space.emplace(r, std::move(space));
    return map_space.at(r);
  }

  std::vector<MyMatrix<T>> const &get_repr(int r) {
    auto iter = map_repr.find(r);
    if (iter == map_repr.end()) {
      std::cerr << "IGUSA_FACE_ENUM: no vertex representatives of rank " << r
                << "\n";
      throw TerminalException{1};
    }
    return iter->second;
  }
};

// The face of the whole polyhedron P (rank 0)
template <typename T> IgusaFace<T> igusa_whole_face(int n) {
  IgusaFace<T> face;
  face.pair = {ZeroMatrix<T>(n, n), T(0)};
  face.rank = 0;
  face.B = MyMatrix<T>(0, n);
  face.F0 = MyMatrix<T>(0, 0);
  face.sigma = MyMatrix<T>(0, 0);
  face.M = ZeroMatrix<T>(n, n);
  face.dim0 = 0;
  return face;
}

template <typename T> IgusaPair<T> igusa_normalize_pair(IgusaPair<T> const &p) {
  std::pair<MyMatrix<T>, T> pp = igusa_normalize_facet(p.F, p.b);
  return {pp.first, pp.second};
}

template <typename T, typename Tint, typename Tgroup>
IgusaFace<T> igusa_face_from_pair(IgusaFaceContext<T, Tint, Tgroup> &ctx,
                                  IgusaPair<T> const &pair_in) {
  ctx.n_face_built++;
  IgusaPair<T> pair = igusa_normalize_pair(pair_in);
  IgusaImageFactorization<T> fact = igusa_image_factorization(pair.F);
  IgusaFace<T> face;
  face.pair = pair;
  face.B = fact.B;
  face.F0 = fact.F0;
  face.rank = fact.B.rows();
  IgusaSpace<T, Tint> &space_r = ctx.get_space(face.rank);
  IgusaFaceVertices<T> info;
  std::vector<std::pair<MyMatrix<T>, size_t>> l_pair =
      igusa_bounded_face_vertices<T, Tint, Tgroup>(
          space_r, face.F0, pair.b, ctx.get_repr(face.rank), info, ctx.os);
  if (l_pair.empty()) {
    std::cerr << "IGUSA_FACE_ENUM: a face without vertices, the pair is not "
                 "valid or the representatives are incomplete\n";
    throw TerminalException{1};
  }
  int r = face.rank;
  face.sigma = ZeroMatrix<T>(r, r);
  for (auto &pv : l_pair) {
    face.l_vert.push_back(pv.first);
    face.l_type.push_back(pv.second);
    face.sigma += pv.first;
  }
  face.sigma /= T(static_cast<int>(l_pair.size()));
  face.M = face.B.transpose() * Inverse(face.sigma) * face.B;
  face.dim0 = igusa_affine_dimension(face.l_vert);
  return face;
}

// X -> g X g^T applied to the face
template <typename T>
IgusaFace<T> igusa_transform_face(IgusaFace<T> const &face,
                                  MyMatrix<T> const &g) {
  MyMatrix<T> ginv = Inverse(g);
  IgusaFace<T> ret = face;
  ret.pair.F = ginv.transpose() * face.pair.F * ginv;
  ret.B = face.B * ginv;
  ret.M = ginv.transpose() * face.M * ginv;
  return ret;
}

// F' with M = B^T F' B, for M vanishing on the kernel of B
template <typename T>
MyMatrix<T> igusa_restrict_functional(MyMatrix<T> const &M,
                                      MyMatrix<T> const &B) {
  MyMatrix<T> BBt_inv = Inverse(MyMatrix<T>(B * B.transpose()));
  MyMatrix<T> Fr = BBt_inv * B * M * B.transpose() * BBt_inv;
#ifdef SANITY_CHECK_IGUSA_FACE_ENUM
  if (B.transpose() * Fr * B != M) {
    std::cerr << "IGUSA_FACE_ENUM: the functional does not factor through B\n";
    throw TerminalException{1};
  }
#endif
  return Fr;
}

// A unimodular matrix whose first rows are B
template <typename T> MyMatrix<T> igusa_completion(MyMatrix<T> const &B, int n) {
  if (B.rows() == n) {
    return B;
  }
  return Concatenate(B, SubspaceCompletionInt(B, n));
}

/*
  The g in GL_n(Z) with B2 g = h B1, h acting on the forms of W in the bases
  B1, B2: it maps the lift through W1 of a face Phi0 onto the lift through W2
  of h Phi0 h^T.
 */
template <typename T>
MyMatrix<T> igusa_lift_transformation(MyMatrix<T> const &B1,
                                      MyMatrix<T> const &B2,
                                      MyMatrix<T> const &h, int n) {
  int r = h.rows();
  MyMatrix<T> U1 = igusa_completion(B1, n);
  MyMatrix<T> U2 = igusa_completion(B2, n);
  MyMatrix<T> big = IdentityMat<T>(n);
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < r; j++) {
      big(i, j) = h(i, j);
    }
  }
  return Inverse(U2) * big * U1;
}

template <typename T> struct IgusaFaceKey {
  std::string key;
  // Bc with Bc L Bc^T canonical for the list L
  MyMatrix<T> Bc;
};

/*
  The key of a face Psi under the common stabilizer of the faces of
  l_super (given by their canonical functionals M), all containing Psi.
 */
template <typename T, typename Tint, typename Tgroup>
IgusaFaceKey<T> igusa_face_key(IgusaFace<T> const &face,
                               std::vector<MyMatrix<T>> const &l_super_M,
                               std::ostream &os) {
  std::vector<MyMatrix<T>> ListM{RemoveFractionMatrix(Inverse(face.sigma))};
  for (auto &M : l_super_M) {
    MyMatrix<T> Mr = igusa_restrict_functional(M, face.B);
    if (IsZeroMatrix(Mr)) {
      ListM.push_back(Mr);
    } else {
      ListM.push_back(RemoveFractionMatrix(Mr));
    }
  }
  MyMatrix<Tint> Bc = ComputeCanonicalFormMultiple<T, Tint, Tgroup>(ListM, os);
  MyMatrix<T> Bc_T = UniversalMatrixConversion<T, Tint>(Bc);
  std::ostringstream os_key;
  os_key << face.rank << " " << face.l_vert.size();
  for (auto &L : ListM) {
    os_key << " " << StringMatrixGAP(MyMatrix<T>(Bc_T * L * Bc_T.transpose()));
  }
  return {os_key.str(), Bc_T};
}

/*
  For faces with the same key under GL_n(Z): the g mapping face1 onto face2.
 */
template <typename T>
MyMatrix<T> igusa_face_equivalence(IgusaFace<T> const &face1,
                                   IgusaFaceKey<T> const &key1,
                                   IgusaFace<T> const &face2,
                                   IgusaFaceKey<T> const &key2, int n) {
  // k sigma1^{-1} k^T = sigma2^{-1} up to the common scaling, h = k^{-T}
  MyMatrix<T> k = Inverse(key2.Bc) * key1.Bc;
  MyMatrix<T> h = Inverse(k).transpose();
  MyMatrix<T> g = igusa_lift_transformation(face1.B, face2.B, h, n);
#ifdef SANITY_CHECK_IGUSA_FACE_ENUM
  if (igusa_transform_face(face1, g).M != face2.M) {
    std::cerr << "IGUSA_FACE_ENUM: the equivalence does not map the faces\n";
    throw TerminalException{1};
  }
#endif
  return g;
}

/*
  The minimum over I_n of tr(H X), H positive semidefinite, and a
  minimizer. The functional is reduced to the rank s of its image, where it
  is positive definite, and the minimizer of I_s is extended to I_n.
 */
template <typename T, typename Tint, typename Tgroup>
std::pair<T, MyMatrix<T>> igusa_min_psd(IgusaFaceContext<T, Tint, Tgroup> &ctx,
                                        MyMatrix<T> const &H) {
  ctx.n_minimization++;
  int n = ctx.n;
  if (IsZeroMatrix(H)) {
    return {T(0), IdentityMat<T>(n)};
  }
  // A positive scaling to a primitive integral matrix keeps the integer
  // programs small; the value is scaled back.
  FractionMatrix<T> fr = RemoveFractionMatrixPlusCoeff(H);
  IgusaImageFactorization<T> fact = igusa_image_factorization(fr.TheMat);
  int s = fact.B.rows();
  IgusaSpace<T, Tint> &space_s = ctx.get_space(s);
  MyVector<T> f = igusa_trace_function(space_s, fact.F0);
  igusa_insert_seed_cuts(space_s, fact.F0, ctx.os);
  MyMatrix<T> NoIneq(0, space_s.dim + 1);
  MyMatrix<T> NoEqua(0, space_s.dim + 1);
  std::optional<IgusaMinimum<T>> opt =
      igusa_integral_minimization(space_s, f, NoIneq, NoEqua, ctx.os);
  IgusaMinimum<T> res = unfold_opt(opt, "The minimum over I_s");
  MyMatrix<T> Y = igusa_matrix(space_s, res.x);
  MyMatrix<T> X = igusa_extend_form(fact.B, Y);
  T value = res.value / fr.TheMult;
#ifdef SANITY_CHECK_IGUSA_FACE_ENUM
  if (frobenius_inner(H, X) != value) {
    std::cerr << "IGUSA_FACE_ENUM: the extended minimizer has another value\n";
    throw TerminalException{1};
  }
#endif
  return {value, X};
}

// A point of the face (for P, the identity)
template <typename T>
MyMatrix<T> igusa_face_point(IgusaFace<T> const &face, int n) {
  if (face.rank == 0) {
    return IdentityMat<T>(n);
  }
  return igusa_extend_form(face.B, face.l_vert[0]);
}

// A basis (rows) of the integral kernel of F, the identity for F = 0
template <typename T> MyMatrix<T> igusa_kernel_basis(MyMatrix<T> const &F) {
  int n = F.rows();
  if (IsZeroMatrix(F)) {
    return IdentityMat<T>(n);
  }
  return NullspaceIntMat(F);
}

template <typename T> struct IgusaValidity {
  bool valid;
  // a point of the face where h < 0
  MyMatrix<T> X;
  // h + mu f_Phi >= 0 on I_n
  T mu;
};

/*
  Whether h(X) = tr(M_h X) - c_h is nonnegative on the face Phi. If yes, a
  mu with h + mu f_Phi >= 0 on I_n; if no, a point of Phi where h < 0.
 */
template <typename T, typename Tint, typename Tgroup>
IgusaValidity<T> igusa_face_validity(IgusaFaceContext<T, Tint, Tgroup> &ctx,
                                     IgusaFace<T> const &phi,
                                     MyMatrix<T> const &M_h, T const &c_h) {
  int n = ctx.n;
  std::ostream &os = ctx.os;
  MyMatrix<T> const &F_phi = phi.pair.F;
  T const &b_phi = phi.pair.b;
  auto h_eval = [&](MyMatrix<T> const &X) -> T {
    return frobenius_inner(M_h, X) - c_h;
  };
  // h decreasing along a recession direction of Phi
  MyMatrix<T> K = igusa_kernel_basis(F_phi);
  MyMatrix<T> QK = K * M_h * K.transpose();
  if (K.rows() > 0 && !IsPositiveSemiDefinite(QK, os)) {
    MyVector<Tint> vK = igusa_non_positive_vector<T, Tint>(QK, T(0), os);
    MyVector<T> v = K.transpose() * UniversalVectorConversion<T, Tint>(vK);
    MyMatrix<T> vvT = v * v.transpose();
    MyMatrix<T> X = igusa_face_point(phi, n);
    T t(1);
    while (h_eval(X + t * vvT) >= 0) {
      t *= 2;
    }
    return {false, MyMatrix<T>(X + t * vvT), T(0)};
  }
  T mu(0);
  if (!IsZeroMatrix(F_phi)) {
    mu = T(1);
    T mu_max(1);
    for (int i = 0; i < 60; i++) {
      mu_max *= 2;
    }
    while (!IsPositiveSemiDefinite(MyMatrix<T>(M_h + mu * F_phi), os)) {
      mu *= 2;
      if (mu > mu_max) {
        std::cerr << "IGUSA_FACE_ENUM: the functional is positive "
                     "semidefinite on the kernel of the face but no multiple "
                     "of the face makes it positive semidefinite\n";
        throw TerminalException{1};
      }
    }
    // twice a valid mu: the kernel of the sum is in the kernel of F_phi
    mu *= 2;
  }
  while (true) {
    MyMatrix<T> H = M_h + mu * F_phi;
    std::pair<T, MyMatrix<T>> pmin = igusa_min_psd(ctx, H);
    T val = pmin.first - (c_h + mu * b_phi);
    if (val >= 0) {
      return {true, MyMatrix<T>(), mu};
    }
    MyMatrix<T> const &X = pmin.second;
    T fX = frobenius_inner(F_phi, X) - b_phi;
    if (fX == 0) {
      return {false, X, mu};
    }
    T mu_new = (c_h - frobenius_inner(M_h, X)) / fX;
    if (mu_new <= mu) {
      std::cerr << "IGUSA_FACE_ENUM: the multiplier does not increase\n";
      throw TerminalException{1};
    }
    mu = mu_new;
  }
}

/*
  The pair of the face Phi cap {h = 0}, for h >= 0 on Phi: h + (mu+1) f_Phi,
  whose zero set on P is exactly that face.
 */
template <typename T, typename Tint, typename Tgroup>
IgusaPair<T> igusa_extend_to_pair(IgusaFaceContext<T, Tint, Tgroup> &ctx,
                                  IgusaFace<T> const &phi,
                                  MyMatrix<T> const &M_h, T const &c_h) {
  IgusaValidity<T> v = igusa_face_validity(ctx, phi, M_h, c_h);
  if (!v.valid) {
    std::cerr << "IGUSA_FACE_ENUM: the functional is not valid on the face\n";
    throw TerminalException{1};
  }
  T mu = v.mu + 1;
  return {MyMatrix<T>(M_h + mu * phi.pair.F), c_h + mu * phi.pair.b};
}

// The facets of the polytope Phi0, as faces of P (vertical facets of Phi).
template <typename T, typename Tint, typename Tgroup>
std::vector<IgusaFace<T>>
igusa_polytope_facets(IgusaFaceContext<T, Tint, Tgroup> &ctx,
                      IgusaFace<T> const &phi) {
  std::vector<IgusaFace<T>> l_ret;
  if (phi.rank == 0 || phi.dim0 == 0) {
    return l_ret;
  }
  int r = phi.rank;
  IgusaSpace<T, Tint> &space_r = ctx.get_space(r);
  std::vector<MyVector<T>> l_x;
  for (auto &Y : phi.l_vert) {
    l_x.push_back(igusa_coordinates(space_r, Y));
  }
  std::vector<MyMatrix<T>> ListM{RemoveFractionMatrix(Inverse(phi.sigma))};
  std::vector<MyMatrix<Tint>> l_g =
      igusa_functional_stabilizer<T, Tint, Tgroup>(ListM, ctx.os);
  Tgroup GRP;
  std::vector<IgusaRidge<T>> l_ridge =
      igusa_facet_ridges<T, Tint, Tgroup>(space_r, l_x, l_g, GRP, ctx.os);
  for (auto &ridge : l_ridge) {
    MyMatrix<T> Gr = igusa_functional_matrix(space_r, ridge.g);
    MyMatrix<T> M_h = phi.B.transpose() * Gr * phi.B;
    T c_h = -ridge.g(0);
    IgusaPair<T> pair = igusa_extend_to_pair(ctx, phi, M_h, c_h);
    IgusaFace<T> face = igusa_face_from_pair(ctx, pair);
#ifdef SANITY_CHECK_IGUSA_FACE_ENUM
    if (face.l_vert.size() != ridge.face.count() || face.rank != r) {
      std::cerr << "IGUSA_FACE_ENUM: the lifted facet of the polytope has "
                << face.l_vert.size() << " vertices and rank " << face.rank
                << " instead of " << ridge.face.count() << " and " << r
                << "\n";
      throw TerminalException{1};
    }
#endif
    l_ret.push_back(face);
  }
  return l_ret;
}

/*
  The smallest s with Qg + s Qf positive semidefinite, Qf positive
  semidefinite: an interval [lo, hi] of length at most 2^-40 containing it,
  or none if it is -infinity.
 */
template <typename T>
std::optional<std::pair<T, T>> igusa_psd_threshold(MyMatrix<T> const &Qg,
                                                   MyMatrix<T> const &Qf,
                                                   std::ostream &os) {
  if (Qg.rows() == 0) {
    return {};
  }
  auto is_psd = [&](T const &s) -> bool {
    return IsPositiveSemiDefinite(MyMatrix<T>(Qg + s * Qf), os);
  };
  T big(1);
  for (int i = 0; i < 60; i++) {
    big *= 2;
  }
  T hi(1);
  while (!is_psd(hi)) {
    hi *= 2;
    if (hi > big) {
      std::cerr << "IGUSA_FACE_ENUM: the pencil is never positive "
                   "semidefinite\n";
      throw TerminalException{1};
    }
  }
  T lo(-1);
  while (is_psd(lo)) {
    lo *= 2;
    if (-lo > big) {
      return {};
    }
  }
  T eps(1);
  for (int i = 0; i < 40; i++) {
    eps /= 2;
  }
  while (hi - lo > eps) {
    T mid = (lo + hi) / 2;
    if (is_psd(mid)) {
      hi = mid;
    } else {
      lo = mid;
    }
  }
  return std::make_pair(lo, hi);
}

/*
  The flip: for Theta a facet of Psi and Psi a facet of Phi, the other facet
  Psi' of Phi containing Theta. With g, f the functionals of Theta and Psi,
  it is the zero set on Phi of g + s f for the smallest s with g + s f >= 0
  on Phi.
 */
template <typename T, typename Tint, typename Tgroup>
IgusaFace<T> igusa_face_flip(IgusaFaceContext<T, Tint, Tgroup> &ctx,
                             IgusaFace<T> const &phi, IgusaFace<T> const &psi,
                             IgusaFace<T> const &theta) {
  ctx.n_flip++;
  int n = ctx.n;
  std::ostream &os = ctx.os;
  MyMatrix<T> const &G = theta.pair.F;
  T const &c_g = theta.pair.b;
  MyMatrix<T> const &Ff = psi.pair.F;
  T const &c_f = psi.pair.b;
  auto g_eval = [&](MyMatrix<T> const &X) -> T {
    return frobenius_inner(G, X) - c_g;
  };
  auto f_eval = [&](MyMatrix<T> const &X) -> T {
    return frobenius_inner(Ff, X) - c_f;
  };
  // Points of Phi off Psi: the lifted vertices of Phi0 and their translates
  // by v v^T, v in the kernel of Phi.
  MyMatrix<T> K = igusa_kernel_basis(phi.pair.F);
  std::vector<MyMatrix<T>> l_base;
  if (phi.rank == 0) {
    l_base.push_back(IdentityMat<T>(n));
  } else {
    for (auto &Y : phi.l_vert) {
      l_base.push_back(igusa_extend_form(phi.B, Y));
    }
  }
  std::optional<T> s_c;
  auto consider = [&](MyMatrix<T> const &X) -> void {
    T fX = f_eval(X);
    if (fX > 0) {
      T ratio = -g_eval(X) / fX;
      if (!s_c || ratio > *s_c) {
        s_c = ratio;
      }
    }
  };
  for (auto &X : l_base) {
    consider(X);
    for (int i = 0; i < K.rows(); i++) {
      MyVector<T> v = GetMatrixRow(K, i);
      consider(MyMatrix<T>(X + v * v.transpose()));
    }
  }
  if (!s_c) {
    std::cerr << "IGUSA_FACE_ENUM: no point of the face off the facet\n";
    throw TerminalException{1};
  }
  MyMatrix<T> Qg = K * G * K.transpose();
  MyMatrix<T> Qf = K * Ff * K.transpose();
  std::optional<std::pair<T, T>> thr = igusa_psd_threshold(Qg, Qf, os);
  auto test = [&](T const &s) -> IgusaValidity<T> {
    MyMatrix<T> M_h = G + s * Ff;
    T c_h = c_g + s * c_f;
    return igusa_face_validity(ctx, phi, M_h, c_h);
  };
  auto finish = [&](T const &s) -> IgusaFace<T> {
    MyMatrix<T> M_h = G + s * Ff;
    T c_h = c_g + s * c_f;
    IgusaPair<T> pair = igusa_extend_to_pair(ctx, phi, M_h, c_h);
    IgusaFace<T> face = igusa_face_from_pair(ctx, pair);
#ifdef SANITY_CHECK_IGUSA_FACE_ENUM
    MyMatrix<T> Xt = igusa_face_point(theta, n);
    if (igusa_pair_eval(face.pair, Xt) != 0) {
      std::cerr << "IGUSA_FACE_ENUM: the flipped face does not contain the "
                   "codimension 2 face\n";
      throw TerminalException{1};
    }
    if (face.M == psi.M) {
      std::cerr << "IGUSA_FACE_ENUM: the flip returned the same facet\n";
      throw TerminalException{1};
    }
#endif
    return face;
  };
  // Lower bound s_c from the points: if above the threshold, iterate on
  // the violating points.
  T s = *s_c;
  bool use_points = !thr || s > thr->second;
  if (use_points) {
    while (true) {
      IgusaValidity<T> v = test(s);
      if (v.valid) {
        return finish(s);
      }
      s = -g_eval(v.X) / f_eval(v.X);
    }
  }
  // Just above the threshold, where the functionals are coercive on the
  // kernel, and closer if no point violates there.
  T delta(1);
  T delta_min(1);
  for (int i = 0; i < 30; i++) {
    delta_min /= 2;
  }
  T hi = thr->second;
  s = hi + delta;
  bool updated = false;
  while (true) {
    IgusaValidity<T> v = test(s);
    if (!v.valid) {
      s = -g_eval(v.X) / f_eval(v.X);
      updated = true;
      continue;
    }
    if (updated) {
      return finish(s);
    }
    delta /= 2;
    if (delta < delta_min) {
      break;
    }
    s = hi + delta;
  }
  // At the threshold, with a singular functional on the kernel
  T q = igusa_simplest_rational(thr->first, thr->second);
#ifdef DEBUG_IGUSA_FACE_ENUM
  os << "IGUSA_FACE_ENUM: flip at the threshold s=" << q << "\n";
#endif
  IgusaValidity<T> vq = test(q);
  if (!vq.valid) {
    std::cerr << "IGUSA_FACE_ENUM: degenerate flip, the threshold is not "
                 "valid\n";
    throw TerminalException{1};
  }
  return finish(q);
}

/*
  The facets of the face Psi (from Enum(Psi), up to Stab(Psi)) split into
  orbits under Stab(Psi) cap Stab(Phi): the finite part G of Stab(Psi) acts
  on the restriction N of M_Phi to W(Psi), its stabilizer is the finite
  part H of the intersection, and the classes are given by g^{-1} Theta for
  g with g N running over the orbit of N.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<IgusaFace<T>>
igusa_orbit_splitting(IgusaFaceContext<T, Tint, Tgroup> &ctx,
                      IgusaFace<T> const &phi, IgusaFace<T> const &psi,
                      std::vector<IgusaFace<T>> const &l_theta) {
  if (phi.rank == 0) {
    return l_theta;
  }
  int n = ctx.n;
  int r = psi.rank;
  std::vector<MyMatrix<T>> ListM{RemoveFractionMatrix(Inverse(psi.sigma))};
  std::vector<MyMatrix<Tint>> l_gen =
      igusa_functional_stabilizer<T, Tint, Tgroup>(ListM, ctx.os);
  MyMatrix<T> N = igusa_restrict_functional(phi.M, psi.B);
  std::map<MyMatrix<T>, MyMatrix<T>> map_orbit;
  std::vector<MyMatrix<T>> l_pt{N};
  map_orbit[N] = IdentityMat<T>(r);
  size_t pos = 0;
  while (pos < l_pt.size()) {
    MyMatrix<T> p = l_pt[pos];
    MyMatrix<T> e = map_orbit.at(p);
    pos++;
    for (auto &gen : l_gen) {
      MyMatrix<T> g = UniversalMatrixConversion<T, Tint>(gen);
      MyMatrix<T> ginv = Inverse(g);
      MyMatrix<T> q = ginv.transpose() * p * ginv;
      if (map_orbit.count(q) == 0) {
        map_orbit[q] = g * e;
        l_pt.push_back(q);
      }
    }
  }
  if (l_pt.size() == 1) {
    return l_theta;
  }
  std::vector<MyMatrix<T>> l_super{psi.M, phi.M};
  std::vector<IgusaFace<T>> l_ret;
  for (auto &theta : l_theta) {
    std::map<std::string, size_t> map_key;
    for (auto &p : l_pt) {
      MyMatrix<T> k = Inverse(map_orbit.at(p));
      MyMatrix<T> gk = igusa_lift_transformation(psi.B, psi.B, k, n);
      IgusaFace<T> th = igusa_transform_face(theta, gk);
      std::string key =
          igusa_face_key<T, Tint, Tgroup>(th, l_super, ctx.os).key;
      if (map_key.count(key) == 0) {
        map_key[key] = l_ret.size();
        l_ret.push_back(th);
      }
    }
  }
#ifdef DEBUG_IGUSA_FACE_ENUM
  ctx.os << "IGUSA_FACE_ENUM: orbit splitting [G:H]=" << l_pt.size() << " "
         << l_theta.size() << " -> " << l_ret.size() << "\n";
#endif
  return l_ret;
}

template <typename T> struct IgusaBankEntry {
  IgusaFace<T> face;
  IgusaFaceKey<T> key;
  std::vector<IgusaFace<T>> l_facet;
};

template <typename T> struct IgusaFaceBank {
  std::vector<IgusaBankEntry<T>> l_entry;
  std::map<std::string, size_t> map_key;
};

template <typename T> std::string igusa_face_summary(IgusaFace<T> const &f, int n) {
  std::ostringstream os;
  os << "rank=" << f.rank << " |vert|=" << f.l_vert.size()
     << " dim=" << igusa_face_dim(f, n);
  return os.str();
}

/*
  The facets of the face phi up to Stab(phi). For phi = P (rank 0) use
  igusa_enumerate_facets instead, which seeds the lower rank facets.
 */
template <typename T, typename Tint, typename Tgroup>
std::vector<IgusaFace<T>>
igusa_face_facets(IgusaFaceContext<T, Tint, Tgroup> &ctx, IgusaFaceBank<T> &bank,
                  IgusaFace<T> const &phi,
                  std::optional<IgusaFace<T>> const &seed, int depth) {
  int n = ctx.n;
  std::ostream &os = ctx.os;
  IgusaFaceKey<T> key_phi = igusa_face_key<T, Tint, Tgroup>(phi, {}, os);
  auto iter = bank.map_key.find(key_phi.key);
  if (iter != bank.map_key.end()) {
    IgusaBankEntry<T> const &entry = bank.l_entry[iter->second];
    MyMatrix<T> g =
        igusa_face_equivalence(entry.face, entry.key, phi, key_phi, n);
    std::vector<IgusaFace<T>> l_ret;
    for (auto &facet : entry.l_facet) {
      l_ret.push_back(igusa_transform_face(facet, g));
    }
    return l_ret;
  }
#ifdef TIMINGS_IGUSA_FACE_ENUM
  MicrosecondTime time;
#endif
#ifdef DEBUG_IGUSA_FACE_ENUM
  os << "IGUSA_FACE_ENUM: depth=" << depth << " enter "
     << igusa_face_summary(phi, n) << "\n";
#endif
  std::vector<IgusaFace<T>> l_facet = igusa_polytope_facets(ctx, phi);
  if (phi.rank < n) {
    // The walk on the non vertical facets
    std::vector<IgusaFace<T>> l_nv;
    std::vector<std::optional<IgusaFace<T>>> l_seed;
    std::vector<uint8_t> l_done;
    std::map<std::string, size_t> map_nv;
    std::vector<MyMatrix<T>> l_super{phi.M};
    auto insert = [&](IgusaFace<T> const &psi,
                      std::optional<IgusaFace<T>> const &psi_seed) -> void {
      std::string key = igusa_face_key<T, Tint, Tgroup>(psi, l_super, os).key;
      if (map_nv.count(key) > 0) {
        return;
      }
      map_nv[key] = l_nv.size();
      l_nv.push_back(psi);
      l_seed.push_back(psi_seed);
      l_done.push_back(0);
    };
    IgusaInitialFacet<T, Tint> init = igusa_initial_facet<T, Tint, Tgroup>(
        ctx.params, phi.pair.F, phi.pair.b, ctx.get_repr(phi.rank), os);
    if (init.codim == 1) {
      insert(igusa_face_from_pair(ctx, IgusaPair<T>{init.Fnew, init.bnew}),
             {});
    } else if (seed && seed->rank > phi.rank) {
      insert(*seed, {});
    } else {
      std::cerr << "IGUSA_FACE_ENUM: no initial non vertical facet: C(Phi0) "
                   "does not span and no seed\n";
      throw TerminalException{1};
    }
    while (true) {
      int i_sel = -1;
      for (size_t i = 0; i < l_nv.size(); i++) {
        if (!l_done[i] && (i_sel == -1 ||
                           l_nv[i].l_vert.size() < l_nv[i_sel].l_vert.size())) {
          i_sel = i;
        }
      }
      if (i_sel == -1) {
        break;
      }
      l_done[i_sel] = 1;
      IgusaFace<T> psi = l_nv[i_sel];
      std::vector<IgusaFace<T>> l_theta =
          igusa_face_facets(ctx, bank, psi, l_seed[i_sel], depth + 1);
      std::vector<IgusaFace<T>> l_split =
          igusa_orbit_splitting(ctx, phi, psi, l_theta);
      for (auto &theta : l_split) {
        IgusaFace<T> psi2 = igusa_face_flip(ctx, phi, psi, theta);
        if (psi2.rank == phi.rank) {
          // a vertical facet, already given by the polytope
          continue;
        }
        insert(psi2, theta);
      }
#ifdef DEBUG_IGUSA_FACE_ENUM
      os << "IGUSA_FACE_ENUM: depth=" << depth << " processed "
         << igusa_face_summary(psi, n) << " |facets|=" << l_theta.size()
         << " split=" << l_split.size() << " |non vertical|=" << l_nv.size()
         << "\n";
#endif
    }
    for (auto &psi : l_nv) {
      l_facet.push_back(psi);
    }
  }
#ifdef TIMINGS_IGUSA_FACE_ENUM
  os << "|IGUSA_FACE_ENUM: depth=" << depth << " " << igusa_face_summary(phi, n)
     << " |facets|=" << l_facet.size() << "|=" << time << "\n";
#endif
  bank.map_key[key_phi.key] = bank.l_entry.size();
  bank.l_entry.push_back({phi, key_phi, l_facet});
  return l_facet;
}

template <typename T> struct IgusaFacetEnumeration {
  std::vector<IgusaFace<T>> l_facet;
  // for each facet: the number of orbits of ridges and of flips done
  std::vector<size_t> l_n_ridge;
  std::vector<size_t> l_n_flip;
};

/*
  The orbits of facets of P_n (Corollary COR_full_rank_reached): the facets
  of lower rank are given; for them only the bounded ridges are flipped, the
  other ridges leading to facets of lower rank. For the full rank facets all
  the ridges are flipped.
 */
template <typename T, typename Tint, typename Tgroup>
IgusaFacetEnumeration<T>
igusa_enumerate_facets(IgusaFaceContext<T, Tint, Tgroup> &ctx,
                       std::vector<IgusaPair<T>> const &l_lower,
                       std::vector<IgusaPair<T>> const &l_initial,
                       int max_facet) {
  int n = ctx.n;
  std::ostream &os = ctx.os;
  IgusaFace<T> whole = igusa_whole_face<T>(n);
  IgusaFaceBank<T> bank;
  IgusaFacetEnumeration<T> res;
  std::map<std::string, size_t> map_key;
  bool lower_phase = true;
  auto insert = [&](IgusaFace<T> const &G) -> void {
    std::string key = igusa_face_key<T, Tint, Tgroup>(G, {}, os).key;
    if (map_key.count(key) > 0) {
      return;
    }
    if (G.rank < n && !lower_phase) {
      std::cerr << "IGUSA_FACE_ENUM: a facet of lower rank that is not in "
                   "the given list, rank="
                << G.rank << "\n";
      throw TerminalException{1};
    }
    map_key[key] = res.l_facet.size();
    res.l_facet.push_back(G);
    res.l_n_ridge.push_back(0);
    res.l_n_flip.push_back(0);
#ifdef DEBUG_IGUSA_FACE_ENUM
    os << "IGUSA_FACE_ENUM: new facet " << res.l_facet.size() << " "
       << igusa_face_summary(G, n) << " rhs=" << G.pair.b << "\n";
#endif
  };
  for (auto &pair : l_lower) {
    insert(igusa_face_from_pair(ctx, pair));
  }
  lower_phase = false;
  for (auto &pair : l_initial) {
    insert(igusa_face_from_pair(ctx, pair));
  }
  size_t pos = 0;
  while (pos < res.l_facet.size()) {
    if (max_facet > 0 && static_cast<int>(pos) >= max_facet) {
      break;
    }
    IgusaFace<T> G = res.l_facet[pos];
#ifdef TIMINGS_IGUSA_FACE_ENUM
    MicrosecondTime time;
#endif
    std::vector<IgusaFace<T>> l_ridge = igusa_face_facets(ctx, bank, G, {}, 1);
    size_t n_flip = 0;
    for (auto &R : l_ridge) {
      if (G.rank < n && R.rank < n) {
        continue;
      }
      n_flip++;
      insert(igusa_face_flip(ctx, whole, G, R));
    }
    res.l_n_ridge[pos] = l_ridge.size();
    res.l_n_flip[pos] = n_flip;
#ifdef TIMINGS_IGUSA_FACE_ENUM
    os << "|IGUSA_FACE_ENUM: facet " << pos << " " << igusa_face_summary(G, n)
       << " |ridges|=" << l_ridge.size() << " flips=" << n_flip
       << " |facets|=" << res.l_facet.size() << "|=" << time << "\n";
#endif
    pos++;
  }
  return res;
}

inline FullNamelist NAMELIST_GetStandard_IGUSA_FACE_ENUMERATION() {
  std::map<std::string, SingleBlock> ListBlock;
  {
    std::map<std::string, std::string> ListStringValues;
    std::map<std::string, int> ListIntValues;
    std::map<std::string, bool> ListBoolValues;
    ListStringValues["arithmetic"] = "gmp";
    ListStringValues["IlpMethod"] = "default";
    ListStringValues["OutFile"] = "stderr";
    // The vertex representatives of the P_r, r <= n, all ranks in one
    // ListMatrix file
    ListStringValues["FileVertexRepr"] = "unset";
    // The facets of lower rank of P_n (lifts of the full rank facets of
    // the P_r, r < n): n x n matrices and right hand sides
    ListStringValues["FileLowerFacets"] = "unset";
    ListStringValues["FileLowerRhs"] = "unset";
    // Optional full rank facets to start from
    ListStringValues["FileInitialFacets"] = "unset";
    ListStringValues["FileInitialRhs"] = "unset";
    // Stop after processing that many facets (0: no limit)
    ListIntValues["MaxFacet"] = 0;
    SingleBlock BlockDATA;
    BlockDATA.setListStringValues(ListStringValues);
    BlockDATA.setListIntValues(ListIntValues);
    BlockDATA.setListBoolValues(ListBoolValues);
    ListBlock["DATA"] = BlockDATA;
  }
  return FullNamelist(ListBlock);
}

// clang-format off
#endif  // SRC_IGUSA_IGUSA_FACE_ENUM_H_
// clang-format on
