// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_CMORDER_H_
#define SRC_GENERA_CMORDER_H_

// clang-format off
#include "RealifiedRing.h"
#include "ScalarNumberField.h"
#include "MaximalOrders.h"
#include "Shvec_exact.h"
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_CM_ORDER
#endif

/*
  The order O_E of a CM field E = Q(theta), and the Hermitian lattices over
  it realified (RealifiedRing.h) for any degree d = [E : Q].

  E is given by the minimal polynomial of theta, monic and integral, and
  the complex conjugate of theta on the powers of theta. O_E is the maximal
  order (Round 2 of MaximalOrders.h, from Z[theta]), or Z[theta] itself
  when asked, as the quadratic rings Z[omega] of the description files
  are. The elements are coordinate row vectors on a Z-basis e_0, ..., e_{d-1}
  of O_E with e_0 = 1, and an element acts by its multiplication matrix on
  the right, x -> x M, as in ScalarNumberField.h. The roots of unity are
  the x of O_E with Tr(x conj(x)) = d: that sum of the |sigma(x)|^2 is at
  least d for a nonzero integer, with equality exactly when every
  |sigma(x)| = 1 (Kronecker).

  The realified lattice: the vector (x_1, ..., x_m) of E^m, x_j = sum_a
  c_{a, j} e_a, is the row vector with c_{a, j} at the position a m + j,
  which is the layout of RealifiedRing.h for O_K = Z[omega] (e_0 = 1, e_1 =
  omega). The form is B = Tr_{E/Q} h. A Z-linear map is O_E-linear iff it
  commutes with the multiplication by theta, which generates E, so the
  pair (B, B W^T), W the multiplication by theta, is the family whose
  canonical forms and automorphism groups are those of the Hermitian
  lattice: an E-linear map preserving Tr h preserves h, the trace form
  being nondegenerate.
 */

template <typename T, typename Tint> struct CMOrder {
  NumberField<T> field;
  int deg;
  // the rows: the basis e_a of O_E on the powers of theta
  MyMatrix<T> basis_theta;
  // the multiplication by e_a, and the conjugation, on the basis e
  std::vector<MyMatrix<Tint>> mult;
  MyMatrix<Tint> conj;
  // theta on the basis e
  MyVector<Tint> theta;
  // the trace of e_a
  std::vector<Tint> trace;
  std::vector<MyVector<Tint>> roots_of_unity;
  // The multiplication matrix of the element x (on the basis e).
  MyMatrix<Tint> mult_matrix(MyVector<Tint> const &x) const {
    MyMatrix<Tint> M = ZeroMatrix<Tint>(deg, deg);
    for (int a = 0; a < deg; a++) {
      if (x(a) != 0) {
        M += x(a) * mult[a];
      }
    }
    return M;
  }
  MyVector<Tint> product(MyVector<Tint> const &x,
                         MyVector<Tint> const &y) const {
    return mult_matrix(x).transpose() * y;
  }
  MyVector<Tint> conjugate(MyVector<Tint> const &x) const {
    return conj.transpose() * x;
  }
  Tint trace_of(MyVector<Tint> const &x) const {
    Tint t(0);
    for (int a = 0; a < deg; a++) {
      t += x(a) * trace[a];
    }
    return t;
  }
  // The form Tr(x conj(y)) on the basis e, positive definite.
  MyMatrix<Tint> trace_norm_form() const {
    MyMatrix<Tint> Q(deg, deg);
    for (int a = 0; a < deg; a++) {
      for (int b = 0; b < deg; b++) {
        MyVector<Tint> ea = ZeroVector<Tint>(deg), eb = ZeroVector<Tint>(deg);
        ea(a) = 1;
        eb(b) = 1;
        Q(a, b) = trace_of(product(ea, conjugate(eb)));
      }
    }
    return Q;
  }
};

/*
  A Z-basis of the lattice spanned by the rows of A (full rank, rational)
  whose first vector is v, a primitive vector of it: the coordinates of v
  on the rows are completed to a unimodular matrix.
 */
template <typename T, typename Tint>
MyMatrix<T> CMOrderBasisStartingWith(MyMatrix<T> const &A,
                                     MyVector<T> const &v) {
  int n = A.rows();
  std::optional<MyVector<T>> opt = SolutionMat(A, v);
  if (!opt || !IsIntegralVector(*opt)) {
    std::cerr << "CMORDER: 1 is not in the order\n";
    throw TerminalException{1};
  }
  MyVector<Tint> c = UniversalVectorConversion<Tint, T>(*opt);
  // the coordinates of v, then n - 1 vectors completing them to a basis
  MyMatrix<Tint> compl_c = ComplementToBasis(c);
  MyMatrix<Tint> U(n, n);
  for (int j = 0; j < n; j++) {
    U(0, j) = c(j);
  }
  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < n; j++) {
      U(i + 1, j) = compl_c(i, j);
    }
  }
  if (T_abs(DeterminantMat(U)) != 1) {
    std::cerr << "CMORDER: 1 is not primitive in the order\n";
    throw TerminalException{1};
  }
  return MyMatrix<T>(UniversalMatrixConversion<T, Tint>(U) * A);
}

template <typename T, typename Tint>
CMOrder<T, Tint> CMOrderOfField(NumberField<T> const &K, bool maximal,
                                std::ostream &os) {
  int d = K.deg;
  CMOrder<T, Tint> ord;
  ord.field = K;
  ord.deg = d;
  // Z[theta]: the multiplication matrices of the powers of theta
  std::vector<MyMatrix<T>> O0;
  for (int k = 0; k < d; k++) {
    O0.push_back(EMulMatrix(EZetaPower(k, K), K));
  }
  std::vector<MyMatrix<T>> O = maximal ? CommutativeMaximalOrder<T, Tint>(O0, d)
                                       : O0;
  MyMatrix<T> A(d, d);
  for (int a = 0; a < d; a++) {
    // 1 * M is the element of the matrix M
    for (int i = 0; i < d; i++) {
      A(a, i) = O[a](0, i);
    }
  }
  MyVector<T> one = ZeroVector<T>(d);
  one(0) = 1;
  ord.basis_theta = CMOrderBasisStartingWith<T, Tint>(A, one);
  MyMatrix<T> Binv = Inverse(ord.basis_theta);
  auto to_int = [&](MyMatrix<T> const &M) -> MyMatrix<Tint> {
    if (!IsIntegralMatrix(M)) {
      std::cerr << "CMORDER: the basis is not that of an order\n";
      throw TerminalException{1};
    }
    return UniversalMatrixConversion<Tint, T>(M);
  };
  for (int a = 0; a < d; a++) {
    ElemE<T> e_a(d, T(0));
    for (int i = 0; i < d; i++) {
      e_a[i] = ord.basis_theta(a, i);
    }
    MyMatrix<T> M = EMulMatrix(e_a, K);
    ord.mult.push_back(to_int(MyMatrix<T>(ord.basis_theta * M * Binv)));
  }
  MyMatrix<T> Ctheta(d, d);
  for (int j = 0; j < d; j++) {
    ElemE<T> c = EConj(EZetaPower(j, K), K);
    for (int i = 0; i < d; i++) {
      Ctheta(j, i) = c[i];
    }
  }
  ord.conj = to_int(MyMatrix<T>(ord.basis_theta * Ctheta * Binv));
  MyVector<T> theta_t = ZeroVector<T>(d);
  if (d > 1) {
    theta_t(1) = 1;
  } else {
    theta_t(0) = -UniversalScalarConversion<T, int64_t>(K.Phi[0]);
  }
  MyVector<T> theta_e = Binv.transpose() * theta_t;
  if (!IsIntegralVector(theta_e)) {
    std::cerr << "CMORDER: theta is not in the order\n";
    throw TerminalException{1};
  }
  ord.theta = UniversalVectorConversion<Tint, T>(theta_e);
  for (int a = 0; a < d; a++) {
    Tint t(0);
    for (int i = 0; i < d; i++) {
      t += ord.mult[a](i, i);
    }
    ord.trace.push_back(t);
  }
  // the roots of unity: Tr(x conj(x)) = d
  MyMatrix<Tint> Q = ord.trace_norm_form();
  MyMatrix<T> Q_T = UniversalMatrixConversion<T, Tint>(Q);
  CVPSolver<T, Tint> solver(Q_T, os);
  // one of each pair x, -x
  for (auto &v : solver.fixed_norm_vectors(T(d))) {
    ord.roots_of_unity.push_back(v);
    ord.roots_of_unity.push_back(MyVector<Tint>(-v));
  }
#ifdef SANITY_CHECK_CM_ORDER
  for (int a = 0; a < d; a++) {
    for (int b = 0; b < d; b++) {
      if (ord.mult[a] * ord.mult[b] != ord.mult[b] * ord.mult[a]) {
        std::cerr << "CMORDER: the order is not commutative\n";
        throw TerminalException{1};
      }
    }
  }
  if (ord.conj * ord.conj != IdentityMat<Tint>(d)) {
    std::cerr << "CMORDER: the conjugation is not an involution\n";
    throw TerminalException{1};
  }
  for (int a = 0; a < d; a++) {
    MyVector<Tint> ea = ZeroVector<Tint>(d);
    ea(a) = 1;
    // conj(e_a x) = conj(e_a) conj(x): C M_a = M_{conj(e_a)} C
    if (ord.conj * ord.mult_matrix(ord.conjugate(ea)) !=
        ord.mult[a] * ord.conj) {
      std::cerr << "CMORDER: the conjugation is not a ring morphism\n";
      throw TerminalException{1};
    }
  }
  if (!IsPositiveDefinite(Q_T, os)) {
    std::cerr << "CMORDER: Tr(x conj(x)) is not positive definite, E is not "
                 "CM or totally real\n";
    throw TerminalException{1};
  }
  if (ord.roots_of_unity.size() % 2 != 0 || ord.roots_of_unity.empty()) {
    std::cerr << "CMORDER: the roots of unity are not closed under -1\n";
    throw TerminalException{1};
  }
#endif
  return ord;
}

/*
  The realified multiplication by the element x of O_E on E^m: the entry
  (a m + j, b m + j) is M_x(a, b).
 */
template <typename T, typename Tint>
MyMatrix<Tint> CMOrderRealified(MyMatrix<Tint> const &M, int m) {
  int d = M.rows();
  MyMatrix<Tint> R = ZeroMatrix<Tint>(d * m, d * m);
  for (int a = 0; a < d; a++) {
    for (int b = 0; b < d; b++) {
      if (M(a, b) != 0) {
        for (int j = 0; j < m; j++) {
          R(a * m + j, b * m + j) = M(a, b);
        }
      }
    }
  }
  return R;
}

/*
  The realified ring of E^m with the Hermitian form h = sum_k theta^k h_k
  (h_k on the powers of theta): B(a m + i, b m + j) = Tr(e_a h_ij
  conj(e_b)), which must be integral. W is the multiplication by theta, the
  units are the roots of unity; for d = 2, s and p are the trace and the
  norm of e_1, so that O_E = Z[e_1] is the ring of RealifiedRing.h.
 */
template <typename T, typename Tint>
RealifiedRing<T, Tint> CMOrderRealifiedRing(CMOrder<T, Tint> const &ord, int m,
                                            std::vector<MyMatrix<T>> const &h) {
  int d = ord.deg;
  RealifiedRing<T, Tint> ring;
  ring.d = d;
  ring.m = m;
  if (d == 2) {
    MyVector<Tint> e1 = ZeroVector<Tint>(2);
    e1(1) = 1;
    ring.s = static_cast<long>(
        UniversalScalarConversion<int64_t, Tint>(ord.trace_of(e1)));
    ring.p = static_cast<long>(UniversalScalarConversion<int64_t, Tint>(
        ord.trace_of(ord.product(e1, ord.conjugate(e1))) / Tint(2)));
  }
  ring.W = CMOrderRealified<T, Tint>(ord.mult_matrix(ord.theta), m);
  for (int a = 0; a < d; a++) {
    ring.l_ring_basis.push_back(CMOrderRealified<T, Tint>(ord.mult[a], m));
  }
  ring.Conj = CMOrderRealified<T, Tint>(ord.conj, m);
  for (auto &u : ord.roots_of_unity) {
    ring.l_unit.push_back(CMOrderRealified<T, Tint>(ord.mult_matrix(u), m));
  }
  // the entries h_ij on the basis e
  MyMatrix<T> Binv = Inverse(ord.basis_theta);
  MyMatrix<T> B_T(d * m, d * m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      MyVector<T> h_theta = ZeroVector<T>(d);
      for (int k = 0; k < static_cast<int>(h.size()); k++) {
        h_theta(k) = h[k](i, j);
      }
      MyVector<T> h_e = Binv.transpose() * h_theta;
      // Tr(e_a h conj(e_b)) = Tr(h conj(e_b) e_a), over Q
      for (int a = 0; a < d; a++) {
        for (int b = 0; b < d; b++) {
          MyMatrix<T> Ma = UniversalMatrixConversion<T, Tint>(ord.mult[a]);
          MyVector<Tint> eb = ZeroVector<Tint>(d);
          eb(b) = 1;
          MyMatrix<T> Mcb =
              UniversalMatrixConversion<T, Tint>(ord.mult_matrix(ord.conjugate(eb)));
          MyVector<T> prod = Mcb.transpose() * (Ma.transpose() * h_e);
          T tr(0);
          for (int c = 0; c < d; c++) {
            tr += prod(c) * UniversalScalarConversion<T, Tint>(ord.trace[c]);
          }
          B_T(a * m + i, b * m + j) = tr;
        }
      }
    }
  }
  if (!IsIntegralMatrix(B_T)) {
    std::cerr << "CMORDER: the trace form of h is not integral\n";
    throw TerminalException{1};
  }
  ring.B = UniversalMatrixConversion<Tint, T>(B_T);
  return ring;
}

/*
  The description of a lattice over a CM field, as read_realified_ring
  for Z and the quadratic rings:

    ring NumberField deg c_0 ... c_deg conj t_0 ... t_{deg-1}
    dimension m
    form [h_0] ... [h_{deg-1}]

  with theta a root of c_0 + c_1 x + ... + c_deg x^deg (c_deg = 1), conj(theta)
  = t_0 + t_1 theta + ..., and h = sum_k theta^k h_k (m x m matrices).
  O_E is the maximal order. The other rings are passed to
  read_realified_ring.
 */
template <typename T, typename Tint>
std::pair<RealifiedRing<T, Tint>, std::optional<CMOrder<T, Tint>>>
read_cm_realified_ring(DescriptionTokenReader<Tint> &tr, std::ostream &os) {
  size_t pos0 = tr.pos;
  tr.expect("ring");
  std::string ring_name = tr.next();
  if (ring_name != "NumberField") {
    tr.pos = pos0;
    return {read_realified_ring<T, Tint>(tr), {}};
  }
  NumberField<T> K;
  K.d = 0;
  K.deg = tr.next_int();
  for (int k = 0; k <= K.deg; k++) {
    K.Phi.push_back(ParseScalar<int64_t>(tr.next()));
  }
  if (K.Phi[K.deg] != 1) {
    std::cerr << "CMORDER: the polynomial of theta must be monic\n";
    throw TerminalException{1};
  }
  tr.expect("conj");
  K.conj = ElemE<T>(K.deg, T(0));
  for (int k = 0; k < K.deg; k++) {
    K.conj[k] = ParseScalar<T>(tr.next());
  }
  tr.expect("dimension");
  int m = tr.next_int();
  tr.expect("form");
  std::vector<MyMatrix<T>> h;
  for (int k = 0; k < K.deg; k++) {
    h.push_back(UniversalMatrixConversion<T, Tint>(tr.next_matrix()));
  }
  CMOrder<T, Tint> ord = CMOrderOfField<T, Tint>(K, true, os);
  return {CMOrderRealifiedRing<T, Tint>(ord, m, h), ord};
}

// clang-format off
#endif  // SRC_GENERA_CMORDER_H_
// clang-format on
