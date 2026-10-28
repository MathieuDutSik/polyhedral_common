// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_REALIFIEDRING_H_
#define SRC_GENERA_REALIFIEDRING_H_

// clang-format off
#include "ImagQuadRing.h"
#include "MAT_Matrix.h"
#include "Basic_string.h"
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
// clang-format on

/*
  A lattice over R = Z or R = O_K = Z[omega] (omega^2 = s omega - p, K
  imaginary quadratic) with a form h, realified: the vector (x_1, ..., x_m)
  of R^m, x_j = a_j + b_j omega, is the row vector (a_1, ..., a_m, b_1, ...,
  b_m) of Z^{d m} (d = 1 for Z, 2 for O_K), a matrix of GL_m(R) is the
  matrix of the induced Z-linear map, acting on the right, and the form is
  replaced by the real form B = Tr_{K/Q} h. The group of h is then

    G = { g in GL_{d m}(Z) : g W = W g and g B g^T = B },

  W the matrix of the multiplication by omega. The form h is sesquilinear,
  linear in the first argument: h(u, v) = u h v^* on the row vectors.

  This is shared by src_projective_systems (the cell complexes) and
  src_genera (the genera of Hermitian lattices and the algebraic modular
  forms).
 */

template <typename T, typename Tint> struct RealifiedRing {
  int d = 1;
  int m = 0;
  long s = 0;
  long p = 0;
  MyMatrix<Tint> W;
  std::vector<MyMatrix<Tint>> l_unit;
  std::vector<MyMatrix<Tint>> l_ring_basis;
  MyMatrix<Tint> B;
  int n() const { return d * m; }
  ImagQuadRing<T> quad() const { return ImagQuadRing<T>{T(s), T(p)}; }
  // Whether g is in G: integral, R-linear and preserving B.
  template <typename Tf> bool in_group(MyMatrix<Tf> const &g) const {
    if (!IsIntegralMatrix(g)) {
      return false;
    }
    MyMatrix<Tf> B_f = UniversalMatrixConversion<Tf, Tint>(B);
    if (d == 2) {
      MyMatrix<Tf> W_f = UniversalMatrixConversion<Tf, Tint>(W);
      if (g * W_f != W_f * g) {
        return false;
      }
    }
    return g * B_f * g.transpose() == B_f;
  }
  // The element a + b omega of O_K for the coordinate j of v.
  std::pair<T, T> entry(MyVector<Tint> const &v, int const &j) const {
    return {T(v(j)), T(v(m + j))};
  }
  /*
    The primitive vector of the line of v: v divided by the gcd of its
    entries (over R).
   */
  MyVector<Tint> primitive(MyVector<T> const &v_in) const {
    MyVector<T> v = RemoveFractionVector(v_in);
    if (d == 1) {
      return UniversalVectorConversion<Tint, T>(v);
    }
    ImagQuadRing<T> ring = quad();
    std::vector<std::pair<T, T>> l_x;
    for (int j = 0; j < m; j++) {
      std::pair<T, T> x{v(j), v(m + j)};
      if (!ring.is_zero(x)) {
        l_x.push_back(x);
      }
    }
    std::optional<std::pair<T, T>> opt = ring.gcd(l_x);
    if (!opt) {
      std::cerr << "RING: primitive, the gcd of the entries is not principal\n";
      throw TerminalException{1};
    }
    MyVector<Tint> w(n());
    for (int j = 0; j < m; j++) {
      std::pair<T, T> y = ring.div({v(j), v(m + j)}, *opt);
      w(j) = UniversalScalarConversion<Tint, T>(y.first);
      w(m + j) = UniversalScalarConversion<Tint, T>(y.second);
    }
    return w;
  }
};

/*
  The matrix W of the multiplication by omega on the realified R^m:
  (a, b) -> omega (a + b omega) = -p b + (a + s b) omega.
 */
template <typename Tint>
MyMatrix<Tint> get_omega_matrix(int const &m, long const &s, long const &p) {
  MyMatrix<Tint> W = ZeroMatrix<Tint>(2 * m, 2 * m);
  for (int j = 0; j < m; j++) {
    W(j, m + j) = 1;
    W(m + j, j) = -p;
    W(m + j, m + j) = s;
  }
  return W;
}

/*
  The real form B: for R = Z the form h itself, for R = O_K the trace form
  Tr_{K/Q} h(e, f) on the real basis e_j, omega e_j of O_K^m, h = h_a +
  omega h_b.
 */
template <typename T, typename Tint>
MyMatrix<Tint> get_real_form(int const &m, long const &s, long const &p,
                             MyMatrix<Tint> const &h_a,
                             MyMatrix<Tint> const &h_b) {
  ImagQuadRing<T> ring{T(s), T(p)};
  std::vector<std::pair<T, T>> l_basis{{T(1), T(0)}, {T(0), T(1)}};
  MyMatrix<Tint> B(2 * m, 2 * m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      std::pair<T, T> h_ij{T(h_a(i, j)), T(h_b(i, j))};
      for (int ki = 0; ki < 2; ki++) {
        for (int kj = 0; kj < 2; kj++) {
          std::pair<T, T> val =
              ring.mul(ring.mul(l_basis[ki], h_ij), ring.conj(l_basis[kj]));
          B(ki * m + i, kj * m + j) =
              UniversalScalarConversion<Tint, T>(ring.trace(val));
        }
      }
    }
  }
  return B;
}

/*
  The tokens of a description file, the comments (from # to the end of
  the line) removed.
 */
inline std::vector<std::string> read_description_tokens(std::string const &file) {
  std::ifstream is(file);
  if (!is.good()) {
    std::cerr << "RING: the description file " << file << " cannot be read\n";
    throw TerminalException{1};
  }
  std::vector<std::string> l_tok;
  std::string line;
  while (std::getline(is, line)) {
    size_t pos = line.find('#');
    if (pos != std::string::npos) {
      line = line.substr(0, pos);
    }
    std::istringstream iss(line);
    std::string tok;
    while (iss >> tok) {
      l_tok.push_back(tok);
    }
  }
  return l_tok;
}

template <typename Tint> struct DescriptionTokenReader {
  std::vector<std::string> l_tok;
  size_t pos = 0;
  bool at_end() const { return pos >= l_tok.size(); }
  std::string next() {
    if (pos >= l_tok.size()) {
      std::cerr << "RING: the description file ends too early\n";
      throw TerminalException{1};
    }
    return l_tok[pos++];
  }
  void expect(std::string const &key) {
    std::string tok = next();
    if (tok != key) {
      std::cerr << "RING: expected " << key
                << " in the description file, found " << tok << "\n";
      throw TerminalException{1};
    }
  }
  int next_int() { return ParseScalar<int>(next()); }
  MyMatrix<Tint> next_matrix() {
    int n_row = next_int();
    int n_col = next_int();
    MyMatrix<Tint> M(n_row, n_col);
    for (int i = 0; i < n_row; i++) {
      for (int j = 0; j < n_col; j++) {
        M(i, j) = ParseScalar<Tint>(next());
      }
    }
    return M;
  }
};

/*
  The ring and the form, the head of a description file:

    ring Z                     or   ring QuadraticImag s p
    dimension m
    form [m x m matrix]        or   form [h_a] [h_b]     (h = h_a + omega h_b)
 */
template <typename T, typename Tint>
RealifiedRing<T, Tint> read_realified_ring(DescriptionTokenReader<Tint> &tr) {
  RealifiedRing<T, Tint> ring;
  tr.expect("ring");
  std::string ring_name = tr.next();
  if (ring_name == "Z") {
    ring.d = 1;
  } else if (ring_name == "QuadraticImag") {
    ring.d = 2;
    ring.s = tr.next_int();
    ring.p = tr.next_int();
  } else {
    std::cerr << "RING: the ring should be Z or QuadraticImag, not "
              << ring_name << "\n";
    throw TerminalException{1};
  }
  tr.expect("dimension");
  ring.m = tr.next_int();
  int n = ring.n();
  tr.expect("form");
  if (ring.d == 1) {
    ring.B = tr.next_matrix();
  } else {
    MyMatrix<Tint> h_a = tr.next_matrix();
    MyMatrix<Tint> h_b = tr.next_matrix();
    ring.B = get_real_form<T, Tint>(ring.m, ring.s, ring.p, h_a, h_b);
  }
  if (ring.B.rows() != n || ring.B.cols() != n) {
    std::cerr << "RING: the form should be of dimension " << n << "\n";
    throw TerminalException{1};
  }
  MyMatrix<Tint> id = IdentityMat<Tint>(n);
  ring.l_ring_basis.push_back(id);
  if (ring.d == 1) {
    ring.l_unit = {id, MyMatrix<Tint>(-id)};
  } else {
    ring.W = get_omega_matrix<Tint>(ring.m, ring.s, ring.p);
    ring.l_ring_basis.push_back(ring.W);
    ImagQuadRing<T> iq = ring.quad();
    for (auto &u : iq.units()) {
      MyMatrix<Tint> U = UniversalScalarConversion<Tint, T>(u.first) * id +
                         UniversalScalarConversion<Tint, T>(u.second) * ring.W;
      ring.l_unit.push_back(U);
    }
  }
  return ring;
}

/*
  Elements of R named in a comma separated list: an integer N or, over
  O_K, a + b omega written "a+bw", "a-bw", "bw" (w = omega). Used for the
  levels (a generator of the ideal) and the primes.
 */
struct RingElementName {
  std::string name;
  int64_t a;
  int64_t b;
};

inline std::vector<RingElementName> parse_ring_elements(std::string const &str) {
  std::vector<RingElementName> l;
  for (auto &s : STRING_Split(str, ",")) {
    size_t pos_w = s.find('w');
    if (pos_w == std::string::npos) {
      l.push_back({s, ParseScalar<int64_t>(s), 0});
      continue;
    }
    // The position of the sign separating a and b (not the leading one).
    size_t pos_sign = std::string::npos;
    for (size_t i = 1; i < pos_w; i++) {
      if (s[i] == '+' || s[i] == '-') {
        pos_sign = i;
      }
    }
    int64_t a = 0;
    std::string str_b;
    if (pos_sign == std::string::npos) {
      str_b = s.substr(0, pos_w);
    } else {
      a = ParseScalar<int64_t>(s.substr(0, pos_sign));
      str_b = s.substr(pos_sign, pos_w - pos_sign);
      if (str_b[0] == '+') {
        str_b = str_b.substr(1);
      }
    }
    int64_t b = 1;
    if (str_b == "-") {
      b = -1;
    } else if (!str_b.empty()) {
      b = ParseScalar<int64_t>(str_b);
    }
    l.push_back({s, a, b});
  }
  return l;
}

// clang-format off
#endif  // SRC_GENERA_REALIFIEDRING_H_
// clang-format on
