// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_QUADRATICHECKELARB_H_
#define SRC_GENERA_QUADRATICHECKELARB_H_

// clang-format off
#include "ArbReal.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <optional>
#include <vector>
// clang-format on

/*
  A certified enclosure of L(1 - r, chi_{E/K}) with Arb (flint), for the
  exact rounding of the quadratic masses with E not abelian
  (QuadraticMassNonAbelianE). The method is that of QuadHeckeLNegative,
  Lambda(1 - r) = sum_k a_k F(k / A) with F = G_r + G_{1-r}, every
  approximation being bounded:

  - The series of F at the poles z_j = -delta - 2j of gamma, cut after J
    terms. By Cauchy on the circle |w - z_j| = 1 the j-th residue of
    gamma(w) t^{-w} / (s - w) is at most
      pi^{n(j+1/2)} (2 pi / Gamma(j + 1/2))^n |t^{-w}| / (|s - z_j| - 1),
    as (w + delta)/2 = -j + v with |v| = 1/2 and |Gamma(-j + v)| =
    |Gamma(1 + v)| / (|v| prod_{i <= j} |v - i|) <= sqrt(pi) /
    ((1/2) Gamma(j + 1/2) / sqrt(pi)). With t = e^u, u = x + iy, |t^{-w}|
    <= e^{(delta + 2j) x + |x| + |y|}. Consecutive bounds have the ratio
    pi^n e^{2x} / (j + 1/2)^n, so the tail is summed until the ratio is
    1/2 and closed by a geometric series.
  - The Chebyshev interpolation of F(e^u) of degree D on a segment
    [m - h, m + h]: for F analytic in the Bernstein ellipse of parameter rho
    with |F| <= M there, |F - p_D| <= 4 M rho^{-D} / (rho - 1) (Trefethen,
    Approximation theory and approximation practice, Thm 8.2). M is
    bounded by the majorant series: |gamma(r)| e^{-r x} + |gamma(1-r)|
    e^{(r-1) x} + sum_k |u|^k sum_j |q_kj| e^{(delta + 2j) x} plus the tail
    above, on the ellipse.
  - The terms k > X. As phi > 0, G_s(t) = t^{-s} int_t^infty phi(x) x^s
    dx / x <= gamma(sigma) t^{-sigma} for every sigma >= s, and |a_k| is
    at most the k-th coefficient of zeta(s)^n, so that for sigma > 2
      sum_{k > X} |a_k| F(k/A) <= 2 gamma(sigma) A^sigma X^{2-sigma} zeta(2)^n.
  - The arithmetic: Arb balls.
 */

// An exact upper bound of |x| (a ball of radius 0)
inline ArbReal ArbUpper(ArbReal const &x) {
  arf_t u;
  arf_init(u);
  arb_get_abs_ubound_arf(u, x.v, ArbReal::prec());
  ArbReal r;
  arb_set_arf(r.v, u);
  arf_clear(u);
  return r;
}

inline ArbReal QuadHeckeGammaArb(ArbReal const &s, int n, int delta) {
  ArbReal u = (s + ArbReal(delta)) / ArbReal(2);
  ArbReal g = arb_exp(-u * arb_log(arb_pi())) * arb_gamma(u);
  return arb_pow_si(g, n);
}

struct QuadHeckeKernelArb {
  int n;
  int delta;
  int r;
  ArbReal gs;
  ArbReal g1s;
  // q[k][j]: the coefficient of t^{delta + 2j} (log t)^k
  std::vector<std::vector<ArbReal>> q;
  // the bound of the j-th residue without the factor |t^{-w}|, j > J
  ArbReal res_const(int j) const {
    ArbReal pi = arb_pi();
    ArbReal half(1);
    half /= ArbReal(2);
    ArbReal jh = ArbReal(j) + half;
    ArbReal base = ArbReal(2) * pi / arb_gamma(jh);
    ArbReal dist = ArbReal(1 - r + delta + 2 * j) - ArbReal(1);
    return arb_exp(ArbReal(n) * jh * arb_log(pi)) * arb_pow_si(base, n) *
           ArbReal(2) / dist;
  }
  int J() const { return static_cast<int>(q[0].size()) - 1; }
  // an upper bound of the residues j > J for Re u in [., xe], |Re u| <= xa,
  // |Im u| <= y (+infinity when the bounds are not geometric soon enough)
  ArbReal tail(double xe, double xa, double y) const {
    ArbReal pi = arb_pi();
    ArbReal total(0);
    int j = J() + 1;
    ArbReal ex2 = arb_exp(ArbReal(2) * arb_from_double(xe));
    while (true) {
      ArbReal term =
          res_const(j) * arb_exp(ArbReal(delta + 2 * j) * arb_from_double(xe) +
                                 arb_from_double(xa) + arb_from_double(y));
      // the ratio to the next bound
      ArbReal ratio = arb_pow_si(pi, n) * ex2 /
                      arb_pow_si(ArbReal(j) + ArbReal(1) / ArbReal(2), n);
      ArbReal halfb(1);
      halfb /= ArbReal(2);
      if (arb_le(ratio.v, halfb.v)) {
        total += ArbReal(2) * term;
        return ArbUpper(total);
      }
      total += term;
      j++;
      if (j > J() + 2000) {
        // too far to the right (a large ellipse): no useful bound
        ArbReal inf;
        arb_pos_inf(inf.v);
        return inf;
      }
    }
  }
};

inline QuadHeckeKernelArb MakeQuadHeckeKernelArb(int n, int delta, int r,
                                                 int J) {
  QuadHeckeKernelArb ker;
  ker.n = n;
  ker.delta = delta;
  ker.r = r;
  ArbReal s(r);
  ArbReal s1(1 - r);
  ker.gs = QuadHeckeGammaArb(s, n, delta);
  ker.g1s = QuadHeckeGammaArb(s1, n, delta);
  ArbReal pi = arb_pi();
  ArbReal logpi = arb_log(pi);
  ArbReal euler = arb_euler();
  std::vector<ArbReal> zeta(n + 1);
  for (int k = 2; k <= n; k++) {
    zeta[k] = arb_zeta(static_cast<ulong>(k));
  }
  std::vector<ArbReal> H(n + 1, ArbReal(0));
  ArbReal C = arb_pow_si(ArbReal(2), n);
  ArbReal pin = arb_pow_si(pi, n);
  std::vector<ArbReal> fact(n + 1);
  fact[0] = ArbReal(1);
  for (int k = 1; k <= n; k++) {
    fact[k] = fact[k - 1] * ArbReal(k);
  }
  ker.q.assign(n, std::vector<ArbReal>(J + 1, ArbReal(0)));
  ArbReal half = ArbReal(1) / ArbReal(2);
  for (int j = 0; j <= J; j++) {
    if (j > 0) {
      ArbReal ij = ArbReal(1) / ArbReal(j);
      ArbReal pw = ij;
      for (int k = 1; k <= n; k++) {
        H[k] += pw;
        pw *= ij;
      }
      C *= pin * arb_pow_si(-ij, n);
    }
    std::vector<ArbReal> a(n, ArbReal(0));
    ArbReal hk = half;
    for (int k = 1; k < n; k++) {
      ArbReal sk =
          (k == 1) ? H[1] - euler
                   : (((k % 2 == 0) ? zeta[k] : -zeta[k]) + H[k]) / ArbReal(k);
      a[k] = ArbReal(n) * sk * hk;
      hk *= half;
    }
    if (n > 1) {
      a[1] -= logpi * ArbReal(n) / ArbReal(2);
    }
    std::vector<ArbReal> e(n, ArbReal(0));
    e[0] = ArbReal(1);
    for (int k = 1; k < n; k++) {
      ArbReal sum(0);
      for (int i = 1; i <= k; i++) {
        sum += ArbReal(i) * a[i] * e[k - i];
      }
      e[k] = sum / ArbReal(k);
    }
    ArbReal z(-delta - 2 * j);
    ArbReal i1 = ArbReal(1) / (s - z);
    ArbReal i2 = ArbReal(1) / (s1 - z);
    std::vector<ArbReal> h(n);
    ArbReal p1 = i1, p2 = i2;
    for (int k = 0; k < n; k++) {
      h[k] = p1 + p2;
      p1 *= i1;
      p2 *= i2;
    }
    for (int k = 0; k < n; k++) {
      ArbReal g(0);
      for (int i = 0; i <= n - 1 - k; i++) {
        g += e[i] * h[n - 1 - k - i];
      }
      ArbReal sign((k % 2 == 0) ? 1 : -1);
      ker.q[k][j] = C * g * sign / fact[k];
    }
  }
  return ker;
}

// F(e^u) for a real u, with the tail of the series as error
inline ArbReal QuadHeckeKernelValueArb(QuadHeckeKernelArb const &ker,
                                       ArbReal const &u, double ud) {
  ArbReal t2 = arb_exp(ArbReal(2) * u);
  ArbReal res = ker.gs * arb_exp(-ArbReal(ker.r) * u) +
                ker.g1s * arb_exp(ArbReal(ker.r - 1) * u);
  ArbReal corr(0);
  ArbReal Lk(1);
  int J = ker.J();
  for (int k = 0; k < ker.n; k++) {
    ArbReal poly = ker.q[k][J];
    for (int j = J - 1; j >= 0; j--) {
      poly = poly * t2 + ker.q[k][j];
    }
    corr += Lk * poly;
    Lk *= u;
  }
  if (ker.delta == 1) {
    corr *= arb_exp(u);
  }
  res -= corr;
  res.add_error(ker.tail(ud + 1e-9, std::abs(ud) + 1e-9, 0.0));
  return res;
}

// An upper bound of |F(e^u)| for Re u in [xmin, xmax], |u| <= U, |Im u| <= y
inline ArbReal QuadHeckeMajorantArb(QuadHeckeKernelArb const &ker, double xmin,
                                    double xmax, double U, double y) {
  ArbReal xm = arb_from_double(xmin), xM = arb_from_double(xmax);
  ArbReal Ub = arb_from_double(U);
  ArbReal M = arb_abs(ker.gs) * arb_exp(-ArbReal(ker.r) * xm) +
              arb_abs(ker.g1s) * arb_exp(ArbReal(ker.r - 1) * xM);
  ArbReal t2 = arb_exp(ArbReal(2) * xM);
  ArbReal Uk(1);
  int J = ker.J();
  for (int k = 0; k < ker.n; k++) {
    ArbReal poly = ArbUpper(ker.q[k][J]);
    for (int j = J - 1; j >= 0; j--) {
      poly = poly * t2 + ArbUpper(ker.q[k][j]);
    }
    M += Uk * poly * arb_exp(ArbReal(ker.delta) * xM);
    Uk *= Ub;
  }
  M += ker.tail(xmax, std::max(std::abs(xmin), std::abs(xmax)), y);
  return ArbUpper(M);
}

struct QuadHeckeInterpArb {
  std::vector<double> lo;
  std::vector<double> hi;
  std::vector<ArbReal> mid;
  std::vector<ArbReal> half;
  std::vector<std::vector<ArbReal>> coef;
  // the bound of |F - p| on the segment
  std::vector<ArbReal> err;
};

inline std::optional<QuadHeckeInterpArb>
MakeQuadHeckeInterpArb(QuadHeckeKernelArb const &ker, double u0, double u1,
                       ArbReal const &tol, int D) {
  QuadHeckeInterpArb res;
  ArbReal pi = arb_pi();
  // cos(pi k (2i + 1) / (2 (D + 1)))
  std::vector<std::vector<ArbReal>> cs(D + 1, std::vector<ArbReal>(D + 1));
  for (int k = 0; k <= D; k++) {
    for (int i = 0; i <= D; i++) {
      cs[k][i] = arb_cos(pi * ArbReal(k * (2 * i + 1)) / ArbReal(2 * (D + 1)));
    }
  }
  std::vector<double> rhos{1.5, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64, 128};
  bool ok = true;
  std::function<void(double, double, int)> build = [&](double a, double b,
                                                       int depth) -> void {
    if (!ok) {
      return;
    }
    // slightly enlarged, for the rounding of u_k when assigning k
    double eps = 1e-9 * (1 + std::abs(a) + std::abs(b));
    double ae = a - eps, be = b + eps;
    double m = (ae + be) / 2, h = (be - ae) / 2;
    // the best bound over rho
    ArbReal best;
    bool have = false;
    for (double rho : rhos) {
      double ca = (rho + 1 / rho) / 2, cb = (rho - 1 / rho) / 2;
      ArbReal M = QuadHeckeMajorantArb(ker, m - h * ca, m + h * ca,
                                       std::abs(m) + h * ca, h * cb);
      ArbReal rb = arb_from_double(rho);
      ArbReal e = ArbReal(4) * M / (arb_pow_si(rb, D) * (rb - ArbReal(1)));
      e = ArbUpper(e);
      if (!have || arb_lt(e.v, best.v)) {
        best = e;
        have = true;
      }
    }
    if (!arb_le(best.v, tol.v)) {
      if (depth > 40) {
        std::cerr << "QHECKE_ARB: no certified interpolation on [" << a << ", "
                  << b << "]\n";
        ok = false;
        return;
      }
      double mm = (a + b) / 2;
      build(a, mm, depth + 1);
      build(mm, b, depth + 1);
      return;
    }
    ArbReal mb = arb_from_double(m), hb = arb_from_double(h);
    std::vector<ArbReal> g(D + 1);
    for (int i = 0; i <= D; i++) {
      // the node m + h cos(pi (2i+1) / (2(D+1))), exactly as a ball
      ArbReal u = mb + hb * cs[1][i];
      g[i] = QuadHeckeKernelValueArb(ker, u, u.mid_d());
    }
    std::vector<ArbReal> c(D + 1);
    for (int k = 0; k <= D; k++) {
      ArbReal sum(0);
      for (int i = 0; i <= D; i++) {
        sum += g[i] * cs[k][i];
      }
      c[k] = sum * ArbReal(2) / ArbReal(D + 1);
    }
    c[0] /= ArbReal(2);
    res.lo.push_back(a);
    res.hi.push_back(b);
    res.mid.push_back(mb);
    res.half.push_back(hb);
    res.coef.push_back(c);
    res.err.push_back(best);
  };
  double step = 0.25;
  for (double a = u0; a < u1; a += step) {
    build(a, std::min(a + step, u1), 0);
  }
  if (!ok) {
    return {};
  }
  return res;
}

// the interpolant p at u (in the segment seg)
inline ArbReal QuadHeckeInterpArbValue(QuadHeckeInterpArb const &I,
                                       ArbReal const &u, size_t seg) {
  std::vector<ArbReal> const &c = I.coef[seg];
  ArbReal x = (u - I.mid[seg]) / I.half[seg];
  ArbReal x2 = ArbReal(2) * x;
  ArbReal b1(0), b2(0);
  for (int k = static_cast<int>(c.size()) - 1; k >= 1; k--) {
    ArbReal b0 = x2 * b1 - b2 + c[k];
    b2 = b1;
    b1 = b0;
  }
  return x * b1 - b2 + c[0];
}

/*
  The ball of L(1 - r, chi) with an absolute radius of about 2^{-abs_target}
  relative to the size of Lambda, every error being bounded. Nothing when
  the number of terms exceeds the limit or the interpolation cannot be
  certified.
 */
template <typename T>
std::optional<ArbReal> QuadHeckeLNegativeArb(QuadHeckeLData<T> &Ld,
                                             int prec_bits) {
  int r = Ld.m / 2;
  int n = Ld.n;
  int delta = Ld.delta;
  ArbPrecisionGuard g0(64);
  ArbReal A2 = arb_from_integer(Ld.A2);
  double logA = 0.5 * arb_log(A2).mid_d();
  double log2_scale =
      r * logA / std::log(2.0) +
      std::log2(QuadHeckeGammaArb(ArbReal(r), n, delta).abs_ubound());
  double target = std::max(20.0, prec_bits - log2_scale + 20);
  // the number of terms from the tail bound, optimized over sigma
  double lz2 = n * std::log(M_PI * M_PI / 6);
  auto log_tail = [&](double logX) -> double {
    double best = 1e300;
    for (double sigma = std::max(3.0, r + 0.0); sigma < 2000; sigma *= 1.02) {
      double lg = n * (-(sigma + delta) / 2 * std::log(M_PI) +
                       std::lgamma((sigma + delta) / 2));
      double v = std::log(2.0) + lg + sigma * logA + (2 - sigma) * logX + lz2;
      best = std::min(best, v);
    }
    return best;
  };
  double goal = -(target + 8) * std::log(2.0);
  double logX = logA;
  while (log_tail(logX) > goal) {
    logX += 0.05;
  }
  double X_double = std::ceil(std::exp(logX));
  const double X_limit = 4.0e8;
  if (X_double > X_limit) {
    std::cerr << "QHECKE_ARB: L(1-r) needs " << X_double
              << " terms, beyond the limit " << X_limit << "\n";
    return {};
  }
  int64_t X = static_cast<int64_t>(X_double);
  double abs_bits = target + std::log2(X_double) + 24;
  double t_max = X_double / std::exp(logA);
  double u0 = -logA - 0.01;
  double u1 = std::log(X_double) - logA + 0.01;
  // J: the residue bounds at the largest Re u of the ellipses
  double xJ = u1 + 0.5;
  int J = 1;
  while (true) {
    double jh = J + 0.5;
    double lt = n * jh * std::log(M_PI) + n * std::log(2 * M_PI) -
                n * std::lgamma(jh) + (delta + 2 * J) * xJ + std::abs(xJ);
    double ratio = n * std::log(M_PI) + 2 * xJ - n * std::log(jh);
    if (J > 8 && lt / std::log(2.0) < -abs_bits - 20 && ratio < -1) {
      break;
    }
    J++;
  }
  double log2_peak =
      n * M_PI * std::pow(t_max * std::exp(0.5), 2.0 / n) / std::log(2.0);
  double log2_Fmax = std::max(0.0, log2_scale) + 8;
  slong wprec =
      static_cast<slong>(abs_bits + std::max(log2_peak, log2_Fmax) + 64);
  int D =
      std::clamp(static_cast<int>((log2_peak + abs_bits) / 3) + 10, 30, 240);
  ArbPrecisionGuard g1(wprec);
  QuadHeckeKernelArb ker = MakeQuadHeckeKernelArb(n, delta, r, J);
  ArbReal tol = arb_pow_si(ArbReal(2), -static_cast<slong>(abs_bits));
  std::optional<QuadHeckeInterpArb> I =
      MakeQuadHeckeInterpArb(ker, u0, u1, tol, D);
  if (!I) {
    return {};
  }
  QuadHeckeLExtend(Ld, X);
  ArbReal logAb = arb_log(arb_from_integer(Ld.A2)) / ArbReal(2);
  ArbReal Lambda(0);
  std::vector<int64_t> abs_sum(I->lo.size(), 0);
  bool inside = true;
  QuadHeckeForEachTerm(Ld, X, [&](int64_t k, int64_t ak, QuadDD) -> void {
    ArbReal u = arb_log(ArbReal(static_cast<slong>(k))) - logAb;
    double ud = u.mid_d();
    size_t seg = static_cast<size_t>(
        std::upper_bound(I->lo.begin(), I->lo.end(), ud) - I->lo.begin());
    seg = (seg == 0) ? 0 : seg - 1;
    if (ud < I->lo[seg] - 1e-9 || ud > I->hi[seg] + 1e-9) {
      inside = false;
      return;
    }
    Lambda +=
        ArbReal(static_cast<slong>(ak)) * QuadHeckeInterpArbValue(*I, u, seg);
    abs_sum[seg] += (ak < 0) ? -ak : ak;
  });
  if (!inside) {
    std::cerr << "QHECKE_ARB: a term outside the segments\n";
    return {};
  }
  // the interpolation errors
  ArbReal err(0);
  for (size_t i = 0; i < abs_sum.size(); i++) {
    err += ArbReal(static_cast<slong>(abs_sum[i])) * I->err[i];
  }
  // the terms k > X, with sigma optimized as above
  {
    double best = 1e300, sbest = 3;
    for (double sigma = std::max(3.0, r + 0.0); sigma < 2000; sigma *= 1.02) {
      double lg = n * (-(sigma + delta) / 2 * std::log(M_PI) +
                       std::lgamma((sigma + delta) / 2));
      double v = lg + sigma * logA + (2 - sigma) * std::log(X_double);
      if (v < best) {
        best = v;
        sbest = sigma;
      }
    }
    ArbReal sg = arb_from_double(sbest);
    ArbReal tail =
        ArbReal(2) * QuadHeckeGammaArb(sg, n, delta) * arb_exp(sg * logAb) *
        arb_exp((ArbReal(2) - sg) * arb_log(ArbReal(static_cast<slong>(X)))) *
        arb_pow_si(arb_zeta(2), n);
    err += ArbUpper(tail);
  }
  Lambda.add_error(ArbUpper(err));
  ArbReal s1(1 - r);
  ArbReal res =
      Lambda / (arb_exp(logAb * s1) * QuadHeckeGammaArb(s1, n, delta));
#ifdef DEBUG_QHECKE
  std::cerr << "QHECKE_ARB: X=" << X << " J=" << J << " D=" << D
            << " segments=" << I->lo.size() << " wprec=" << wprec
            << " L(1-r)=" << res.str(30) << "\n";
#endif
  return res;
}

// clang-format off
#endif  // SRC_GENERA_QUADRATICHECKELARB_H_
// clang-format on
