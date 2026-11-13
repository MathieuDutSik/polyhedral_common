// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "QuaternionicMass.h"
#include <chrono>
#include <random>
#include <vector>
// clang-format on

/*
  The order of the image of O(A) (or of the isometries commuting with an
  order acting) modulo p^c, from the certified generators of
  LocalOrthogonalImageGenerators (with the stabilizer chain modulo p, the
  modulo 4 check at 2, the independence pruning, the affine filter of the
  lifts and the layers), against the listing of LocalStabilizerImage, an
  independent enumeration:
  - random positive definite Gram matrices of rank 2 to 4 with p | det, p
    = 2, 3, 5, c = 1, 2 (the radical for e = 1, and e >= 2);
  - Q8 acting on the Hurwitz order D4 and on the Lipschitz order 2 Z^4 at
    p = 2, c = 2, 3 (the order acting, its commutant and the forms of the
    trace dual).
 */

using T = mpq_class;
using Tint = mpz_class;

struct Case {
  MyMatrix<T> A;
  std::vector<MyMatrix<T>> Om;   // the order acting (columns)
  std::vector<MyMatrix<T>> Lam;  // its integral commutant (columns)
  bool with_forms;
};

// e, the residues of A and of p^e A^{-1} modulo p^W
void Prepare(MyMatrix<T> const &A, int64_t p, int &e, ModMat &A_mod,
             ModMat &Be) {
  int m = A.rows();
  MyMatrix<T> Ainv = Inverse(A);
  e = 0;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      if (Ainv(i, j) != 0) {
        e = std::max(e, ValuationOfRational(T(GetDenominator(Ainv(i, j))), p));
      }
    }
  }
  int W = MultiplicityNewtonPrecision(p);
  Tint pW = 1;
  for (int i = 0; i < W; i++) {
    pW *= static_cast<long>(p);
  }
  T pe = 1;
  for (int i = 0; i < e; i++) {
    pe *= static_cast<long>(p);
  }
  A_mod.assign(m, std::vector<int64_t>(m));
  Be.assign(m, std::vector<int64_t>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      A_mod[i][j] = UniversalScalarConversion<int64_t, T>(A(i, j));
      T x = Ainv(i, j) * pe;
      Tint den = GetDenominator_z(x);
      Tint inv = mod_inv(ResInt(den, pW), pW);
      Tint v = ResInt(Tint(GetNumerator_z(x) * inv), pW);
      if (v < 0) {
        v += pW;
      }
      Be[i][j] = UniversalScalarConversion<int64_t, Tint>(v);
    }
  }
}

ModMat ToMod(MyMatrix<T> const &Y) {
  int m = Y.rows();
  ModMat Ym(m, std::vector<int64_t>(m));
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      Ym[i][j] = UniversalScalarConversion<int64_t, T>(Y(i, j));
    }
  }
  return Ym;
}

// 1 when the two orders agree, 0 when they differ, -1 when a method is
// beyond its bounds
int Compare(Case const &cs, int64_t p, int c) {
  int m = cs.A.rows();
  int e;
  ModMat A_mod, Be;
  Prepare(cs.A, p, e, A_mod, Be);
  std::vector<ModMat> Om, Lam, Forms;
  for (auto const &X : cs.Om) {
    Om.push_back(ToMod(X));
  }
  for (auto const &X : cs.Lam) {
    Lam.push_back(ToMod(X));
  }
  if (cs.with_forms) {
    std::vector<MyMatrix<T>> Ls;
    for (auto const &W2 : cs.Om) {
      Ls.push_back(MyMatrix<T>(W2.transpose() * cs.A));
    }
    for (auto const &N : QuatIntegralPoints<T, Tint>(Ls, m)) {
      Forms.push_back(ToMod(N));
    }
  }
  // bounds small enough for a test that runs in the CI: the cases beyond
  // them are counted as skipped
  MultiplicityFactorBounds bounds;
  bounds.max_node = 200000;
  bounds.max_solution = 200000;
  LocalImageOrder ord;
  auto t0 = std::chrono::steady_clock::now();
  std::vector<MyMatrix<Tint>> LamInt;
  for (auto const &X : cs.Lam) {
    LamInt.push_back(UniversalMatrixConversion<Tint, T>(X));
  }
  int64_t pW = ResiduePrimePower(p, MultiplicityNewtonPrecision(p));
  std::optional<std::vector<ModMat>> Lie = ResidueMatrices(
      IntegralLieLatticeColumns<T, Tint>(cs.A, LamInt), pW);
  std::optional<std::vector<ModMat>> gens = LocalOrthogonalImageGenerators(
      A_mod, Be, Om, Lam, Forms, Lie, p, c, e, bounds, &ord);
  auto t1 = std::chrono::steady_clock::now();
  // the listing, on the rows: the transposes of the commutant, or all the
  // matrices
  std::vector<ModMat> basis;
  if (cs.Lam.empty()) {
    for (int a = 0; a < m; a++) {
      for (int b = 0; b < m; b++) {
        ModMat E(m, std::vector<int64_t>(m, 0));
        E[a][b] = 1;
        basis.push_back(E);
      }
    }
  } else {
    for (auto const &X : cs.Lam) {
      basis.push_back(ToMod(MyMatrix<T>(X.transpose())));
    }
  }
  std::optional<std::vector<ModMat>> list =
      LocalStabilizerImage(A_mod, basis, p, c, e, bounds);
  auto t2 = std::chrono::steady_clock::now();
  double s1 = std::chrono::duration<double>(t1 - t0).count();
  double s2 = std::chrono::duration<double>(t2 - t1).count();
  if (s1 + s2 > 2) {
    std::cerr << "slow case p=" << p << " c=" << c << " e=" << e << " m=" << m
              << " generators " << s1 << " s (" << (gens ? "ok" : "none")
              << "), listing " << s2 << " s (" << (list ? "ok" : "none")
              << ")\n";
  }
  if (!gens || !list) {
    return -1;
  }
  Tint n_gen = 1;
  for (size_t x : ord.orbits) {
    n_gen *= static_cast<long>(x);
  }
  for (size_t i = 0; i < ord.layer_exp; i++) {
    n_gen *= static_cast<long>(p);
  }
  Tint n_list = static_cast<long>(list->size());
  if (n_gen != n_list) {
    std::cerr << "MISMATCH p=" << p << " c=" << c << " e=" << e << " m=" << m
              << " generators " << n_gen << " listing " << n_list << "\nA=\n";
    WriteMatrix(std::cerr, cs.A);
    return 0;
  }
  return 1;
}

MyMatrix<T> FromInt(std::vector<std::vector<int>> const &M) {
  int m = M.size();
  MyMatrix<T> R(m, m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      R(i, j) = M[i][j];
    }
  }
  return R;
}

int main() {
  size_t n_ok = 0, n_skip = 0, n_bad = 0;
  auto record = [&](int res) {
    if (res == 1) {
      n_ok++;
    } else if (res == 0) {
      n_bad++;
    } else {
      n_skip++;
    }
  };
  // random Gram matrices
  std::mt19937 rng(23);
  std::uniform_int_distribution<int> dist(-2, 2);
  for (int64_t p : {2, 3, 5}) {
    for (int m = 2; m <= 4; m++) {
      int found = 0;
      for (int attempt = 0; attempt < 4000 && found < 4; attempt++) {
        MyMatrix<T> M(m, m);
        for (int i = 0; i < m; i++) {
          for (int j = 0; j < m; j++) {
            M(i, j) = dist(rng);
          }
        }
        MyMatrix<T> A = M.transpose() * M;
        for (int i = 0; i < m; i++) {
          A(i, i) += 1;
        }
        if (p == 2) {
          A = T(2) * A;
        }
        T det = DeterminantMat(A);
        if (det == 0 || ValuationOfRational(det, p) == 0) {
          continue;
        }
        found++;
        for (int c = 1; c <= 2; c++) {
          if (p == 2 && c == 1) {
            continue;
          }
          record(Compare({A, {}, {}, false}, p, c));
        }
      }
    }
  }
  // Q8 on the Hurwitz order (basis h = (1 + i + j + k)/2, i, j, k) and on
  // the Lipschitz order (1, i, j, k)
  std::vector<std::vector<std::vector<int>>> HurL = {
      {{-1, -1, -1, -1}, {1, 1, 0, 1}, {1, 1, 1, 0}, {1, 0, 1, 1}},
      {{-1, -2, 0, 0}, {1, 1, 0, 0}, {0, 1, 0, -1}, {1, 1, 1, 0}},
      {{-1, 0, -2, 0}, {1, 0, 1, 1}, {1, 0, 1, 0}, {0, -1, 1, 0}},
      {{-1, 0, 0, -2}, {0, 0, -1, 1}, {1, 1, 0, 1}, {1, 0, 0, 1}}};
  std::vector<std::vector<std::vector<int>>> HurR = {
      {{-1, -1, -1, -1}, {1, 1, 1, 0}, {1, 0, 1, 1}, {1, 1, 0, 1}},
      {{-1, -2, 0, 0}, {1, 1, 0, 0}, {1, 1, 0, 1}, {0, 1, -1, 0}},
      {{-1, 0, -2, 0}, {0, 0, 1, -1}, {1, 0, 1, 0}, {1, 1, 1, 0}},
      {{-1, 0, 0, -2}, {1, 0, 1, 1}, {0, -1, 0, 1}, {1, 0, 0, 1}}};
  Case hur;
  hur.A = FromInt({{2, 1, 1, 1}, {1, 2, 0, 0}, {1, 0, 2, 0}, {1, 0, 0, 2}});
  for (auto const &X : HurL) {
    hur.Om.push_back(FromInt(X));
  }
  for (auto const &X : HurR) {
    hur.Lam.push_back(FromInt(X));
  }
  hur.with_forms = true;
  // the Lipschitz order: left and right multiplications by 1, i, j, k
  auto lip = [](int w, bool left) -> MyMatrix<T> {
    // quaternion units as signed permutations of (1, i, j, k)
    int tab[4][4] = {{0, 1, 2, 3}, {1, 0, 3, 2}, {2, 3, 0, 1}, {3, 2, 1, 0}};
    int sgn_left[4][4] = {
        {1, 1, 1, 1}, {1, -1, 1, -1}, {1, -1, -1, 1}, {1, 1, -1, -1}};
    MyMatrix<T> R = ZeroMatrix<T>(4, 4);
    for (int b = 0; b < 4; b++) {
      // w * e_b (left) or e_b * w (right)
      int x = left ? w : b;
      int y = left ? b : w;
      R(tab[x][y], b) = sgn_left[x][y];
    }
    return R;
  };
  Case lipc;
  lipc.A = T(2) * IdentityMat<T>(4);
  for (int w = 0; w < 4; w++) {
    lipc.Om.push_back(lip(w, true));
    lipc.Lam.push_back(lip(w, false));
  }
  lipc.with_forms = true;
  for (Case const &cs : {hur, lipc}) {
    for (auto const &X : cs.Om) {
      for (auto const &Y : cs.Lam) {
        if (MyMatrix<T>(X * Y) != MyMatrix<T>(Y * X)) {
          std::cerr << "the commutant does not commute\n";
          return 1;
        }
      }
    }
    for (int c = 2; c <= 3; c++) {
      int res = Compare(cs, 2, c);
      std::cerr << "Q8 on " << (cs.A(0, 1) == 0 ? "2 Z^4" : "D4") << " c=" << c
                << ": " << (res == 1 ? "agree" : (res == 0 ? "DIFFER" : "skipped"))
                << "\n";
      record(res);
    }
  }
  std::cerr << "Test_LocalImage: " << n_ok << " agree, " << n_bad
            << " differ, " << n_skip << " beyond the bounds\n";
  if (n_bad > 0 || n_ok == 0) {
    return 1;
  }
  return 0;
}
