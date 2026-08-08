// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "PeriodicDelaunay.h"
// clang-format on

/*
  The coset-aware closest vector computation is checked against a brute
  force enumeration over a box of translations: same minimal norm and same
  set of realizing points.
 */

using T = mpq_class;
using Tint = mpz_class;

static uint64_t lcg_state = 987654321;
int rand_int(int bound) {
  lcg_state = 6364136223846793005ULL * lcg_state + 1442695040888963407ULL;
  return static_cast<int>((lcg_state >> 33) % bound);
}

void check(bool test, std::string const &context) {
  if (!test) {
    std::cerr << "TEST_PeriodicCVP failed: " << context << "\n";
    throw TerminalException{1};
  }
}

// Brute force over the translations z in [-B, B]^n of every coset.
std::pair<T, std::set<MyVector<Tint>>>
BruteForce(MyMatrix<T> const &GramMat, PeriodicPointSet<Tint> const &pps,
           MyVector<T> const &eV, int B) {
  int n = eV.size();
  int m = pps.cosets_num.rows();
  std::optional<T> min_norm;
  std::set<MyVector<Tint>> set_vect;
  std::vector<int> z(n, -B);
  while (true) {
    for (int k = 0; k < m; k++) {
      MyVector<Tint> u(n);
      for (int j = 0; j < n; j++) {
        u(j) = pps.cosets_num(k, j) + pps.N * Tint(z[j]);
      }
      T norm = PeriodicNormDiff(GramMat, pps.N, u, eV);
      if (!min_norm || norm < *min_norm) {
        min_norm = norm;
        set_vect.clear();
      }
      if (norm == *min_norm) {
        set_vect.insert(u);
      }
    }
    int pos = 0;
    while (pos < n) {
      z[pos]++;
      if (z[pos] <= B) {
        break;
      }
      z[pos] = -B;
      pos++;
    }
    if (pos == n) {
      break;
    }
  }
  return {*min_norm, std::move(set_vect)};
}

/*
  Direct check that the transformation maps the periodic point set into
  itself: the images of the cosets translated over a box must all be
  points of the set. Independent of the congruence reasoning of
  PeriodicCosetPermutation.
 */
bool BruteForcePreserved(PeriodicPointSet<Tint> const &pps,
                         PeriodicAffineTransform<Tint> const &x, int B) {
  int n = pps.cosets_num.cols();
  int m = pps.cosets_num.rows();
  T N_T = UniversalScalarConversion<T, Tint>(pps.N);
  T d_T = UniversalScalarConversion<T, Tint>(x.d);
  std::vector<int> z(n, -B);
  while (true) {
    for (int k = 0; k < m; k++) {
      // The point c_k + z, then its image under the transformation.
      MyVector<T> p(n);
      for (int j = 0; j < n; j++) {
        T num = UniversalScalarConversion<T, Tint>(pps.cosets_num(k, j));
        p(j) = num / N_T + T(z[j]);
      }
      MyVector<T> img(n);
      for (int j = 0; j < n; j++) {
        T eSum = UniversalScalarConversion<T, Tint>(x.w(j)) / d_T;
        for (int i = 0; i < n; i++) {
          eSum += p(i) * UniversalScalarConversion<T, Tint>(x.A(i, j));
        }
        img(j) = eSum;
      }
      // The image must have integral numerators over N and sit on a coset.
      MyVector<Tint> u(n);
      for (int j = 0; j < n; j++) {
        T scaled = img(j) * N_T;
        if (!IsInteger(scaled)) {
          return false;
        }
        u(j) = UniversalScalarConversion<Tint, T>(scaled);
      }
      if (!GetCosetIndex(pps, u)) {
        return false;
      }
    }
    int pos = 0;
    while (pos < n) {
      z[pos]++;
      if (z[pos] <= B) {
        break;
      }
      z[pos] = -B;
      pos++;
    }
    if (pos == n) {
      break;
    }
  }
  return true;
}

int main() {
  try {
    std::ostream &os = std::cerr;
    for (int n = 2; n <= 3; n++) {
      for (int iter = 0; iter < 6; iter++) {
        // A positive definite form, diagonally dominant by construction.
        MyMatrix<T> GramMat = ZeroMatrix<T>(n, n);
        for (int i = 0; i < n; i++) {
          for (int j = i + 1; j < n; j++) {
            T val(rand_int(5) - 2);
            GramMat(i, j) = val;
            GramMat(j, i) = val;
          }
        }
        for (int i = 0; i < n; i++) {
          GramMat(i, i) = T(4 * n + rand_int(4));
        }
        // Cosets over a denominator that varies with the iteration.
        int den = 2 + rand_int(3);
        int m = 1 + rand_int(3);
        MyMatrix<T> Cosets(m, n);
        for (int k = 0; k < m; k++) {
          for (int j = 0; j < n; j++) {
            Cosets(k, j) = T(rand_int(den)) / T(den);
          }
        }
        // Random cosets often form a group, and such a set is a lattice
        // rather than a genuinely periodic one: those are skipped.
        std::optional<PeriodicPointSet<Tint>> opt_pps =
            PeriodicPointSetFromRational_Opt<Tint>(Cosets);
        if (!opt_pps) {
          continue;
        }
        PeriodicPointSet<Tint> pps = *opt_pps;
        CVPSolver<T, Tint> solver(GramMat, os);
        for (int i_pt = 0; i_pt < 4; i_pt++) {
          MyVector<T> eV(n);
          for (int j = 0; j < n; j++) {
            eV(j) = T(rand_int(21) - 10) / T(3);
          }
          PeriodicCVPResult<T, Tint> res =
              PeriodicClosestVectors(solver, pps, eV);
          std::pair<T, std::set<MyVector<Tint>>> bf =
              BruteForce(GramMat, pps, eV, 4);
          check(res.TheNorm == bf.first, "minimal norm against brute force");
          std::set<MyVector<Tint>> set_res;
          for (int i = 0; i < res.ListVectScaled.rows(); i++) {
            set_res.insert(GetMatrixRow(res.ListVectScaled, i));
          }
          check(set_res == bf.second,
                "set of closest points against brute force");
          // Every returned point realizes the announced norm.
          for (auto &u : set_res) {
            check(PeriodicNormDiff(GramMat, pps.N, u, eV) == res.TheNorm,
                  "returned point realizes the norm");
          }
        }
      }
    }
    // The coset-permutation predicate against the direct membership test.
    // The transformations are drawn so that both answers occur.
    {
      int n = 2;
      MyMatrix<T> Cosets(2, n);
      Cosets(0, 0) = 0;
      Cosets(0, 1) = 0;
      Cosets(1, 0) = T(1) / T(3);
      Cosets(1, 1) = 0;
      PeriodicPointSet<Tint> pps =
          PeriodicPointSetFromRational<Tint>(Cosets);
      size_t n_pres = 0, n_not = 0;
      for (int iter = 0; iter < 200; iter++) {
        MyMatrix<Tint> A = IdentityMat<Tint>(n);
        for (int k = 0; k < 4; k++) {
          int i = rand_int(n), j = rand_int(n);
          if (i == j) {
            for (int u = 0; u < n; u++) {
              A(i, u) = -A(i, u);
            }
          } else {
            Tint c(rand_int(3) - 1);
            for (int u = 0; u < n; u++) {
              A(i, u) += c * A(j, u);
            }
          }
        }
        MyVector<Tint> w(n);
        for (int j = 0; j < n; j++) {
          w(j) = Tint(rand_int(9) - 4);
        }
        // Denominators 2 (the natural one) and 3 (never preserving here).
        Tint d(iter % 3 == 0 ? 3 : 2);
        PeriodicAffineTransform<Tint> x{A, w, d};
        bool pred = IsPeriodicPointSetPreserved(pps, x);
        bool bf = BruteForcePreserved(pps, x, 3);
        check(pred == bf, "coset predicate against direct membership");
        if (pred) {
          n_pres++;
        } else {
          n_not++;
        }
        if (pred) {
          // The announced permutation is the one actually realized.
          std::vector<size_t> sigma = *PeriodicCosetPermutation(pps, x);
          for (size_t k = 0; k < sigma.size(); k++) {
            MyVector<Tint> img(n);
            for (int j = 0; j < n; j++) {
              Tint eSum = pps.N * x.w(j) / x.d;
              for (int i = 0; i < n; i++) {
                eSum += pps.cosets_num(k, i) * x.A(i, j);
              }
              img(j) = eSum;
            }
            std::optional<size_t> opt = GetCosetIndex(pps, img);
            check(opt.has_value() && *opt == sigma[k],
                  "announced coset permutation is realized");
          }
        }
      }
      check(n_pres > 0 && n_not > 0,
            "both preserving and non preserving cases occur");
      std::cerr << "coset predicate: preserved=" << n_pres
                << " rejected=" << n_not << "\n";
    }
    std::cerr << "Normal termination of TEST_PeriodicCVP\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of TEST_PeriodicCVP\n";
    exit(e.eVal);
  }
}
