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
        PeriodicPointSet<Tint> pps =
            PeriodicPointSetFromRational<Tint>(Cosets);
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
    std::cerr << "Normal termination of TEST_PeriodicCVP\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of TEST_PeriodicCVP\n";
    exit(e.eVal);
  }
}
