// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MultiplicityLocalFactor.h"
#include <random>
#include <set>
#include <vector>
// clang-format on

/*
  HowellFormPrimePower against brute force on random subgroups of
  (Z / p^C)^r: the form only depends on the subgroup (two generating sets of
  one subgroup give the same form, and different subgroups different
  forms), and the order read from the pivots is the order of the subgroup.
 */

std::set<std::vector<int64_t>> BruteSpan(ModMatrix const &gens, int64_t N,
                                         size_t r) {
  std::set<std::vector<int64_t>> seen{std::vector<int64_t>(r, 0)};
  std::vector<std::vector<int64_t>> queue{std::vector<int64_t>(r, 0)};
  while (!queue.empty()) {
    std::vector<int64_t> x = queue.back();
    queue.pop_back();
    for (auto const &g : gens) {
      std::vector<int64_t> y(r);
      for (size_t a = 0; a < r; a++) {
        y[a] = (x[a] + g[a]) % N;
      }
      if (seen.insert(y).second) {
        queue.push_back(y);
      }
    }
  }
  return seen;
}

int main() {
  std::mt19937 rng(17);
  size_t n_test = 0;
  for (int64_t p : {2, 3, 5}) {
    for (int C = 1; C <= 3; C++) {
      for (size_t r = 1; r <= 4; r++) {
        int64_t N = MyPow<int64_t>(p, C);
        double sz = std::pow(static_cast<double>(N), static_cast<double>(r));
        if (sz > 20000) {
          continue;
        }
        for (int trial = 0; trial < 40; trial++) {
          auto random_gens = [&](size_t k) -> ModMatrix {
            ModMatrix g;
            for (size_t i = 0; i < k; i++) {
              std::vector<int64_t> row(r);
              for (size_t a = 0; a < r; a++) {
                // often divisible by p, to get the interesting subgroups
                int64_t x = static_cast<int64_t>(rng() % N);
                if (rng() % 2 == 0) {
                  x = (x * p) % N;
                }
                row[a] = x;
              }
              g.push_back(row);
            }
            return g;
          };
          ModMatrix g1 = random_gens(1 + rng() % 3);
          std::set<std::vector<int64_t>> S1 = BruteSpan(g1, N, r);
          // another generating set of the same subgroup: random elements
          // of S1 together with g1 shuffled and combined
          std::vector<std::vector<int64_t>> elts(S1.begin(), S1.end());
          ModMatrix g2;
          for (size_t i = 0; i < g1.size(); i++) {
            std::vector<int64_t> row = g1[(i + 1) % g1.size()];
            std::vector<int64_t> const &z = elts[rng() % elts.size()];
            // row + z is in S1 and, with z kept, generates S1 again
            std::vector<int64_t> w(r);
            for (size_t a = 0; a < r; a++) {
              w[a] = (row[a] + z[a]) % N;
            }
            g2.push_back(w);
            g2.push_back(z);
          }
          for (size_t i = 0; i < g1.size(); i++) {
            g2.push_back(g1[i]);
          }
          ModMatrix H1 = HowellFormPrimePower(g1, p, C);
          ModMatrix H2 = HowellFormPrimePower(g2, p, C);
          if (H1 != H2) {
            std::cerr << "Test_HowellForm: two forms for one subgroup (p=" << p
                      << " C=" << C << " r=" << r << ")\n";
            return 1;
          }
          // the order from the pivots
          double ord = 1;
          for (auto const &row : H1) {
            for (size_t a = 0; a < r; a++) {
              if (row[a] != 0) {
                int v = 0;
                for (int64_t x = row[a]; x % p == 0; x /= p) {
                  v++;
                }
                ord *= std::pow(static_cast<double>(p), C - v);
                break;
              }
            }
          }
          if (static_cast<size_t>(std::llround(ord)) != S1.size()) {
            std::cerr << "Test_HowellForm: order " << ord << " against "
                      << S1.size() << " (p=" << p << " C=" << C << " r=" << r
                      << ")\n";
            return 1;
          }
          // a random other subgroup: same form iff same subgroup
          ModMatrix g3 = random_gens(1 + rng() % 3);
          std::set<std::vector<int64_t>> S3 = BruteSpan(g3, N, r);
          bool same_form = (HowellFormPrimePower(g3, p, C) == H1);
          if (same_form != (S3 == S1)) {
            std::cerr << "Test_HowellForm: forms and subgroups disagree (p="
                      << p << " C=" << C << " r=" << r << ")\n";
            return 1;
          }
          n_test++;
        }
      }
    }
  }
  std::cout << "Test_HowellForm: " << n_test << " random subgroups, all good\n";
  return 0;
}
