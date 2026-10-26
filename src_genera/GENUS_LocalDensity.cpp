// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "LocalDensity.h"
#include "LocalDensityBadPrimes.h"
#include "Arithmetic_dispatch_types.h"
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
// clang-format on

/*
  Exact evaluation of the pieces of Siegel's representation formula, one
  task per piece. The input file is a count followed by that many records;
  the output is a GAP list (return [ ... ];) with one value per record.

    bernoulli               record: k          value: B_k (B_1 = -1/2)
    generalized_bernoulli   record: k D        value: B_{k, chi_D}
    euler_value             record: n m detS detT p
                            value: the good-prime density at p, from its
                                   Euler factors
    generic_part            record: n m detS detT
                            value: alpha_oo times all good primes,
                                   without the constant c = 2
    siegel_series           record: p m S (m*m entries, S even)
                            value: the coefficients [c_0, c_1, ...] of
                                   Katsurada's polynomial F_p(S/2, X)
    alpha_unimodular        record: p k xi_T m S
                            value: alpha_p(S, T) for T unimodular at p of
                                   rank 2k and character xi_T at p
    rank_one_density        record: p n T a      (T even, n*n entries)
                            value: alpha_p(a, T), Q(x) = x^t T x / 2
    unit_chain              record: p n T m S    (p odd)
                            value: alpha_p(S, T), or fail when S has two
                                   or more non-unit entries at p
    counting                record: p n T m S    (p odd)
                            value: alpha_p(S, T) by counting on all of S at
                                   once (no unit chain), or fail beyond the
                                   counting budget
    reduction               record: p n T m S    (p odd)
                            value: alpha_p(S, T) by the unit chain and the
                                   reduction formula
    mass_ratio              record: n T m S      (n even, det T odd)
                            value: M_S / mass of the genus of T, or fail
                                   outside the coverage of the unit chain
 */

template <typename T> MyMatrix<T> ReadSquareMatrixEntries(std::istream &is, int m) {
  MyMatrix<T> S(m, m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      long val;
      is >> val;
      S(i, j) = T(val);
    }
  }
  return S;
}

template <typename T> std::string PolynomialToGAP(UnivariatePolynomial<T> const &P) {
  std::ostringstream oss;
  oss << "[";
  for (size_t i = 0; i < P.coeff.size(); i++) {
    oss << (i > 0 ? ", " : "") << P.coeff[i];
  }
  oss << "]";
  return oss.str();
}

template <typename T, typename Tint>
void ProcessTask(std::string const &task, std::string const &FileIn,
                 std::ostream &os) {
  std::ifstream is(FileIn);
  if (!is.good()) {
    std::cerr << "GENUS_LocalDensity: cannot read " << FileIn << "\n";
    throw TerminalException{1};
  }
  size_t n_rec;
  is >> n_rec;
  std::vector<std::string> ListVal;
  LocalDensityCache<T> cache;
  std::map<std::vector<T>, GenusDensityData> genus_data;
  auto push = [&](T const &val) -> void {
    std::ostringstream oss;
    oss << val;
    ListVal.push_back(oss.str());
  };
  for (size_t i_rec = 0; i_rec < n_rec; i_rec++) {
    if (task == "bernoulli") {
      int k;
      is >> k;
      push(BernoulliNumbers<T>(k)[k]);
    } else if (task == "generalized_bernoulli") {
      int k;
      int64_t D;
      is >> k >> D;
      push(GeneralizedBernoulli<T>(k, D));
    } else if (task == "euler_value") {
      int n, m;
      int64_t detS, detT, p;
      is >> n >> m >> detS >> detT >> p;
      push(EulerValueAtPrime<T>(GoodPrimeEulerTerms(n, m, detS, detT), p));
    } else if (task == "generic_part") {
      int n, m;
      int64_t detS, detT;
      is >> n >> m >> detS >> detT;
      push(GenericDensityPart<T>(n, m, detS, detT));
    } else if (task == "siegel_series") {
      int64_t p;
      int m;
      is >> p >> m;
      MyMatrix<T> S = ReadSquareMatrixEntries<T>(is, m);
      ListVal.push_back(
          PolynomialToGAP(SiegelSeriesPolynomial<T, Tint>(S, p, cache.siegel)));
    } else if (task == "rank_one_density") {
      int64_t p, a;
      int n;
      is >> p >> n;
      MyMatrix<T> Tg = ReadSquareMatrixEntries<T>(is, n);
      is >> a;
      SiegelJordanBlocks q = SiegelJordanBlocksOf<T, Tint>(Tg / T(2), p);
      push(RankOneDensity<T>(q, PadicFromInt64(a, p), cache.bad));
    } else if (task == "unit_chain") {
      int64_t p;
      int n, m;
      is >> p >> n;
      MyMatrix<T> Tg = ReadSquareMatrixEntries<T>(is, n);
      is >> m;
      MyMatrix<T> S = ReadSquareMatrixEntries<T>(is, m);
      std::optional<T> opt = UnitChainDensity<T, Tint>(
          S, OddDiagonalClasses<T, Tint>(Tg, p), p, cache.bad);
      if (opt) {
        push(*opt);
      } else {
        ListVal.push_back("fail");
      }
    } else if (task == "counting") {
      int64_t p;
      int n, m;
      is >> p >> n;
      MyMatrix<T> Tg = ReadSquareMatrixEntries<T>(is, n);
      is >> m;
      MyMatrix<T> S = ReadSquareMatrixEntries<T>(is, m);
      std::optional<T> opt = CountingDensity<T>(
          OddDiagonalClasses<T, Tint>(S, p), OddDiagonalClasses<T, Tint>(Tg, p),
          p, cache.bad);
      if (opt) {
        push(*opt);
      } else {
        ListVal.push_back("fail");
      }
    } else if (task == "reduction") {
      int64_t p;
      int n, m;
      is >> p >> n;
      MyMatrix<T> Tg = ReadSquareMatrixEntries<T>(is, n);
      is >> m;
      MyMatrix<T> S = ReadSquareMatrixEntries<T>(is, m);
      push(DensityOddDiagonal<T>(OddDiagonalClasses<T, Tint>(S, p),
                                 OddDiagonalClasses<T, Tint>(Tg, p), p,
                                 cache.bad));
    } else if (task == "mass_ratio") {
      int n, m;
      is >> n;
      MyMatrix<T> Tg = ReadSquareMatrixEntries<T>(is, n);
      is >> m;
      MyMatrix<T> S = ReadSquareMatrixEntries<T>(is, m);
      // the genus data is computed once per distinct T, as a sweep over the
      // root lattices of one genus would
      std::vector<T> key(Tg.data(), Tg.data() + n * n);
      auto iter = genus_data.find(key);
      if (iter == genus_data.end()) {
        iter = genus_data.emplace(key, ComputeGenusDensityData<T, Tint>(Tg))
                   .first;
      }
      std::optional<T> opt =
          MassRatioFromDensities<T, Tint>(S, iter->second, cache);
      if (opt) {
        push(*opt);
      } else {
        ListVal.push_back("fail");
      }
    } else if (task == "alpha_unimodular") {
      int64_t p;
      int k, xi_T, m;
      is >> p >> k >> xi_T >> m;
      MyMatrix<T> S = ReadSquareMatrixEntries<T>(is, m);
      push(AlphaUnimodularTarget<T, Tint>(S, p, k, xi_T, cache.siegel));
    } else {
      std::cerr << "GENUS_LocalDensity: unknown task " << task << "\n";
      throw TerminalException{1};
    }
    if (!is.good() && !is.eof()) {
      std::cerr << "GENUS_LocalDensity: malformed record " << i_rec
                << " in " << FileIn << "\n";
      throw TerminalException{1};
    }
  }
  os << "return [";
  for (size_t i = 0; i < ListVal.size(); i++) {
    os << (i > 0 ? ",\n" : "") << ListVal[i];
  }
  os << "];\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    std::string task;
    std::string FileIn;
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    add_rational_integer_arithmetics_by_library<false>(
        dispatch, [&]<typename T, typename Tint>() -> void {
          ProcessTask<T, Tint>(task, FileIn, *os_ptr);
        });
    if (argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GENUS_LocalDensity [arith] [task] [FileIn] [FileOut]\n";
      std::cerr << "\n";
      std::cerr << "arith values:\n";
      dispatch.print_descriptions(std::cerr);
      std::cerr << "\n";
      std::cerr << "task values (FileIn: a count, then that many records):\n";
      std::cerr << "  bernoulli              records: k\n";
      std::cerr << "  generalized_bernoulli  records: k D\n";
      std::cerr << "  euler_value            records: n m detS detT p\n";
      std::cerr << "  generic_part           records: n m detS detT\n";
      std::cerr << "  siegel_series          records: p m S (m*m entries)\n";
      std::cerr << "  alpha_unimodular       records: p k xi_T m S\n";
      std::cerr << "  rank_one_density       records: p n T a\n";
      std::cerr << "  unit_chain             records: p n T m S\n";
      std::cerr << "  counting               records: p n T m S\n";
      std::cerr << "  reduction              records: p n T m S\n";
      std::cerr << "  mass_ratio             records: n T m S\n";
      std::cerr << "FileOut: a GAP list, or stdout / stderr\n";
      return -1;
    }
    std::string arith = argv[1];
    task = argv[2];
    FileIn = argv[3];
    std::string FileOut = argv[4];
    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileOut, f);
    std::cerr << "Normal termination of GENUS_LocalDensity\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GENUS_LocalDensity\n";
    exit(e.eVal);
  }
  runtime(time);
}
