// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "LocalDensity.h"
#include "Arithmetic_dispatch_types.h"
#include <fstream>
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
 */

template <typename T>
void ProcessTask(std::string const &task, std::string const &FileIn,
                 std::ostream &os) {
  std::ifstream is(FileIn);
  if (!is.good()) {
    std::cerr << "GENUS_LocalDensity: cannot read " << FileIn << "\n";
    throw TerminalException{1};
  }
  size_t n_rec;
  is >> n_rec;
  std::vector<T> ListVal;
  for (size_t i_rec = 0; i_rec < n_rec; i_rec++) {
    if (task == "bernoulli") {
      int k;
      is >> k;
      ListVal.push_back(BernoulliNumbers<T>(k)[k]);
    } else if (task == "generalized_bernoulli") {
      int k;
      int64_t D;
      is >> k >> D;
      ListVal.push_back(GeneralizedBernoulli<T>(k, D));
    } else if (task == "euler_value") {
      int n, m;
      int64_t detS, detT, p;
      is >> n >> m >> detS >> detT >> p;
      ListVal.push_back(
          EulerValueAtPrime<T>(GoodPrimeEulerTerms(n, m, detS, detT), p));
    } else if (task == "generic_part") {
      int n, m;
      int64_t detS, detT;
      is >> n >> m >> detS >> detT;
      ListVal.push_back(GenericDensityPart<T>(n, m, detS, detT));
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
    add_rational_arithmetics_by_library(dispatch, [&]<typename T>() -> void {
      ProcessTask<T>(task, FileIn, *os_ptr);
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
