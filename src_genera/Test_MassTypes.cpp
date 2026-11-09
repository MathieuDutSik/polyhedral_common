// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "NumberTheoryBoostGmpInt.h"
#include "AbelianHermitianMass.h"
#include "AbelianQuadraticMass.h"
#include <string>
#include <vector>
// clang-format on

/*
  The masses of AbelianHermitianMass.h and AbelianQuadraticMass.h computed
  with several rational types: they must agree with the mpq_class values.

  Test_MassTypes herm [Gram] [Omega] [n] [h1,h2,...]
  Test_MassTypes quad [Gram] [Omega] [n] [h1,h2,...]
 */

std::vector<int64_t> ParseList(std::string const &lst) {
  std::vector<int64_t> res;
  std::stringstream ss(lst);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    if (!tok.empty()) {
      res.push_back(std::stoll(tok));
    }
  }
  return res;
}

template <typename T>
std::string ComputeMass(std::string const &mode, std::string const &FileGram,
                        std::string const &FileOmega, int64_t n,
                        std::vector<int64_t> const &gens) {
  MyMatrix<T> B = ReadMatrixFile<T>(FileGram);
  std::vector<MyMatrix<T>> Omega = ReadListMatrixFile<T>(FileOmega);
  std::optional<T> mass;
  if (mode == "herm") {
    mass = HermitianMassAbelianLattice(B, Omega, MakeAbelianCMField(n, gens));
  } else {
    mass = QuadraticMassAbelianLattice(B, Omega, MakeAbelianRealField(n, gens));
  }
  if (!mass) {
    return "fail";
  }
  std::stringstream ss;
  ss << *mass;
  return ss.str();
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 6) {
      std::cerr << "Test_MassTypes [herm|quad] [Gram] [Omega] [n] [h1,h2,...]\n";
      return -1;
    }
    std::string mode = argv[1];
    int64_t n = std::stoll(argv[4]);
    std::vector<int64_t> gens = ParseList(argv[5]);
    std::vector<std::pair<std::string, std::string>> res;
    res.push_back({"mpq_class", ComputeMass<mpq_class>(mode, argv[2], argv[3], n, gens)});
    res.push_back({"mpq_rational",
                   ComputeMass<boost::multiprecision::mpq_rational>(mode, argv[2], argv[3], n, gens)});
#ifdef ENABLE_FLINT_SUPPORT
    res.push_back({"fmpq_class", ComputeMass<fmpq_class>(mode, argv[2], argv[3], n, gens)});
#endif
    bool ok = true;
    for (auto const &r : res) {
      std::cout << r.first << " " << r.second << "\n";
      if (r.second != res[0].second) {
        ok = false;
      }
    }
    std::cout << (ok ? "agree" : "DISAGREE") << "\n";
    return ok ? 0 : 1;
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
