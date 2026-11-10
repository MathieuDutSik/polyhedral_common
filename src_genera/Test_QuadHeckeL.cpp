// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "QuadraticHeckeL.h"
#include <string>
#include <vector>
// clang-format on

/*
  The quadratic Hecke character chi_{E/K} of a quadratic lattice over a
  totally real abelian field: the field K, the discriminant d, |D_K|, the
  norm of the conductor and the Dirichlet coefficients a_1, ..., a_X of
  L(s, chi), and L(1 - r, chi) for m = 2r, for a comparison with PARI.

  Test_QuadHeckeL [Gram] [Omega] [n] [h1,h2,...] [X] [prec_bits]
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
void Process(std::string const &FileGram, std::string const &FileOmega,
             int64_t n, std::vector<int64_t> const &gens, int64_t X,
             int prec_bits) {
  MyMatrix<T> B = ReadMatrixFile<T>(FileGram);
  std::vector<MyMatrix<T>> Omega = ReadListMatrixFile<T>(FileOmega);
  int m = B.rows() / static_cast<int>(Omega.size());
  AbelianCMField FK = MakeAbelianRealField(n, gens);
  std::optional<std::pair<CyclotomicField<T>, ElemE<T>>> Kd =
      QuadDiscriminantK(B, Omega);
  if (!Kd) {
    std::cout << "fail\n";
    return;
  }
  auto f_none = [](int64_t, std::vector<QuadHeckePrime> const &) -> void {};
  std::optional<QuadHeckeCharacter<T>> C = ComputeQuadHeckeCharacter(
      B, Omega, FK, Kd->second, Kd->first, m, X, true, f_none);
  if (!C) {
    std::cout << "fail\n";
    return;
  }
  std::cout << "f:";
  for (auto const &c : Kd->first.Phi) {
    std::cout << " " << c;
  }
  std::cout << "\nd:";
  for (auto const &c : Kd->second) {
    std::cout << " " << c;
  }
  std::cout << "\nDK: " << C->D_K << "\nNf: " << C->norm_conductor
            << "\ndelta: " << C->delta << "\na:";
  for (auto const &x : QuadHeckeCoefficients(*C, X)) {
    std::cout << " " << x;
  }
  std::cout << "\n";
  if (prec_bits > 0) {
    std::optional<QuadHeckeLData<T>> Ld =
        MakeQuadHeckeLData(B, Omega, FK, Kd->second, Kd->first, m);
    std::optional<mpfr::mpreal> L = QuadHeckeLNegative(*Ld, prec_bits, false);
    if (L) {
      std::cout << "L: "
                << L->toString(static_cast<int>(prec_bits * 0.30103) + 5)
                << "\n";
    } else {
      std::cout << "L: fail\n";
    }
#ifdef ENABLE_FLINT_SUPPORT
    std::optional<ArbReal> La = QuadHeckeLNegativeArb(*Ld, prec_bits);
    if (La) {
      std::cout << "Larb: "
                << La->str(static_cast<slong>(prec_bits * 0.30103) + 5) << "\n";
    } else {
      std::cout << "Larb: fail\n";
    }
#endif
  }
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 7) {
      std::cerr
          << "Test_QuadHeckeL [Gram] [Omega] [n] [h1,h2,...] [X] [prec_bits]\n";
      return -1;
    }
    Process<mpq_class>(argv[1], argv[2], std::stoll(argv[3]),
                       ParseList(argv[4]), std::stoll(argv[5]),
                       std::stoi(argv[6]));
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
