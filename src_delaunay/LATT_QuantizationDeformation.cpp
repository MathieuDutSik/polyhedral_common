// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "QuantizationDeformation.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

// Given a positive definite form Q and a general symmetric perturbation H,
// study the deformation Q + t H: compute the Taylor data at t = 0 of
// SecMoment(t) and the derivatives of the normalized quantizer G(t).
template <typename T, typename Tint, typename Tgroup>
void process(std::string const &Qfile, std::string const &Hfile,
             std::ostream &os) {
  MyMatrix<T> Q = ReadMatrixFile<T>(Qfile);
  MyMatrix<T> H = ReadMatrixFile<T>(Hfile);
  int n = Q.rows();
  os << "QDEF: n=" << n << "\n";
  //
  DeformationDerivatives<T> der =
      compute_deformation_derivatives<T, Tint, Tgroup>(Q, H, os);
  os << "QDEF: SecMoment(t) reconstructed as a rational function of degree "
     << der.secmoment_degree << "\n";
  os << "QDEF: SecMoment numerator   = " << StringVectorGAP(der.secmoment_num)
     << "\n";
  os << "QDEF: SecMoment denominator = " << StringVectorGAP(der.secmoment_den)
     << "\n";
  os << "QDEF: SecMoment(0)   = " << der.S0 << "\n";
  os << "QDEF: SecMoment'(0)  = " << der.S1 << "\n";
  os << "QDEF: SecMoment''(0) = " << der.S2 << "\n";
  os << "QDEF: det(Q)         = " << der.det0 << "\n";
  os << "QDEF: det'(0)        = " << der.det1 << "\n";
  os << "QDEF: det''(0)       = " << der.det2 << "\n";
  os << "QDEF: G(0)           = " << der.G0 << "\n";
  os << "QDEF: G'(0)          = " << der.G1 << "\n";
  os << "QDEF: G''(0)         = " << der.G2 << "\n";
}

template <typename T, typename Tint>
void process_B(std::string const &Qfile, std::string const &Hfile) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Tint_grp = mpz_class;
  using Tgroup = permutalib::Group<Telt, Tint_grp>;
  return process<T, Tint, Tgroup>(Qfile, Hfile, std::cerr);
}

void process_C(std::string const &arith, std::string const &Qfile,
               std::string const &Hfile) {
  if (arith == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return process_B<T, Tint>(Qfile, Hfile);
  }
  if (arith == "safe") {
    using T = Rational<SafeInt64>;
    using Tint = SafeInt64;
    return process_B<T, Tint>(Qfile, Hfile);
  }
  std::cerr << "LATT_QuantizationDeformation: no match for arith=" << arith
            << "\n";
  std::cerr << "Available types: gmp, safe\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_QuantizationDeformation [arith] [Qfile] [Hfile]\n";
      std::cerr << "     arith   gmp, safe\n";
      std::cerr << "     Qfile   the positive definite form Q\n";
      std::cerr << "     Hfile   the symmetric perturbation H (any symmetric "
                   "form)\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string Qfile = argv[2];
    std::string Hfile = argv[3];
    process_C(arith, Qfile, Hfile);
    std::cerr << "Normal termination of LATT_QuantizationDeformation\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_QuantizationDeformation\n";
    exit(e.eVal);
  }
  runtime(time);
}
