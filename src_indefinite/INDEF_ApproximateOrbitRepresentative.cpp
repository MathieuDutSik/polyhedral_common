// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "ApproximateModels.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &MatFile, std::string const &XnormStr,
             std::string const &OutFormat, std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  // The group order is counted in TintGroup, which is unrelated to the
  // arithmetic of the coordinates: fixing it to mpz_class keeps the
  // permutation group out of whichever arithmetic Tint happens to be.
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  MyMatrix<T> Qmat = ReadMatrixFile<T>(MatFile);
  T Xnorm = ParseScalar<T>(XnormStr);
  ApproximateModel<T, Tint> approx =
      INDEF_FORM_EichlerCriterion_TwoHyperplanesEven<T, Tint, Tgroup>(Qmat);
  std::vector<MyVector<Tint>> LVect =
      approx.GetCoveringOrbitRepresentatives(Xnorm, std::cerr);
  if (OutFormat == "GAP") {
    if (LVect.empty()) {
      os_out << "return rec(LVect:=[]);\n";
    } else {
      MyMatrix<Tint> MatVect = MatrixFromVectorFamily(LVect);
      os_out << "return rec(LVect:=" << StringMatrixGAP(MatVect) << ");\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr
          << "INDEF_ApproximateOrbitRepresentative [arith] [MatFile] [X]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_ApproximateOrbitRepresentative [arith] [MatFile] [X] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "        --------\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Allowed values for arith: gmp, flint\n";
#else
      std::cerr << "Allowed values for arith: gmp (build with "
                << "ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string MatFile = argv[2];
    std::string XnormStr = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(MatFile, XnormStr, OutFormat, os);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return process<T, Tint>(MatFile, XnormStr, OutFormat, os);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return process<T, Tint>(MatFile, XnormStr, OutFormat, os);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        using Tint = fmpz_class;
        return process<T, Tint>(MatFile, XnormStr, OutFormat, os);
      }
#endif
      std::cerr << "Failed to find matching type for arith\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Allowed values: gmp, flint (and, with "
                << "ENABLE_BOOST_TYPES, gmp_boost and "
                << "multi_boost)\n";
#else
      std::cerr << "Allowed values: gmp (build with "
                << "ENABLE_FLINT_SUPPORT=1 for flint; with "
                << "ENABLE_BOOST_TYPES for the boost ones)\n";
#endif
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(OutFile, f);
    //
    std::cerr << "Normal termination of INDEF_ApproximateOrbitRepresentative\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in INDEF_ApproximateOrbitRepresentative runtime="
              << time << "\n";
    exit(e.eVal);
  }
  runtime(time);
}
