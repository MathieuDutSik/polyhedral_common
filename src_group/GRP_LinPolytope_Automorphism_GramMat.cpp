// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "GRP_GroupFct.h"
#include "Group.h"
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "Permutation.h"
#include "Temp_PolytopeEquiStab.h"
// clang-format on

template <typename T, typename Tgroup>
void full_process_A(std::string const &eFileEXT, std::string const &eFileGram,
                    std::string const &OutFormat, std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileEXT);
  MyMatrix<T> GramMat = ReadMatrixFile<T>(eFileGram);
  int nbCol = EXT.cols();
  int nbRow = EXT.rows();
  std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
  //
  Tgroup GRP = LinPolytope_Automorphism_GramMat<T, Tgroup>(
      EXT, GramMat, std::cerr);
  std::cerr << "|GRP|=" << GRP.size() << "\n";
  if (OutFormat == "GAP") {
    os << "return " << GRP.GapString() << ";\n";
    return;
  }
  if (OutFormat == "Oscar") {
    WriteGroup(os, GRP);
    return;
  }
  std::cerr << "GRP_LinPolytope_Automorphism : Failed to find matching entry "
               "for OutFormat="
            << OutFormat << "\n";
  throw TerminalException{1};
}

template <typename Tgroup>
void full_process_B(std::string const &arith, std::string const &eFileEXT,
                    std::string const &eFileGram, std::string const &OutFormat,
                    std::ostream &os) {
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
  using Trat = boost::multiprecision::mpq_rational;
#else
  using Trat = mpq_class;
#endif
  if (arith == "safe_rational") {
    using T = Rational<SafeInt64>;
    return full_process_A<T, Tgroup>(eFileEXT, eFileGram, OutFormat, os);
  }
  if (arith == "rational") {
    using T = Trat;
    return full_process_A<T, Tgroup>(eFileEXT, eFileGram, OutFormat, os);
  }
  if (arith == "Qsqrt5") {
    using T = QuadField<Trat, 5>;
    return full_process_A<T, Tgroup>(eFileEXT, eFileGram, OutFormat, os);
  }
  if (arith == "Qsqrt2") {
    using T = QuadField<Trat, 2>;
    return full_process_A<T, Tgroup>(eFileEXT, eFileGram, OutFormat, os);
  }
  std::optional<std::string> opt_realalgebraic =
      get_postfix(arith, "RealAlgebraic=");
  if (opt_realalgebraic) {
    std::string const &FileAlgebraicField = *opt_realalgebraic;
    if (!IsExistingFile(FileAlgebraicField)) {
      std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                << " is missing\n";
      throw TerminalException{1};
    }
    HelperClassRealField<Trat> hcrf(FileAlgebraicField);
    int const idx_real_algebraic_field = 1;
    insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
    using T = RealField<idx_real_algebraic_field>;
    return full_process_A<T, Tgroup>(eFileEXT, eFileGram, OutFormat, os);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_LinPolytope_Automorphism_GramMat Arith [FileEXT] "
                   "[FileGram] [OutFormat] [OutGroup]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_LinPolytope_Automorphism_GramMat Arith [FileEXT] "
                   "[FileGram]\n";
      std::cerr << "\n";
      std::cerr << "FileEXT   : The list of vertices\n";
      std::cerr << "FileGram  : The Gram matrix\n";
      std::cerr << "OutFormat : The output format (GAP or Oscar)\n";
      std::cerr << "OutGroup  : The automorphism group file\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr
          << "safe_rational          : rational based on int64_t that fails\n";
      std::cerr << "    gracefully in overflow\n";
      std::cerr << "rational               : rational arithmetic on input\n";
      std::cerr
          << "Qsqrt2                 : arithmetic over the field Q(sqrt(2))\n";
      std::cerr
          << "Qsqrt5                 : arithmetic over the field Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc : For the real algebraic case of a ";
      std::cerr << "    field whose description is in FileDesc\n";
      return -1;
    }
    //
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Tint = boost::multiprecision::mpz_int;
#else
    using Tint = mpz_class;
#endif
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string arith = argv[1];
    std::string eFileEXT = argv[2];
    std::string eFileGram = argv[3];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      FileOut = argv[5];
    }
    //
    if (FileOut == "stderr") {
      full_process_B<Tgroup>(arith, eFileEXT, eFileGram, OutFormat, std::cerr);
    } else {
      if (FileOut == "stdout") {
        full_process_B<Tgroup>(arith, eFileEXT, eFileGram, OutFormat,
                               std::cout);
      } else {
        std::ofstream os(FileOut);
        full_process_B<Tgroup>(arith, eFileEXT, eFileGram, OutFormat, os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
