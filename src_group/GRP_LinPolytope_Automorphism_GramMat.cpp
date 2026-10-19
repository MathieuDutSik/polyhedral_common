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
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Permutation.h"
#include "PolytopeEquiStab.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T, typename Tgroup>
void full_process_A(std::string const &eFileEXT, std::string const &eFileGram,
                    std::string const &OutFormat, std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileEXT);
  MyMatrix<T> GramMat = ReadMatrixFile<T>(eFileGram);
  //
  Tgroup GRP =
      LinPolytope_Automorphism_GramMat<T, Tgroup>(EXT, GramMat, std::cerr);
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

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Tint = boost::multiprecision::mpz_int;
    using Trat = boost::multiprecision::mpq_rational;
#else
    using Tint = mpz_class;
    using Trat = mpq_class;
#endif
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    std::string eFileEXT;
    std::string eFileGram;
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      full_process_A<T, Tgroup>(eFileEXT, eFileGram, OutFormat, *os_ptr);
    };
    add_rational_arithmetics<Trat>(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Automorphism_GramMat Arith [FileEXT] "
                   "[FileGram] [OutFormat] [OutGroup]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Automorphism_GramMat Arith [FileEXT] "
                   "[FileGram]\n";
      std::cerr << "\n";
      std::cerr << "FileEXT   : The list of vertices\n";
      std::cerr << "FileGram  : The Gram matrix\n";
      std::cerr << "OutFormat : The output format (GAP or Oscar)\n";
      std::cerr << "OutGroup  : The automorphism group file\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    eFileEXT = argv[2];
    eFileGram = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      FileOut = argv[5];
    }
    //
    auto f = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileOut, f);
    std::cerr << "Normal termination of GRP_LinPolytope_Automorphism_GramMat\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Automorphism_GramMat\n";
    exit(e.eVal);
  }
  runtime(time1);
}
