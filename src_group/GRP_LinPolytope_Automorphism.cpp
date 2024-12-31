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
#include "PolytopeEquiStab.h"
// clang-format on

template <typename T, typename Tgroup>
void full_process_A(std::string const &eFile, std::string const &OutFormat,
                    std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFile);
  //
  Tgroup GRP = LinPolytope_Automorphism<T, Tgroup>(EXT, std::cerr);
  std::cerr << "|GRP|=" << GRP.size() << "\n";
  if (OutFormat == "GAP") {
    os << "return " << GRP.GapString() << ";\n";
    return;
  }
  if (OutFormat == "ListMatrixFile") {
    int rnk = RankMat(EXT);
    if (rnk != EXT.cols()) {
      std::cerr << "We have rnk=" << rnk << " but nbCol=" << EXT.cols() << "\n";
      std::cerr << "Therefore, we cannot build the matrix generators\n";
      throw TerminalException{1};
    }
    std::vector<MyMatrix<T>> ListGenMatr;
    for (auto &eGen : GRP.SmallGeneratingSet()) {
      MyMatrix<T> eGenMatr = RepresentVertexPermutation(EXT, EXT, eGen);
      ListGenMatr.push_back(eGenMatr);
    }
    WriteListMatrix(os, ListGenMatr);
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
void full_process_B(std::string const &arith, std::string const &eFile,
                    std::string const &OutFormat, std::ostream &os) {
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
  using Trat = boost::multiprecision::mpq_rational;
#else
  using Trat = mpq_class;
#endif
  if (arith == "safe_rational") {
    using T = Rational<SafeInt64>;
    return full_process_A<T, Tgroup>(eFile, OutFormat, os);
  }
  if (arith == "rational") {
    using T = Trat;
    return full_process_A<T, Tgroup>(eFile, OutFormat, os);
  }
  if (arith == "Qsqrt5") {
    using T = QuadField<Trat, 5>;
    return full_process_A<T, Tgroup>(eFile, OutFormat, os);
  }
  if (arith == "Qsqrt2") {
    using T = QuadField<Trat, 2>;
    return full_process_A<T, Tgroup>(eFile, OutFormat, os);
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
    return full_process_A<T, Tgroup>(eFile, OutFormat, os);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_LinPolytope_Automorphism arith [FileEXT] [OutFormat] "
                   "[FileGRP]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Automorphism arith [FileEXT]\n";
      std::cerr << "\n";
      std::cerr << "FileEXT   : The list of vectors\n";
      std::cerr << "OutFormat : The format of output (GAP, ListMatrixFile, or "
                   "Oscar)\n";
      std::cerr << "FileGRP   : The file for outputting the group\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr
          << "safe_rational           : safe_rational arithmetic on uint64_t\n";
      std::cerr << "     that fails gracefully\n";
      std::cerr << "rational                : rational arithmetic on input\n";
      std::cerr << "Qsqrt2                  : arithmetic over the field\n";
      std::cerr << "                          Q(sqrt(2))\n";
      std::cerr << "Qsqrt5                  : arithmetic over the field\n";
      std::cerr << "                          Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc  : For the real algebraic case\n";
      std::cerr << "          of a field whose description is in FileDesc\n";
      std::cerr << "\n";
      std::cerr << "        --- OutFormat ---\n";
      std::cerr << "\n";
      std::cerr
          << "GAP            : Output to a GAP readable format (default)\n";
      std::cerr
          << "ListMatrixFile : To be used by other software of polyhedral\n";
      std::cerr << "Oscar          : Output to a CPP like format\n";
      std::cerr << "\n";
      std::cerr << "        --- FileGRP ---\n";
      std::cerr << "\n";
      std::cerr << "stderr : standard error (default)\n";
      std::cerr << "stdout : standard output\n";
      std::cerr << "otherwise output to the named file\n";
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
    std::string FileEXT = argv[2];
    std::string OutFormat = "GAP";
    std::string FileGRP = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileGRP = argv[4];
    }
    //
    if (FileGRP == "stderr") {
      full_process_B<Tgroup>(arith, FileEXT, OutFormat, std::cerr);
    } else {
      if (FileGRP == "stdout") {
        full_process_B<Tgroup>(arith, FileEXT, OutFormat, std::cout);
      } else {
        std::ofstream os(FileGRP);
        full_process_B<Tgroup>(arith, FileEXT, OutFormat, os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
