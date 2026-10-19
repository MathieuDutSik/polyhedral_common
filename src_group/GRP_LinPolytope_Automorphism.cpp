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
#include "MatrixGroup.h"
#include "Arithmetic_dispatch_types.h"
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
  if (OutFormat == "RecGAP") {
    std::string strGAPmatr = get_matrs_as_string(EXT, GRP.GeneratorsOfGroup());
    os << "return rec(GAPperm:=" << GRP.GapString()
       << ", GAPmatr:=" << strGAPmatr << ");";
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
  if (OutFormat == "Oscar" || OutFormat == "CPP") {
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
    std::string FileEXT;
    std::string OutFormat = "GAP";
    std::string FileGRP = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      full_process_A<T, Tgroup>(FileEXT, OutFormat, *os_ptr);
    };
    add_rational_arithmetics<Trat>(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
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
                   "Oscar/CPP)\n";
      std::cerr << "FileGRP   : The file for outputting the group\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
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
    std::string arith = argv[1];
    FileEXT = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      FileGRP = argv[4];
    }
    //
    auto f = [&](std::ostream &os_out) -> void {
      os_ptr = &os_out;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileGRP, f);
    std::cerr << "Normal termination of GRP_LinPolytope_Automorphism\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Automorphism\n";
    exit(e.eVal);
  }
  runtime(time1);
}
