// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryGmp.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "POLY_RecursiveDualDesc.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

template <typename T, typename Tidx>
void Process_EXT(MyMatrix<T> const& EXT, std::string const& GRPfile, std::string const& OutFormat, std::string const& OutFile) {
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  Tgroup GRP = ReadGroupFile<Tgroup>(GRPfile);
  MyMatrix<T> EXTred = ColumnReduction(EXT);
  int dimEXT = EXTred.cols();
  PolyHeuristicSerial<TintGroup> AllArr =
    AllStandardHeuristicSerial<T, TintGroup>(dimEXT, std::cerr);
  vectface TheOutput = DualDescriptionStandard(EXTred, GRP, AllArr, std::cerr);
  OutputFacets_file(EXTred, GRP, TheOutput, OutFile,
                    OutFormat, std::cerr);
}

template <typename T>
void Process(std::string const& EXTfile, std::string const& GRPfile, std::string const& OutFormat, std::string const& OutFile) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(EXTfile);
  //
  if (size_t(EXT.rows()) < std::numeric_limits<uint8_t>::max())
    return Process_EXT<T, uint8_t>(EXT, GRPfile, OutFormat, OutFile);
  if (size_t(EXT.rows()) < std::numeric_limits<uint16_t>::max())
    return Process_EXT<T, uint16_t>(EXT, GRPfile, OutFormat, OutFile);
  if (size_t(EXT.rows()) < std::numeric_limits<uint32_t>::max())
    return Process_EXT<T, uint32_t>(EXT, GRPfile, OutFormat, OutFile);
#if !defined __APPLE__
  if (size_t(EXT.rows()) < std::numeric_limits<uint64_t>::max())
    return Process_EXT<T, uint64_t>(EXT, GRPfile, OutFormat, OutFile);
#endif
  std::cerr << "Failed to find a numeric type that matches\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_DirectSerialDualDesc [arith] [EXTfile] [GRPfile] [OutFormat] [OutFile]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string EXTfile = argv[2];
    std::string GRPfile = argv[3];
    std::string OutFormat = argv[4];
    std::string OutFile = argv[5];
    //
    auto f = [&]() -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return Process<T>(EXTfile, GRPfile, OutFormat, OutFile);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return Process<T>(EXTfile, GRPfile, OutFormat, OutFile);
      }
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return Process<T>(EXTfile, GRPfile, OutFormat, OutFile);
      }
      if (arith == "mpq_rational") {
        using T = boost::multiprecision::mpq_rational;
        return Process<T>(EXTfile, GRPfile, OutFormat, OutFile);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return Process<T>(EXTfile, GRPfile, OutFormat, OutFile);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return Process<T>(EXTfile, GRPfile, OutFormat, OutFile);
      }
      std::cerr << "Failed to find a matching type entry arith=" << arith << "\n";
      throw TerminalException{1};
    };
    f();
    //
    std::cerr << "Normal termination of POLY_DirectSerialDualDesc\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_DirectSerialDualDesc\n";
    exit(e.eVal);
  }
  runtime(time1);
}
