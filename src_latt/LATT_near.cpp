/* sv.c  simple driver for shvec                              */
/* Version July 11, 2005                                      */
/* Copyright: Frank Vallentin 2005, frank.vallentin@gmail.com */

// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "ShortestUniversal.h"
#include "Shvec_exact.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &choice, std::string const &FileGram,
             std::string const &FileVect, std::string const &OutFormat,
             std::ostream &os) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileGram);
  int n = GramMat.rows();
  auto get_vect = [&]() -> MyVector<T> {
    if (FileVect == "zero") {
      return ZeroVector<T>(n);
    }
    return ReadVectorFile<T>(FileVect);
  };
  MyVector<T> eV = get_vect();
  auto get_result = [&]() -> MyMatrix<Tint> {
    if (choice == "nearest") {
      resultCVP<T, Tint> res =
          CVPVallentinProgram_exact<T, Tint>(GramMat, eV, std::cerr);
      return res.ListVect;
    }
    std::optional<std::string> opt_near = get_postfix(choice, "near=");
    if (opt_near) {
      T norm = ParseScalar<T>(*opt_near);
      MyVector<T> fV = -eV;
      std::vector<MyVector<Tint>> ListVect =
        FindAtMostNormVectors<T, Tint>(GramMat, fV, norm, std::cerr);
      return MatrixFromVectorFamilyDim(n, ListVect);
    }
    std::optional<std::string> opt_fixed = get_postfix(choice, "fixed=");
    if (opt_fixed) {
      T norm = ParseScalar<T>(*opt_fixed);
      MyVector<T> fV = -eV;
      std::vector<MyVector<Tint>> ListVect =
        FindFixedNormVectors<T, Tint>(GramMat, fV, norm, std::cerr);
      return MatrixFromVectorFamilyDim(n, ListVect);
    }
    std::cerr << "Failed to find a matching entry for choice=" << choice
              << "\n";
    throw TerminalException{1};
  };
  MyMatrix<Tint> result = get_result();
  auto write_result = [&]() -> void {
    if (OutFormat == "GAP") {
      os << "return ";
      WriteMatrixGAP(os, result);
      os << ";\n";
      return;
    }
    if (OutFormat == "Oscar") {
      WriteMatrix(os, result);
      return;
    }
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
              << "\n";
    throw TerminalException{1};
  };
  write_result();
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Trat = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
#else
    using Trat = mpq_class;
    using Tint = mpz_class;
#endif
    if (argc != 5 && argc != 7) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_near arith choice [FileGram] [FileV] [OutFormat] "
                   "[FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "LATT_near arith choice [FileGram] [FileV]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith     : the chosen arithmetic (see below)\n";
      std::cerr << "choice    : The choose option (see below)\n";
      std::cerr << "FileGram  : The list of inequalities\n";
      std::cerr << "FileV     : The vector for which we want to\n";
      std::cerr << "       compute the distance\n";
      std::cerr << "OutFormat : The format of output, GAP or Oscar\n";
      std::cerr << "FileOut   : The file of output (if present, otherwise "
                   "std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "rational  : rational arithmetic on input\n";
      std::cerr << "Qsqrt2    : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5    : arithmetic over the field Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc  : For the real algebraic case of a";
      std::cerr << "  field whose description is in FileDesc\n";
      std::cerr << "\n";
      std::cerr << "        --- choice ---\n";
      std::cerr << "\n";
      std::cerr << "nearest   : The nearest points to that vector\n";
      std::cerr << "near=dist : The vector up to distance near from\n";
      std::cerr << "     that vector\n";
      std::cerr << "fixed=dist : The vector at an exact distance\n";
      std::cerr << "     from that vector\n";

      return -1;
    }
    //
    std::string arith = argv[1];
    std::string choice = argv[2];
    std::string FileGram = argv[3];
    std::string FileVect = argv[4];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      FileOut = argv[6];
    }

    auto call_SV = [&](std::ostream &os) -> void {
      if (arith == "rational") {
        using T = Trat;
        return process<T, Tint>(choice, FileGram, FileVect, OutFormat, os);
      }
      /*
      if (arith == "Qsqrt5") {
        using T = QuadField<Trat, 5>;
        return process<T,T>(choice, FileGram, FileVect, OutFormat, os);
      }
      if (arith == "Qsqrt2") {
        using T = QuadField<Trat, 2>;
        return process<T,T>(choice, FileGram, FileVect, OutFormat, os);
      }
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        std::string FileAlgebraicField = *opt_realalgebraic;
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<Trat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return process<T,T>(choice, FileGram, FileVect, OutFormat, os);
      }
      */
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileOut, call_SV);
    std::cerr << "Normal termination of LATT_near\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_near\n";
    exit(e.eVal);
  }
  runtime(time);
}
