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
#include "QuadField.h"
#include "ShortestUniversal.h"
#include "Shvec_exact.h"
// clang-format on


template<typename T, typename Tint>
void process(std::string const& choice, std::string const& FileGram, std::string const& FileVect, std::string const& OutFormat, std::ostream &os) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(FileGram);
  MyVector<T> eV = ReadVectorFile<T>(FileVect);
  int n = GramMat.rows();
  auto get_result=[&]() -> MyMatrix<Tint> {
    if (choice == "nearest") {
      resultCVP<T, Tint> res = CVPVallentinProgram_exact<T,Tint>(GramMat, eV);
      return res.ListVect;
    }
    std::optional<std::string> opt_near = get_postfix(choice, "near=");
    if (opt_near) {
      T norm = ParseScalar<T>(*opt_near);
      MyVector<T> fV = -eV;
      std::vector<MyVector<Tint>> ListVect = FindAtMostNormVectors<T,Tint>(GramMat, fV, norm);
      if (ListVect.size() == 0)
        return ZeroMatrix<Tint>(0,n);
      return MatrixFromVectorFamily(ListVect);
    }
    std::cerr << "Failed to find a matching entry for choice=" << choice << "\n";
    throw TerminalException{1};
  };
  MyMatrix<Tint> result = get_result();
  auto write_result=[&]() -> void {
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
    std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
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
      std::cerr << "POLY_cdd_lp2 arith choice [FileGram] [FileV] [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_cdd_lp2 arith choice [FileGram] [FileV]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith     : the chosen arithmetic (see below)\n";
      std::cerr << "choice    : The choose option (see below)\n";
      std::cerr << "FileGram  : The list of inequalities\n";
      std::cerr << "FileV     : The vector for which we want the neainequality to be minimized\n";
      std::cerr << "OutFormat : The format of output, GAP or Oscar\n";
      std::cerr << "FileOut   : The file of output (if present, otherwise std::cerr)\n";
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
      std::cerr << "near=dist : The vector up to distance near from that vector\n";
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
        return process<T,Tint>(choice, FileGram, FileVect, OutFormat, os);
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
    if (FileOut == "stderr") {
      call_SV(std::cerr);
    } else {
      if (FileOut == "stdout") {
        call_SV(std::cout);
      } else {
        std::ofstream os(FileOut);
        call_SV(os);
      }
    }
    std::cerr << "Normal termination of the program sv_near\n";
  } catch (TerminalException const &e) {
    std::cerr << "Raised exception led to premature end of sv_near\n";
    exit(e.eVal);
  }
  runtime(time);
}
