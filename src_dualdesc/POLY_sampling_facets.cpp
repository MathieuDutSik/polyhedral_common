// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "POLY_SamplingFacet.h"
#include "POLY_RecursiveDualDesc.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::string const &ansSamp,
             std::string const &OutFormat, std::ostream &os_out) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  std::vector<int> eList = ColumnReductionSet(EXT);
  MyMatrix<T> EXTred = SelectColumn(EXT, eList);
  vectface vf = DirectComputationInitialFacetSet(EXTred, ansSamp, std::cerr);
  OutputFacets_stream(EXT, vf, os_out, OutFormat, os_out);
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat
            << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr
          << "POLY_sampling_facets arith command [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_sampling_facets arith command [FileI]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith     : the chosen arithmetic\n";
      std::cerr << "command   : the command used for sampling the facets\n";
      std::cerr
          << "FileI     : The polyhedral cone with the extreme rays in EXT\n";
      std::cerr << "OutFormat : The format of output, GAP or Oscar\n";
      std::cerr << "FileO     : The file of output (if present, otherwise "
                   "std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "safe_rational : rational arithmetic based on int64_t that "
                   "fails\n";
      std::cerr << "    gracefully on overflow\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr
          << "RealAlgebraic=FileDesc  : For the real algebraic case of a\n";
      std::cerr << "    field whose description is in FileDesc\n";
      std::cerr << "\n";
      std::cerr << "        --- command ---\n";
      std::cerr << "\n";
      std::cerr << "lp_cdd                   : 10 facets being sampled by "
                   "using linear programming\n";
      std::cerr << "lp_cdd:iter_X            : X facets being sampled by using "
                   "linear programming\n";
      std::cerr << "lp_cdd_min               : 10 facets (or less) of minimum "
                   "incidence being sampled\n";
      std::cerr << "lp_cdd_min:iter_X        : X facets (or less) of minimum "
                   "incidence being sampled\n";
      std::cerr
          << "sampling                 : The recursive sampling algorithm\n";
      std::cerr << "lrs_limited              : lrs limited to the first 100 "
                   "vertices being found\n";
      std::cerr << "lrs_limited:upperlimit_X : lrs limited to the first X "
                   "vertices being found\n";
      std::cerr << "\n";
      std::cerr << "        --- OutFormat ---\n";
      std::cerr << "\n";
      std::cerr << "GAP                   : The incidence in a file readable "
                   "in GAP via ReadAsFunction\n";
      std::cerr << "Oscar                 : The incidence in a 0/1 matrix "
                   "representing the obtained rays\n";
      std::cerr << "FacetInequalities     : The inequalities in a file. If EXT "
                   "is not full dimensional, then it is underdefined\n";
      return -1;
    }
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Trat = boost::multiprecision::mpq_rational;
#else
    using Trat = mpq_class;
#endif
    //
    std::string arith = argv[1];
    std::string command = argv[2];
    std::string eFileI = argv[3];
    std::string OutFormat = "GAP";
    std::string eFileO = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      eFileO = argv[5];
    }
    auto call_lrs = [&](std::ostream &os) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(eFileI, command, OutFormat, os);
      }
      if (arith == "rational") {
        using T = Trat;
        return process<T>(eFileI, command, OutFormat, os);
      }
      if (arith == "Qsqrt5") {
        using T = QuadField<Trat, 5>;
        return process<T>(eFileI, command, OutFormat, os);
      }
      if (arith == "Qsqrt2") {
        using T = QuadField<Trat, 2>;
        return process<T>(eFileI, command, OutFormat, os);
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
        return process<T>(eFileI, command, OutFormat, os);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(eFileO, call_lrs);
    std::cerr << "Normal termination of POLY_sampling_facets\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_sampling_facets\n";
    exit(e.eVal);
  }
  runtime(time1);
}
