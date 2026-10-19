// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
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
#include "POLY_SamplingFacet.h"
#include "POLY_RecursiveDualDesc.h"
#include "Arithmetic_dispatch_types.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::string const &ansSamp,
             std::string const &OutFormat, std::ostream &os_out) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  std::vector<int> eList = ColumnReductionSet(EXT);
  MyMatrix<T> EXTred = SelectColumn(EXT, eList);
  vectface vf = DirectComputationInitialFacetSet(EXTred, ansSamp, std::cerr);
  OutputFacets_stream(EXT, vf, os_out, OutFormat, os_out);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
    using Trat = boost::multiprecision::mpq_rational;
#else
    using Trat = mpq_class;
#endif
    //
    std::string command;
    std::string eFileI;
    std::string OutFormat = "GAP";
    std::string eFileO = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(eFileI, command, OutFormat, *os_ptr);
    };
    add_rational_arithmetics<Trat>(dispatch, f_arith);
    add_algebraic_field_arithmetics(dispatch, f_arith);
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
      dispatch.print_descriptions(std::cerr);
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
    std::string arith = argv[1];
    command = argv[2];
    eFileI = argv[3];
    if (argc == 6) {
      OutFormat = argv[4];
      eFileO = argv[5];
    }
    auto call_lrs = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(eFileO, call_lrs);
    std::cerr << "Normal termination of POLY_sampling_facets\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_sampling_facets\n";
    exit(e.eVal);
  }
  runtime(time1);
}
