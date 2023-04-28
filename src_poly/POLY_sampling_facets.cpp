// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#ifdef OSCAR_USE_BOOST_GMP_BINDINGS
# include "NumberTheoryBoostGmpInt.h"
#else
# include "NumberTheory.h"
#endif
#include "NumberTheoryRealField.h"
#include "QuadField.h"
#include "POLY_SamplingFacet.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::string const& ansSamp, std::string const& OutFormat, std::ostream &os) {
  MyMatrix<T> EXT = ReadMatrixFile<T>(eFileI);
  vectface vf = DirectComputationInitialFacetSet(EXT, ansSamp, std::cerr);
  if (OutFormat == "GAP") {
    os << "return ";
    VectVectInt_Gap_Print(os, vf);
    os << ";\n";
    return;
  }
  if (OutFormat == "Oscar") {
    size_t n_ent = vf.size();
    size_t n = vf.size();
    MyMatrix<int> M(n_ent, n);
    for (size_t i_ent=0; i_ent<n_ent; i_ent++) {
      Face f = vf[i_ent];
      for (size_t i=0; i<n; i++)
        M(i_ent, i) = f[i];
    }
    WriteMatrix(os, M);
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat=" << OutFormat << "\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_sampling_facets arith command [FileI] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "POLY_sampling_facets arith command [FileI]\n";
      std::cerr << "\n";
      std::cerr << "with:\n";
      std::cerr << "arith     : the chosen arithmetic\n";
      std::cerr << "command   : the program used for computing the dual description\n";
      std::cerr << "FileI     : The polyhedral cone inequalities\n";
      std::cerr << "OutFormat : The format of output, GAP or Oscar\n";
      std::cerr << "FileO     : The file of output (if present, otherwise std::cerr)\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr << "RealAlgebraic=FileDesc  : For the real algebraic case of a "
                   "field whose description is in FileDesc\n";
      std::cerr << "\n";
      std::cerr << "        --- command ---\n";
      std::cerr << "\n";
      std::cerr << "lp_cdd                   : 10 facets being sampled by using linear programming\n";
      std::cerr << "lp_cdd:iter_X            : X facets being sampled by using linear programming\n";
      std::cerr << "lp_cdd_min               : 10 facets (or less) of minimum incidence being sampled\n";
      std::cerr << "lp_cdd_min:iter_X        : X facets (or less) of minimum incidence being sampled\n";
      std::cerr << "sampling                 : The recursive sampling algorithm\n";
      std::cerr << "lrs_limited              : lrs limited to the first 100 vertices being found\n";
      std::cerr << "lrs_limited:upperlimit_X : lrs limited to the first X vertices being found\n";
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
      if (arith == "rational") {
        return process<Trat>(eFileI, command, OutFormat, os);
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
    if (eFileO == "stderr") {
      call_lrs(std::cerr);
    } else {
      if (eFileO == "stdout") {
        call_lrs(std::cout);
      } else {
        std::ofstream os(eFileO);
        call_lrs(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_dual_description\n";
    exit(e.eVal);
  }
  runtime(time1);
}
