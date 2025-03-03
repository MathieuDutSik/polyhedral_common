// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "POLY_LinearProgramming.h"
// clang-format on

template <typename T>
void process(std::string const &eFileI, std::ostream &os_out,
             std::ostream &os) {
  MyMatrix<T> FAC = ReadMatrixFile<T>(eFileI);
  Face f_adj = ComputeSkeletonClarkson(FAC, os);
  int n_fac = FAC.rows();
  for (int i_fac = 0; i_fac < n_fac; i_fac++) {
    int n_adj = 0;
    for (int j_fac = 0; j_fac < n_fac; j_fac++) {
      int val = f_adj[j_fac + i_fac * n_fac];
      os_out << " " << val;
      n_adj += val;
    }
    os_out << " n_adj=" << n_adj << "\n";
  }
}

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "POLY_SkelettonClarkson arith [DATAIN]\n";
      std::cerr << "\n";
      std::cerr << "        --- arith ---\n";
      std::cerr << "\n";
      std::cerr << "safe_rational : rational arithmetic based on int64_t that "
                   "fails\n";
      std::cerr << "   gracefully if overflowing\n";
      std::cerr << "rational : rational arithmetic on input\n";
      std::cerr << "Qsqrt2   : arithmetic over the field Q(sqrt(2))\n";
      std::cerr << "Qsqrt5   : arithmetic over the field Q(sqrt(5))\n";
      std::cerr
          << "RealAlgebraic=FileDesc  : For the real algebraic case of a\n";
      std::cerr << "   field whose description is in FileDesc\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string eFileI = argv[2];
    auto compute_skeleton = [&](std::ostream &os_out) -> void {
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(eFileI, os_out, std::cerr);
      }
      if (arith == "rational") {
        using T = mpq_class;
        return process<T>(eFileI, os_out, std::cerr);
      }
      if (arith == "Qsqrt5") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 5>;
        return process<T>(eFileI, os_out, std::cerr);
      }
      if (arith == "Qsqrt2") {
        using Trat = mpq_class;
        using T = QuadField<Trat, 2>;
        return process<T>(eFileI, os_out, std::cerr);
      }
      std::optional<std::string> opt_realalgebraic =
          get_postfix(arith, "RealAlgebraic=");
      if (opt_realalgebraic) {
        using T_rat = mpq_class;
        std::string FileAlgebraicField = *opt_realalgebraic;
        if (!IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        using T = RealField<idx_real_algebraic_field>;
        return process<T>(eFileI, os_out, std::cerr);
      }
      std::cerr << "Failed to find a matching field for arith=" << arith
                << "\n";
      std::cerr << "Available possibilities: rational, Qsqrt5, Qsqrt2, "
                   "RealAlgebraic\n";
      throw TerminalException{1};
    };
    compute_skeleton(std::cerr);
    std::cerr << "Normal termination of POLY_SkelettonClarkson\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLY_SkelettonClarkson\n";
    exit(e.eVal);
  }
  runtime(time1);
}
