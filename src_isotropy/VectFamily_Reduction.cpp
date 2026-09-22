// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "ClassicLLL.h"
#include "VectFamilyReduction.h"
#include "norms.h"
// clang-format on

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time1;
  try {
    using T = mpq_class;
    if (argc != 3 && argc != 5) {
      std::cerr
          << "VectFamily_Reduction [FileI] [method] [OutFormat] [FileO]\n";
      std::cerr << "or\n";
      std::cerr << "VectFamily_Reduction [FileI] [method]\n";
      std::cerr << "\n";
      std::cerr << " ------- method --------\n";
      std::cerr << "Possible values for method\n";
      std::cerr << "direct      : the direct LLL method\n";
      std::cerr << "dual        : the dual LLL method, through the "
                   "adjugate\n";
      std::cerr << "seysen      : Seysen reduction, minimising\n";
      std::cerr << "              sum_i |b_i|^2 |dual b_i|^2. This treats "
                   "all\n";
      std::cerr << "              the coordinates symmetrically and is the "
                   "one\n";
      std::cerr << "              to try first here: the consumer of the "
                   "output\n";
      std::cerr << "              pays for every coefficient, not for the "
                   "shortest\n";
      std::cerr << "              vector, which is what LLL optimises\n";
      std::cerr << "seysen_best : Seysen with the best move of each sweep, "
                   "slower\n";
      std::cerr << "seysen_lll  : Seysen and LLL alternated\n";
      std::cerr << "deep        : Schnorr-Euchner deep insertion, "
                   "unrestricted\n";
      std::cerr << "deep5       : deep insertion at depth 5\n";
      std::cerr << "deep10      : deep insertion at depth 10\n";
      std::cerr << "best        : run all of the above and keep whichever "
                   "minimises\n";
      std::cerr << "              the facet coefficient estimate, the input "
                   "being\n";
      std::cerr << "              among the candidates so the result is "
                   "never\n";
      std::cerr << "              worse than it. Recommended unless the cost "
                   "of\n";
      std::cerr << "              the reduction itself matters\n";
      std::cerr << "\n";
      std::cerr << " ------- OutFormat --------\n";
      std::cerr << "Possible values for OutFormat\n";
      std::cerr << "GAP : for writing in GAP readable file\n";
      std::cerr << "CPP : for writing in CPP polyhedral readable file\n";
      std::cerr << "\n";
      std::cerr << " ------- FileO --------\n";
      std::cerr << "Possible values for FileO\n";
      std::cerr << "stderr : for writing to std::cerr\n";
      std::cerr << "stdout : for writing in std::cout\n";
      std::cerr << "otherwise written to the named file in output\n";
      throw TerminalException{1};
    }
    std::string FileInput = argv[1];
    MyMatrix<T> M = ReadMatrixFile<T>(FileInput);
    std::string method = argv[2];
    std::cerr << "method=" << method << "\n";
    std::string OutFormat = "CPP";
    std::string FileO = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileO = argv[4];
    }
    auto matrix_measure = [&](MyMatrix<T> const &Minp) -> void {
      T max_coeff = sqr_estimate_facet_coefficients(Minp);
      int64_t max_int64 = std::numeric_limits<int64_t>::max();
      double max_coeff_d = UniversalScalarConversion<double, T>(max_coeff);
      double sqr_max = sqrt(max_coeff_d);
      double max_int64_d = UniversalScalarConversion<double, int64_t>(max_int64);
      T l1_norm = L1_norm_mat(Minp);
      T linf_norm = Linfinity_norm_mat(Minp);
      std::cerr << "max_coeff=" << max_coeff << "\n";
      std::cerr << "max_coeff_d=" << max_coeff_d << "\n";
      std::cerr << "sqr_max=" << sqr_max << "\n";
      std::cerr << "max_int64_d=" << max_int64_d << "\n";
      std::cerr << "L1(M)=" << l1_norm << "\n";
      std::cerr << "Linf(M)=" << linf_norm << "\n";
    };
    std::cerr << "Original complexity measures\n";
    matrix_measure(M);
    VectFamilyReductionResult<T> res =
        ReduceVectorFamilyGeneral(M, method, std::cerr);
    std::cerr << "Output complexity measures\n";
    matrix_measure(res.Mred);
    std::cerr << "method_used=" << res.method << "\n";
    auto print_mat = [&](std::ostream &os_out) -> void {
      if (OutFormat == "GAP") {
        os_out << "return ";
        WriteMatrixGAP(os_out, res.Mred);
        os_out << ";\n";
        return;
      }
      if (OutFormat == "CPP") {
        return WriteMatrix(os_out, res.Mred);
      }
      std::cerr << "No matching format in print_mat. Allowed options: GAP, "
                << "CPP\n";
      throw TerminalException{1};
    };
    FILE_PrintStderrStdoutFile(FileO, print_mat);
    std::cerr << "Normal termination of VectFamily_Reduction\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in VectFamily_Reduction\n";
    exit(e.eVal);
  }
  runtime(time1);
}
