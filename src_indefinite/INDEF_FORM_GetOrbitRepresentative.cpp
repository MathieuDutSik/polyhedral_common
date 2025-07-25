// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "CombinedAlgorithms.h"
#include "Group.h"
#include "Permutation.h"
// clang-format on

template <typename T, typename Tint, typename Tgroup>
void process(std::string const &MatFile, std::string const &XnormStr,
             std::string const &OutFormat, std::ostream &os_out) {
  MyMatrix<T> Qmat = ReadMatrixFile<T>(MatFile);
  T Xnorm = ParseScalar<T>(XnormStr);
  IndefiniteCombinedAlgo<T, Tint, Tgroup> comb(std::cerr);
  std::vector<MyVector<Tint>> l_repr =
      comb.INDEF_FORM_GetOrbitRepresentative(Qmat, Xnorm);
  for (auto & e_repr: l_repr) {
    T norm = EvaluationQuadForm(Qmat, e_repr);
    if (norm != Xnorm) {
      std::cerr << "Wrong norm. norm=" << norm << " Xnorm=" << Xnorm << "\n";
      throw TerminalException{1};
    }
  }
  if (OutFormat == "PYTHON") {
    if (l_repr.size() == 0) {
      os_out << "[]";
    } else {
      WriteMatrixPYTHON(os_out, MatrixFromVectorFamily(l_repr));
    }
    return;
  }
  if (OutFormat == "GAP") {
    os_out << "return ";
    if (l_repr.size() == 0) {
      os_out << "[]";
    } else {
      WriteMatrixGAP(os_out, MatrixFromVectorFamily(l_repr));
    }
    os_out << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr
          << "INDEF_FORM_GetOrbitRepresentative [arith] [MatFile] [Xnorm] \n";
      std::cerr << "or\n";
      std::cerr << "INDEF_FORM_GetOrbitRepresentative [arith] [MatFile] "
                   "[Xnorm] [OutFormat] [OutFile]\n";
      std::cerr << "     -----\n";
      std::cerr << "arith can be: gmp\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string MatFile = argv[2];
    std::string XnormStr = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using TintGroup = mpz_class;
    using Tgroup = permutalib::Group<Telt, TintGroup>;
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint, Tgroup>(MatFile, XnormStr, OutFormat, os);
      }
      std::cerr << "Failed to find matching type for arith\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(OutFile, f);
    std::cerr << "Normal termination of INDEF_FORM_GetOrbitRepresentative\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in INDEF_FORM_GetOrbitRepresentative, runtime=" << time << "\n";
    exit(e.eVal);
  }
  runtime(time);
}
