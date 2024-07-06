// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "ApproximateModels.h"
// clang-format on

template <typename T, typename Tgroup>
void process(std::string const &MatFile, std::string const& XnormStr, std::string const &OutFormat,
             std::ostream &os_out) {
  MyMatrix<T> Qmat = ReadMatrixFile<T>(MatFile);
  std::cerr << "We have Q\n";
  T Xnorm = ParseScalar<T>(XnormStr);
  std::cerr << "We have Xnorm\n";
  ApproximateModel<T> model = INDEF_FORM_EichlerCriterion_TwoHyperplanesEven(Qmat);
  std::vector<MyVector<T>> LVect = model.GetCoveringOrbitRepresentatives(Xnorm);
  if (OutFormat == "GAP") {
    if (LVect.size() == 0) {
      os_out << "return rec(LVect:=[]);\n";
    } else {
      MyMatrix<T> MatVect = MatrixFromVectorFamily(LVect);
      os_out << "return rec(LVect:=" << StringVectorGAP(MatVect) << ");\n";
    }
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "INDEF_ApproximateOrbitRepresentative [arith] [MatFile] [X]\n";
      std::cerr << "or\n";
      std::cerr << "INDEF_ApproximateOrbitRepresentative [arith] [MatFile] [X] [OutFormat] [OutFile]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string MatFile = argv[2];
    std::string Xstr = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    auto f = [&](std::ostream &os) -> void {
      if (arith == "rational") {
        using T = mpq_class;
        return process<T,Tgroup>(FileI, OutFormat, os);
      }
      std::cerr << "Failed to find matching type for arith\n";
      throw TerminalException{1};
    };
    if (FileO == "stderr") {
      f(std::cerr);
    } else {
      if (FileO == "stdout") {
        f(std::cout);
      } else {
        std::ofstream os(FileO);
        f(os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_FindIsotropic\n";
    exit(e.eVal);
  }
  runtime(time1);
}
