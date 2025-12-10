// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MatrixGroupSimplification.h"
#include "Permutation.h"
// clang-format on

template <typename T>
void process(std::string const &FileMatrGroup, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  std::vector<MyMatrix<T>> ListM = ReadListMatrixFile<T>(FileMatrGroup);
  size_t n_gen = ListM.size();
  std::vector<Telt> ListPerm;
  for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
    Telt elt;
    ListPerm.push_back(elt);
  }
  std::pair<std::vector<MyMatrix<T>>, std::vector<Telt>> retpair =
      ExhaustiveReductionComplexityGroupMatrixPerm<T, Telt>(ListM, ListPerm,
                                                            std::cerr);
  std::vector<MyMatrix<T>> const &ListMred = retpair.first;
  //
  if (OutFormat == "GAP") {
    os_out << "return ";
    WriteListMatrixGAP(os_out, ListMred);
    os_out << ";\n";
    return;
  }
  if (OutFormat == "CPP") {
    WriteListMatrix(os_out, ListMred);
    return;
  }
  std::cerr << "Failed to find a matching OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_MatrixGroupSimplification [Arith] [FileMatrGroup]\n";
      std::cerr << "        [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_LinPolytope_Invariant [Arith] [FileMatrGroup]\n";
      std::cerr << "\n";
      std::cerr << "EXTIN  : The list of vertices\n";
      std::cerr << "OutCan : The file for the hash\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string FileMatrGroup = argv[2];
    std::string OutFormat = "CPP";
    std::string FileOut = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }

    auto f = [&](std::ostream &os) -> void {
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(FileMatrGroup, OutFormat, os);
      }
      if (arith == "mpz_class") {
        using T = mpz_class;
        return process<T>(FileMatrGroup, OutFormat, os);
      }
      std::cerr << "Failed to find a matching arith\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileOut, f);
    std::cerr << "Normal termination of GRP_MatrixGroupSimplification\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_LinPolytope_Invariant\n";
    exit(e.eVal);
  }
  runtime(time);
}
