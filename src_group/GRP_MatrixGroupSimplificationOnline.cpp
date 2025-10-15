// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MatrixGroupSimplification.h"
#include "OnlineExhaustiveReduction.h"
// clang-format on

#ifdef TRACK_INFO
#define	TRACK_INFO_ONLINE_INSERTION_SIZES
#endif

template <typename T>
void process(std::string const &FileMatrGroup, std::string const &OutFormat,
             std::ostream &os_out) {
  using Ttype = std::pair<MyMatrix<T>, MyMatrix<T>>;

  std::vector<MyMatrix<T>> ListM = ReadListMatrixFile<T>(FileMatrGroup);
  size_t n_gen = ListM.size();

  // Create the online reduction kernel
  OnlineExhaustiveReductionComplexityMatrixInfinite<T> online_kernel(std::cerr);

  // Insert generators one by one
  for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
    MyMatrix<T> M = ListM[i_gen];
#ifdef TRACK_INFO_ONLINE_INSERTION_SIZES
    bool test = IsIntegralMatrix(M);
    T det = DeterminantMat(M);
    //    std::cerr << "i_gen=" << i_gen << " test=" << test << " det=" << det << "\n";
#endif
    MyMatrix<T> Minv = Inverse(M);
    Ttype pair{M, Minv};

    (void)online_kernel.insert_generator(pair);

#ifdef TRACK_INFO_ONLINE_INSERTION_SIZES
    std::vector<MyMatrix<T>> l_gens = online_kernel.get_current_matrix_t();
    std::cerr << "i_gen=" << i_gen << " " << compute_complexity_listmat(l_gens)	<< "\n";
#endif
  }

  // Extract the final reduced set
  std::vector<MyMatrix<T>> ListMred = online_kernel.get_current_matrix_t();

  // Output the results
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
      std::cerr << "GRP_MatrixGroupPermSimplificationOnline [Arith] [FileMatrGroup]\n";
      std::cerr << "        [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_MatrixGroupPermSimplificationOnline [Arith] [FileMatrGroup]\n";
      std::cerr << "\n";
      std::cerr << "FileMatrGroup : The list of group generators as matrices\n";
      std::cerr << "OutFormat     : CPP or GAP\n";
      std::cerr << "FileOut       : Output file (or stderr)\n";
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
    std::cerr << "Normal termination of GRP_MatrixGroupSimplificationOnline\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_MatrixGroupSimplificationOnline\n";
    exit(e.eVal);
  }
  runtime(time);
}
