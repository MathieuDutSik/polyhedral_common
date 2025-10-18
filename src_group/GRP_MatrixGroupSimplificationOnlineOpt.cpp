// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "OnlineExhaustiveReduction.h"
#include "Permutation.h"
// clang-format on

#ifdef TRACK_INFO
#define TRACK_INFO_ONLINE_INSERTION_SIZES
#endif

template <typename T>
void process(std::string const &FileMatrGroup, std::string const &OutFormat,
             std::ostream &os_out) {
  std::vector<MyMatrix<T>> ListM = ReadListMatrixFile<T>(FileMatrGroup);
  size_t n_gen = ListM.size();

  if (n_gen == 0) {
    std::cerr << "No generators found in input file\n";
    return;
  }

  int n = ListM[0].rows();

  // Create the hierarchical online reduction system
  OnlineHierarchicalMatrixReduction<T> hierarchical_reducer(n, std::cerr);

  // Insert generators one by one
  for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
    hierarchical_reducer.insert_generator(ListM[i_gen]);
#ifdef TRACK_INFO_ONLINE_INSERTION_SIZES_DISABLE
    std::cerr << "i_gen=" << i_gen;
    hierarchical_reducer.print_invariants(std::cerr);
#endif
  }

  // Extract the final reduced set
  std::vector<MyMatrix<T>> ListMred = hierarchical_reducer.get_current_matrix_t();

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
      std::cerr << "GRP_MatrixGroupPermSimplificationOnlineOpt [Arith] [FileMatrGroup]\n";
      std::cerr << "        [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_MatrixGroupPermSimplificationOnlineOpt [Arith] [FileMatrGroup]\n";
      std::cerr << "\n";
      std::cerr << "This is the optimized hierarchical online version that automatically\n";
      std::cerr << "switches between int16_t -> int32_t -> int64_t -> T numerics\n";
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
    std::cerr << "Normal termination of GRP_MatrixGroupSimplificationOnlineOpt\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_MatrixGroupSimplificationOnlineOpt\n";
    exit(e.eVal);
  }
  runtime(time);
}
