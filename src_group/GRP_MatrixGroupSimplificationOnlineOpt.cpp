// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "OnlineExhaustiveReduction.h"
#include "Arithmetic_dispatch_types.h"
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
  }

  // Extract the final reduced set
  std::vector<MyMatrix<T>> ListMred =
      hierarchical_reducer.get_current_matrix_t();

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
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    std::string FileMatrGroup;
    std::string OutFormat = "CPP";
    std::string FileOut = "stderr";
    std::ostream *os_ptr = nullptr;
    ArithmeticDispatch dispatch;
    auto f_arith = [&]<typename T>() -> void {
      process<T>(FileMatrGroup, OutFormat, *os_ptr);
    };
    add_rational_arithmetics(dispatch, f_arith);
    add_integer_arithmetics(dispatch, f_arith);
    if (argc != 3 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRP_MatrixGroupSimplificationOnlineOpt [Arith] "
                   "[FileMatrGroup]\n";
      std::cerr << "        [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "GRP_MatrixGroupSimplificationOnlineOpt [Arith] "
                   "[FileMatrGroup]\n";
      std::cerr << "\n";
      std::cerr << "This is the optimized hierarchical online version that "
                   "automatically\n";
      std::cerr
          << "switches between int16_t -> int32_t -> int64_t -> T numerics\n";
      std::cerr << "\n";
      std::cerr << "FileMatrGroup : The list of group generators as matrices\n";
      std::cerr << "OutFormat     : CPP or GAP\n";
      std::cerr << "FileOut       : Output file (or stderr)\n";
      std::cerr << "\n";
      std::cerr << "     ------- Arith -------\n";
      std::cerr << "\n";
      dispatch.print_descriptions(std::cerr);
      return -1;
    }
    //
    std::string arith = argv[1];
    FileMatrGroup = argv[2];
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }

    auto f = [&](std::ostream &os) -> void {
      os_ptr = &os;
      dispatch.run(arith);
    };
    FILE_PrintStderrStdoutFile(FileOut, f);
    std::cerr
        << "Normal termination of GRP_MatrixGroupSimplificationOnlineOpt\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in GRP_MatrixGroupSimplificationOnlineOpt\n";
    exit(e.eVal);
  }
  runtime(time);
}
