// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MatrixGroupSimplification.h"
#include "OnlineExhaustiveReduction.h"
#include "Permutation.h"
// clang-format on

template <typename T>
void process(std::string const &FileMatrGroup, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using Ttype = std::pair<MyMatrix<T>, Telt>;
  using Tnorm = T;

  std::vector<MyMatrix<T>> ListM = ReadListMatrixFile<T>(FileMatrGroup);
  size_t n_gen = ListM.size();

  // Define the functions needed for the online algorithm
  auto f_complexity=[&](Ttype const& pair) -> T {
    return get_ell1_complexity_measure(pair.first);
  };
  auto f_product=[](Ttype const& p1, Ttype const& p2) -> Ttype {
    return {p1.first * p2.first, p1.second * p2.second};
  };
  auto f_check=[&]([[maybe_unused]] Ttype const& p) -> bool {
    return true;
  };

  // Create the online reduction kernel
  OnlineExhaustiveReductionComplexityKernel<Ttype, Tnorm, decltype(f_complexity), decltype(f_product), decltype(f_check)> online_kernel(f_complexity, f_product, f_check, std::cerr);

  // Insert generators one by one
  std::cerr << "Inserting " << n_gen << " generators one by one...\n";
  for (size_t i_gen = 0; i_gen < n_gen; i_gen++) {
    Telt elt; // Default permutation element
    Ttype pair{ListM[i_gen], elt};

    // Create the TcombPair structure needed by the online kernel
    auto f_invers = [](Ttype const& pair) -> Ttype {
      return {Inverse(pair.first), Inverse(pair.second)};
    };
    Ttype pair_inv = f_invers(pair);
    std::pair<Ttype,Ttype> pair_di{pair, pair_inv};

    bool success = online_kernel.insert_generator(pair_di);
    if (!success) {
      std::cerr << "Failed to insert generator " << i_gen << "\n";
      throw TerminalException{1};
    }

    if (i_gen % 10 == 9 || i_gen == n_gen - 1) {
      std::cerr << "After inserting " << (i_gen + 1) << " generators: "
                << "current set size = " << online_kernel.size()
                << ", total complexity = " << online_kernel.get_total_complexity() << "\n";
    }
  }

  // Extract the final reduced set
  std::vector<TcombPair<Ttype, Tnorm>> final_set = online_kernel.get_current_set();
  std::vector<MyMatrix<T>> ListMred;

  for (auto const& comb_pair : final_set) {
    ListMred.push_back(comb_pair.pair.first.first); // Extract the matrix from pair.first
  }

  std::cerr << "Original generators: " << n_gen << ", Reduced generators: " << ListMred.size() << "\n";

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
