// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "InvariantVectorFamily.h"
// clang-format on

template <typename T, typename Tint>
void process_C(std::string choice, std::string MatFile, std::string const& OutFormat, std::string const& OutFile) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(MatFile);
  auto f = [&]() -> MyMatrix<Tint> {
    if (choice == "shortest") {
      Tshortest<T, Tint> rec = T_ShortestVector<T,Tint>(GramMat, std::cerr);
      return rec.SHV;
    }
    if (choice == "relevant_voronoi") {
      return ComputeVoronoiRelevantVector<T, Tint>(GramMat, std::cerr);
    }
    if (choice == "filtered_relevant_voronoi") {
      MyMatrix<Tint> M =
          ComputeVoronoiRelevantVector<T, Tint>(GramMat, std::cerr);
      return FilterByNorm(GramMat, M, std::cerr);
    }
    if (choice == "fullrank") {
      return ExtractInvariantVectorFamilyFullRank<T, Tint>(GramMat, std::cerr);
    }
    if (choice == "spanning") {
      return ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat, std::cerr);
    }
    std::cerr << "Failed to find a matching entry for choice\n";
    std::cerr << "Possible choices: shortest, relevant_voronoi, filtered_relevant_voronoi, fullrank, spanning\n";
    throw TerminalException{1};
  };
  MyMatrix<Tint> M = f();
#ifdef SANITY_CHECK
  check_antipodality_mymatrix(M);
#endif
  auto f_print=[&](std::ostream& osf) -> void {
    if (OutFormat == "norms") {
      int n_vect = M.rows();
      int n = GramMat.cols();
      std::cerr << "|M|=" << n_vect << " / " << n << "\n";
      std::map<T, size_t> map;
      for (int i_vect = 0; i_vect < n_vect; i_vect++) {
        MyVector<Tint> V = GetMatrixRow(M, i_vect);
        T norm = EvaluationQuadForm<T,Tint>(GramMat, V);
        map[norm] += 1;
      }
      osf << "norms =";
      for (auto &kv : map) {
        osf << " [" << kv.first << " : " << kv.second << " ]";
      }
      osf << "\n";
      return;
    }
    if (OutFormat == "GAP") {
      osf << "return ";
      WriteMatrixGAP(osf, M);
      osf << ";\n";
      return;
    }
    if (OutFormat == "CPP") {
      return WriteMatrix(osf, M);
    }
    std::cerr << "Failed to find a matching entry for OutFormat\n";
    std::cerr << "Allowed choices: norms, GAP, CPP\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(OutFile, f_print);
}

template <typename T>
void process_B(std::string const &arithmetic_vec, std::string choice,
               std::string MatFile, std::string OutFormat, std::string OutFile) {
  if (arithmetic_vec == "mpz_class") {
    using Tint = mpz_class;
    return process_C<T, Tint>(choice, MatFile, OutFormat, OutFile);
  }
  std::cerr << "process_A failure: No matching entry for arithmetic_mat\n";
  throw TerminalException{1};
}

void process_A(std::string const &arithmetic_mat,
               std::string const &arithmetic_vec, std::string choice,
               std::string MatFile, std::string OutFormat, std::string OutFile) {
  if (arithmetic_mat == "mpq_class") {
    using T = mpq_class;
    return process_B<T>(arithmetic_vec, choice, MatFile, OutFormat, OutFile);
  }
  std::cerr << "process_A failure: No matching entry for arithmetic_mat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 7 && argc != 5) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_GenerateCharacteristicVectorSet arithmetic_mat "
                   "arithmetic_vect choice [MaFile] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "LATT_GenerateCharacteristicVectorSet arithmetic_mat "
                   "arithmetic_vect choice [MaFile]\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "arithmetic_mat: mpq_class\n";
      std::cerr << "arithmetic_vect: mpz_class\n";
      std::cerr << "choice: shortest, relevant_voronoi, filtered_relevant_voronoi, fullrank, spanning\n";
      std::cerr << "OutFormat: norms, GAP, CPP\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arithmetic_mat = argv[1];
    std::string arithmetic_vec = argv[2];
    std::string choice = argv[3];
    std::string MatFile = argv[4];
    std::string OutFormat = "norms";
    std::string OutFile = "stderr";
    if (argc == 7) {
      OutFormat = argv[5];
      OutFile = argv[6];
    }
    process_A(arithmetic_mat, arithmetic_vec, choice, MatFile, OutFormat, OutFile);
    std::cerr << "Normal termination of LATT_GenerateCharacteristicVectorSet\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_GenerateCharacteristicVectorSet\n";
    exit(e.eVal);
  }
  runtime(time);
}
