// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "InvariantVectorFamily.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string choice, std::string MatFile,
             std::string const &OutFormat, std::string const &OutFile) {
  MyMatrix<T> GramMat = ReadMatrixFile<T>(MatFile);
  auto f = [&]() -> MyMatrix<Tint> {
    if (choice == "shortest") {
      Tshortest<T, Tint> rec = T_ShortestVector<T, Tint>(GramMat, std::cerr);
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
    std::cerr << "Possible choices: shortest, relevant_voronoi, "
                 "filtered_relevant_voronoi, fullrank, spanning\n";
    throw TerminalException{1};
  };
  MyMatrix<Tint> M = f();
#ifdef SANITY_CHECK
  check_antipodality_mymatrix(M);
#endif
  auto f_print = [&](std::ostream &osf) -> void {
    if (OutFormat == "norms") {
      int n_vect = M.rows();
      int n = GramMat.cols();
      std::cerr << "|M|=" << n_vect << " / " << n << "\n";
      std::map<T, size_t> map;
      for (int i_vect = 0; i_vect < n_vect; i_vect++) {
        MyVector<Tint> V = GetMatrixRow(M, i_vect);
        T norm = EvaluationQuadForm<T, Tint>(GramMat, V);
        map[norm] += 1;
      }
      osf << "norms =";
      for (auto &[norm, multiplicity] : map) {
        osf << " [" << norm << " : " << multiplicity << " ]";
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

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 6 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_GenerateCharacteristicVectorSet [arith] choice [MatFile] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "LATT_GenerateCharacteristicVectorSet [arith] choice [MatFile]\n";
      std::cerr << "allowed choices:\n";
      std::cerr << "[arith]: gmp, gmp_boost, multi_boost\n";
      std::cerr << "choice: shortest, relevant_voronoi, "
                   "filtered_relevant_voronoi, fullrank, spanning\n";
      std::cerr << "OutFormat: norms, GAP, CPP\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string choice = argv[2];
    std::string MatFile = argv[3];
    std::string OutFormat = "norms";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    auto f=[&]() -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T,Tint>(choice, MatFile, OutFormat, OutFile);
      }
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        using Tint = boost::multiprecision::mpz_int;
        return process<T,Tint>(choice, MatFile, OutFormat, OutFile);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        using Tint = boost::multiprecision::cpp_int;
        return process<T,Tint>(choice, MatFile, OutFormat, OutFile);
      }
      if (arith == "safe") {
        using T = Rational<SafeInt64>;
        using Tint = SafeInt64;
        return process<T,Tint>(choice, MatFile, OutFormat, OutFile);
      }
      std::cerr << "process_A failure: No matching entry for arith\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of LATT_GenerateCharacteristicVectorSet\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_GenerateCharacteristicVectorSet\n";
    exit(e.eVal);
  }
  runtime(time);
}
