// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "PerfectForm.h"
// clang-format on

template <typename T, typename Tint>
void compute_eutacticity_kernel(std::string const &eFile,
                                std::string const &eutacticity) {
  MyMatrix<T> eGram = ReadMatrixFile<T>(eFile);
  Tshortest<T, Tint> rec_shv = T_ShortestVector<T, Tint>(eGram, std::cerr);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(rec_shv.SHV);
  std::optional<MyVector<T>> opt =
      IsEutactic(eGram, SHV_T, eutacticity, std::cerr);
  if (opt) {
    std::cerr << "The matrix is " << eutacticity << "\n";
  } else {
    std::cerr << "The matrix is not " << eutacticity << "\n";
  }
}

void compute_eutacticity(std::string const &arithmetic,
                         std::string const &eFile,
                         std::string const &eutacticity) {
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return compute_eutacticity_kernel<T, Tint>(eFile, eutacticity);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "DeterminantMat [arithmetic] [inputMat] [Eutacticity]\n";
      std::cerr << "\n";
      std::cerr << "Eutacticity: Eutactic, WeaklyEutactic\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string eFile = argv[2];
    std::string eutacticity = argv[3];
    compute_eutacticity(arithmetic, eFile, eutacticity);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
