// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "eutacticity.h"
#include "Arithmetic_dispatch_types.h"
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
  ArithmeticDispatch dispatch;
  add_rational_integer_arithmetics_by_library(
      dispatch, [&]<typename T, typename Tint>() -> void {
        compute_eutacticity_kernel<T, Tint>(eFile, eutacticity);
      });
  dispatch.run(arithmetic);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "IsEutactic [arithmetic] [inputMat] [Eutacticity]\n";
      std::cerr << "\n";
      ArithmeticDispatch dispatch_list;
      add_rational_integer_arithmetics_by_library(
          dispatch_list, []<typename T, typename Tint>() -> void {});
      std::cerr << "arithmetic: " << dispatch_list.choices() << "\n";
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
