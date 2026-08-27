// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "Shvec_exact.h"
#include "strongly_semi_eutactic.h"
// clang-format on

template <typename T, typename Tint>
void test_strongly_semi_eutactic_kernel(std::string const &eFile,
                                        size_t max_node) {
  MyMatrix<T> eGram = ReadMatrixFile<T>(eFile);
  Tshortest<T, Tint> rec_shv = T_ShortestVector<T, Tint>(eGram, std::cerr);
  MyMatrix<T> SHV_T = UniversalMatrixConversion<T, Tint>(rec_shv.SHV);
  std::cerr << "min=" << rec_shv.min << " |SHV|=" << SHV_T.rows() << "\n";
  StronglySemiEutacticResult<T> result =
      TestStronglySemiEutactic(eGram, SHV_T, max_node, std::cerr);
  std::cerr << "n_node=" << result.n_node << "\n";
  if (!result.resolved) {
    std::cerr << "The determination is UNRESOLVED: the node budget of "
              << max_node << " was exhausted\n";
    return;
  }
  if (result.cert) {
    size_t siz = result.cert->subset.count();
    std::cerr << "The matrix is strongly semi-eutactic\n";
    std::cerr << "mu=" << result.cert->mu << " |S|=" << siz << " (out of "
              << SHV_T.rows() << " vectors)\n";
    std::cerr << "S={";
    bool is_first = true;
    for (size_t i = 0; i < result.cert->subset.size(); i++) {
      if (result.cert->subset[i] == 1) {
        if (!is_first)
          std::cerr << ",";
        std::cerr << i;
        is_first = false;
      }
    }
    std::cerr << "}\n";
  } else {
    std::cerr << "The matrix is not strongly semi-eutactic\n";
  }
}

void test_strongly_semi_eutactic(std::string const &arithmetic,
                                 std::string const &eFile, size_t max_node) {
  if (arithmetic == "gmp") {
    using T = mpq_class;
    using Tint = mpz_class;
    return test_strongly_semi_eutactic_kernel<T, Tint>(eFile, max_node);
  }
  if (arithmetic == "gmp_boost") {
    using T = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    return test_strongly_semi_eutactic_kernel<T, Tint>(eFile, max_node);
  }
  if (arithmetic == "multi_boost") {
    using T = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    return test_strongly_semi_eutactic_kernel<T, Tint>(eFile, max_node);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 3 && argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "IsStronglySemiEutactic [arithmetic] [inputMat] "
                   "[max_node]\n";
      std::cerr << "\n";
      std::cerr << "arithmetic: gmp, gmp_boost, multi_boost\n";
      std::cerr << "max_node (optional): the node budget of the search, "
                   "default 1000000\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string eFile = argv[2];
    size_t max_node = 1000000;
    if (argc == 4) {
      std::string max_node_str = argv[3];
      max_node = ParseScalar<size_t>(max_node_str);
    }
    test_strongly_semi_eutactic(arithmetic, eFile, max_node);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
