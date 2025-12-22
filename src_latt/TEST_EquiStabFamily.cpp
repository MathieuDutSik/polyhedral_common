// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "Group.h"
#include "Permutation.h"
#include "EquiStabMemoization.h"
#include "LatticeStabEquiCan.h"
#include "COMB_Combinatorics.h"
#include "GRAPH_GraphicalBasic.h"
#include <unordered_set>
// clang-format on

template <typename T, typename Tint>
void process(std::string const &ListMatFile, std::string const &OutFormat,
             std::ostream &os_out) {
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = mpz_class;
  using Tgroup = permutalib::Group<Telt, TintGroup>;
  using Tgr = GraphListAdj;
  std::vector<MyMatrix<T>> ListMat = ReadListMatrixFile<T>(ListMatFile);
  size_t nMat = ListMat.size();
  //
  // Doing the classic processing and building the test scheme.
  //
  MicrosecondTime time;
  size_t miss_val = std::numeric_limits<size_t>::max();
  std::unordered_set<std::pair<size_t, size_t>> set_equiv;
  std::vector<std::pair<size_t, size_t>> list_cases;
  std::vector<std::pair<size_t, size_t>> l_pair;
  for (size_t iMat = 0; iMat < nMat; iMat++) {
    std::pair<size_t, size_t> pair_stab{iMat, miss_val};
    list_cases.push_back(pair_stab);
    (void)ArithmeticAutomorphismGroup<T, Tint, Tgroup>(ListMat[iMat],
                                                       std::cerr);
    std::cerr << "STAB: iMat=" << iMat << "\n";
    for (size_t jMat = iMat + 1; jMat < nMat; jMat++) {
      std::cerr << "EQUI: Before iMat=" << iMat << " jMat=" << jMat << "\n";
      std::optional<MyMatrix<Tint>> opt = ArithmeticEquivalence<T, Tint, Tgroup>(
          ListMat[iMat], ListMat[jMat], std::cerr);
      std::cerr << "EQUI: After iMat=" << iMat << " jMat=" << jMat << "\n";
      std::pair<size_t, size_t> pair_equiv{iMat, jMat};
      list_cases.push_back(pair_equiv);
      if (opt) {
        set_equiv.insert(pair_equiv);
        l_pair.push_back(pair_equiv);
      }
    }
  }
  Tgr eGR(l_pair, nMat);
  std::vector<std::vector<size_t>> vect_cone = ConnectedComponents_set(eGR);
  std::cerr << "Total runtime for direct computation=" << time << "\n";
  //
  // Now the scheme with the enhanced equivalence
  //
  DatabaseResultEquiStab<MyMatrix<T>, MyMatrix<Tint>> database;
  auto get_stab_inner =
      [&](MyMatrix<T> const &eMat) -> std::vector<MyMatrix<Tint>> {
    std::optional<std::vector<MyMatrix<Tint>>> opt =
        database.attempt_stabilizer(eMat);
    if (opt) {
      return *opt;
    } else {
      std::vector<MyMatrix<Tint>> ListGen =
          ArithmeticAutomorphismGroup<T, Tint, Tgroup>(eMat, std::cerr);
      database.insert_stab(eMat, ListGen);
      return ListGen;
    }
  };
  auto get_stab = [&](size_t const &xMat) -> void {
    MyMatrix<T> eMat = ListMat[xMat];
    std::vector<MyMatrix<Tint>> ListGen = get_stab_inner(eMat);
    for (auto &eEquiv : ListGen) {
      MyMatrix<T> eEquiv_T = UniversalMatrixConversion<T, Tint>(eEquiv);
      MyMatrix<T> eProd = eEquiv_T * eMat * eEquiv_T.transpose();
      if (eProd != eMat) {
        std::cerr << "The matrix does not belong to the stabilizer\n";
        throw TerminalException{1};
      }
    }
  };
  auto get_equiv_inner =
      [&](MyMatrix<T> const &eMat1,
          MyMatrix<T> const &eMat2) -> std::optional<MyMatrix<Tint>> {
    std::optional<std::optional<MyMatrix<Tint>>> opt =
        database.attempt_equiv(eMat1, eMat2);
    if (opt) {
      std::cerr << "get_equiv_inner, from database\n";
      return *opt;
    } else {
      std::cerr << "get_equiv_inner, from direct computation\n";
      std::optional<MyMatrix<Tint>> optB =
        ArithmeticEquivalence<T, Tint, Tgroup>(eMat1, eMat2, std::cerr);
      database.insert_equi(eMat1, eMat2, optB);
      return optB;
    }
  };
  auto get_equiv = [&](size_t const &iMat, size_t const &jMat) -> void {
    MyMatrix<T> eMat1 = ListMat[iMat];
    MyMatrix<T> eMat2 = ListMat[jMat];
    std::optional<MyMatrix<Tint>> opt = get_equiv_inner(eMat1, eMat2);
    std::pair<size_t, size_t> pair{iMat, jMat};
    if (opt) {
      MyMatrix<Tint> const &Eq = *opt;
      MyMatrix<T> Eq_T = UniversalMatrixConversion<T, Tint>(Eq);
      MyMatrix<T> eProd = Eq_T * eMat1 * Eq_T.transpose();
      if (eProd != eMat2) {
        std::cerr << "The matrix is not an equivalene\n";
        throw TerminalException{1};
      }
      if (set_equiv.count(pair) != 1) {
        std::cerr << "This contradicts the previous computation 1\n";
        throw TerminalException{1};
      }
    } else {
      if (set_equiv.count(pair) != 0) {
        std::cerr << "This contradicts the previous computation 2\n";
        throw TerminalException{1};
      }
    }
  };
  size_t n_case = list_cases.size();
  std::vector<size_t> ePerm = RandomPermutation<size_t>(n_case);
  for (size_t iCase = 0; iCase < n_case; iCase++) {
    size_t jCase = ePerm[iCase];
    std::pair<size_t, size_t> eCase = list_cases[jCase];
    size_t iMat = eCase.first;
    size_t jMat = eCase.second;
    if (jMat == miss_val) {
      get_stab(iMat);
    } else {
      get_equiv(iMat, jMat);
    }
  }
  std::cerr << "Total runtime for memoized computation=" << time << "\n";
  //
  if (OutFormat == "GAP") {
    os_out << "return [";
    bool IsFirst = true;
    for (auto &set : vect_cone) {
      if (!IsFirst)
        os_out << ",\n";
      IsFirst = false;
      //
      bool IsFirstB = true;
      os_out << "[";
      for (auto &val1 : set) {
        if (!IsFirstB)
          os_out << ",";
        IsFirstB = false;
        size_t val2 = val1 + 1;
        os_out << val2;
      }
      os_out << "]";
    }
    os_out << "];\n";
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "TEST_EquiStabFamily [arith] [ListMatFile]\n";
      std::cerr << "or\n";
      std::cerr << "TEST_EquiStabFamily [arith] [ListMatFile] [OutFormat] "
                   "[FileOut]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string ListMatFile = argv[2];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }
    //
    auto f = [&](std::ostream &os_out) -> void {
      if (arith == "rational") {
        using T = mpq_class;
        using Tint = mpz_class;
        return process<T, Tint>(ListMatFile, OutFormat, os_out);
      }
      std::cerr << "Failed to find matching type for arith\n";
      throw TerminalException{1};
    };
    print_stderr_stdout_file(FileOut, f);
    std::cerr << "Normal termination of TEST_EquiStabFamily\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_EquiStabFamily\n";
    exit(e.eVal);
  }
  runtime(time);
}
