// Copyright (C) 202" Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryGmp.h"
#include "ApproximateModels.h"
#include "Group.h"
#include "Permutation.h"
#include "EquiStabMemoization.h"
#include "LatticeStabEquiCan.h"
// clang-format on

template <typename T, typename Tint>
void process(std::string const &ListMatFile) {
  std::vector<MyMatrix<T>> ListMat = ReadListMatrixFile<T>(MatFile);
  size_t nMat = ListMat.size();
  //
  // Doing the classic processing and building the test scheme.
  //
  MicrosecondTime time;
  size_t miss_val = std::numeric_limits<size_t>::max();
  std::unordered_map<std::pair<size_t, size_t>> map_equiv;
  std::vector<std::pair<size_t, size_t>> list_cases;
  for (size_t iMat=0; iMat<nMat; iMat++) {
    std::pair<<size_t, size_t> pair_stab{iMat, miss_val};
    list_cases.push_back(pair_stab);
    (void)ArithmeticAutomorphismGroup(ListMat[iMat], std::cerr);
    for (size_t jMat=iMat+1; jMat<nMat; jMat++) {
      std::optional<MyMatrix<Tint>> opt = ArithmeticEquivalence(ListMat[iMat], ListMat[jMat], std::cerr);
      std::pair<<size_t, size_t> pair_equiv{iMat, jMat};
      list_cases.push_back(pair_equiv);
      if (opt) {
        map_equiv.insert(pair_equiv);
      }
    }
  }
  std::cerr << "Total runtime for direct computation=" << time << "\n";
  //
  // Now the scheme with the enhanced equivalence
  //
  DatabaseResultEquiStab<MyMatrix<T>,MyMatrix<Tint>> database;
  auto get_stab_inner=[&](MyMatrix<T> const& eMat) -> std::vector<MyMatrix<Tint>> {
    std::optional<std::vector<Tequiv>> opt = database.attempt_stabilizer(eMat);
    if (opt) {
      return *opt;
    } else {
      std::vector<MyMatrix<Tint>> ListGen = ArithmeticAutomorphismGroup(eMat, std::cerr);
      database.insert_stab(eMat, ListGen);
      return ListGen;
    }
  };
  auto get_stab=[&](size_t const& xMat) -> void {
    MyMatrix<T> eMat = ListMat[xMat];
    std::vector<MyMatrix<Tint>> ListGen = get_stab_inner(eMat);
    for (auto & eEquiv : ListGen) {
      MyMatrix<T> eEquiv_T = UniversalMatrixConversion<T,Tint>(eEquiv);
      MyMatrix<T> eProd = eEquiv_T * eMat * eEquiv_T.transpose();
      if (eProd != eMat) {
        std::cerr << "The matrix does not belong to the stabilizer\n";
        throw TerminalException{1};
      }
    }
  };
  auto get_equiv_inner=[&](MyMatrix<T> const& eMat1, MyMatrix<T> const& eMat2) -> std::optional<MyMatrix<Tint>> {
    std::optional<std::optional<Tequiv>> opt = database.attempt_equiv(eMat);
    if (opt) {
      return *opt;
    } else {
      std::optional<MyMatrix<Tint>> optB = ArithmeticEquivalence(eMat1, eMat2, std::cerr);
      database.insert_equi(eMat1, eMat2, optB);
      return optB;
    }
  };
  auto get_equiv=[&](size_t const& iMat, size_t const& jMat) -> void {
    MyMatrix<T> eMat1 = ListMat[iMat];
    MyMatrix<T> eMat2 = ListMat[jMat];
    std::optional<MyMatrix<Tint>> opt = get_equiv_inner(eMat1, eMat2);
    std::pair<size_t, size_t> pair{iMat, jMat};
    if (opt) {
      MyMatrix<Tint> const& Eq = *opt;
      MyMatrix<T> Eq_T = UniversalMatrixConversion<T,Tint>(Eq);
      MyMatrix<T> eProd = Eq_T * eMat1 * Eq_T.transpose();
      if (eProd != eMat2) {
        std::cerr << "The matrix is not an equivalene\n";
        throw TerminalException{1};
      }
      if (map_equiv.count(pair) != 1) {
        std::cerr << "This contradicts the previous computation 1\n";
        throw TerminalException{1};
      }
    } else {
      if (map_equiv.count(pair) != 0) {
        std::cerr << "This contradicts the previous computation 2\n";
        throw TerminalException{1};
      }
    }
  }
  for (auto & eCase : list_cases) {
    size_t iMat = eCase.first;
    size_t jMat = eCase.second;
    if (jMat == miss_val) {
      get_stab(iMat);
    } else {
      get_equiv(iMat, jMat);
    }
  }
  std::cerr << "Total runtime for memoized computation=" << time << "\n";
}

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 3) {
      std::cerr << "TEST_EquiStabFamily [arith] [ListMatFile]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    std::string ListMatFile = argv[2];
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tint = mpz_class;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    auto f = [&]() -> void {
      if (arith == "rational") {
        using T = mpq_class;
        return process<T,Tgroup>(MatFile);
      }
      std::cerr << "Failed to find matching type for arith\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_FindIsotropic\n";
    exit(e.eVal);
  }
  runtime(time1);
}
