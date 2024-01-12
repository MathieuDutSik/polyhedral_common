// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "GRP_GroupFct.h"
#include "Group.h"
#include "Permutation.h"
#include "PolytopeEquiStab.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time1;
  try {
    if (argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "TEST_WeightMatrixLimited [EXTIN]\n";
      std::cerr << "\n";
      std::cerr
          << "EXTIN : The list of vertices (or inequalities for that matter)\n";
      return -1;
    }
    //
    using T = mpq_class;
    using Tint = mpz_class;
    using Tidx = uint16_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, Tint>;
    //
    std::string eFile = argv[1];
    MyMatrix<T> EXT = ReadMatrixFile<T>(eFile);
    size_t len = EXT.rows();
    //
    Tgroup GRP =
        LinPolytope_Automorphism<T, Tgroup>(EXT, std::cerr);
    std::cerr << "|GRP|=" << GRP.size() << "\n";
    //
    int n_iter1 = 20;
    int n_iter2 = 20;
    auto test_WMat =
        [&](const WeightMatrixLimited<true, T> &WMatLimited) -> void {
      for (int iter1 = 0; iter1 < n_iter1; iter1++) {
        std::cerr << "  iter1=" << iter1 << " / " << n_iter1 << "\n";
        Face eFace(len);
        for (Tidx i = 0; i < len; i++) {
          int eVal = Tidx(random()) % 2;
          eFace[i] = eVal;
        }
        size_t hash1 = WMatLimited.get_hash(eFace);
        for (int iter2 = 0; iter2 < n_iter2; iter2++) {
          std::cerr << "    iter2=" << iter2 << " / " << n_iter2 << "\n";
          Telt eElt = GRP.rand();
          Face fFace = OnFace(eFace, eElt);
          size_t hash2 = WMatLimited.get_hash(eFace);
          if (hash1 != hash2) {
            std::cerr << "hash1=" << hash1 << " hash2=" << hash2 << "\n";
            throw TerminalException{1};
          }
        }
      }
    };
    //
    size_t total = len * (len - 1) / 2;
    int expo = 10;
    for (int i = 0; i <= expo; i++) {
      size_t max_offdiag =
          size_t(static_cast<double>(total + 1) *
                 (static_cast<double>(i) / static_cast<double>(10)));
      std::cerr << "i=" << i << " expo=" << expo
                << " max_offdiag=" << max_offdiag << " total=" << total << "\n";
      WeightMatrixLimited<true, T> WMatLimited1 =
          GetWeightMatrixLimited(EXT, max_offdiag, std::cerr);
      test_WMat(WMatLimited1);
      WeightMatrixLimited<true, T> WMatLimited2(GRP, max_offdiag, std::cerr);
      test_WMat(WMatLimited2);
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in TEST_WeightMatrixLimited\n";
    exit(e.eVal);
  }
  runtime(time1);
}
