#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "Temp_PolytopeEquiStab.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "CTYP_MakeInitialFile [FileIn] [FileOut]\n";
      throw TerminalException{1};
    }
    using Tmat = int;
    //
    std::string FileIn = argv[1];
    std::string FileOut = argv[2];
    //
    uint32_t seed = 0x1b873540;
    std::ifstream is(FileIn);
    std::ofstream os(FileOut);
    int nbType;
    is >> nbType;
    os << nbType << "\n";
    std::vector<size_t> ListHash(nbType);
    for (int iType = 0; iType < nbType; iType++) {
      std::cerr << "iType : " << iType << " / " << nbType << "\n";
      MyMatrix<Tmat> ePerfect_Tmat = ReadMatrix<Tmat>(is);
      //
      int eStatus = 0;
      //
      std::chrono::time_point<std::chrono::system_clock> start =
          std::chrono::system_clock::now();
      MyMatrix<Tmat> eMatCan_Tmat =
          LinPolytopeAntipodalIntegral_CanonicForm<Tmat>(ePerfect_Tmat);
      std::chrono::time_point<std::chrono::system_clock> end =
          std::chrono::system_clock::now();
      int elapsed_time =
          std::chrono::duration_cast<std::chrono::microseconds>(end - start)
              .count();
      std::cerr << " elapsed_time=" << elapsed_time << "\n";
      //
      os << eStatus << "\n";
      size_t e_hash = Matrix_Hash(eMatCan_Tmat, seed);
      ListHash[iType] = e_hash;
      WriteMatrix(os, eMatCan_Tmat);
    }
    for (size_t n_pes = 1; n_pes < 10; n_pes++) {
      //      std::cerr << "step 1\n";
      std::vector<size_t> ListPart(n_pes, 0);
      //      std::cerr << "step 2\n";
      for (int iType = 0; iType < nbType; iType++) {
        //        std::cerr << "step 3\n";
        size_t res = ListHash[iType] % n_pes;
        //        std::cerr << "step 4\n";
        ListPart[res]++;
      }
      //      std::cerr << "step 5\n";
      std::cerr << "n_pes=" << n_pes << " LPart =";
      for (auto &eV : ListPart)
        std::cerr << " " << eV;
      std::cerr << "\n";
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
