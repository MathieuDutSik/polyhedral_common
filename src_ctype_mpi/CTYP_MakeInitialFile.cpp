#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "CTYP_MakeInitialFile [FileIn] [FileOut]\n";
      throw TerminalException{1};
    }
    using Tmat=int;
    //
    std::string FileIn  = argv[1];
    std::string FileOut = argv[2];
    //
    std::ifstream is(FileIn);
    std::ofstream os(FileOut);
    int nbType;
    is >> nbType;
    os << nbType << "\n";
    for (int iType=0; iType<nbType; iType++) {
      std::cerr << "iType : " << iType << " / " << nbType << "\n";
      MyMatrix<Tmat> ePerfect_Tmat = ReadMatrix<Tmat>(is);
      //
      int eStatus=0;
      //
      std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
      MyMatrix<Tmat> eMatCan_Tmat = LinPolytopeAntipodalIntegral_CanonicForm<Tmat>(ePerfect_Tmat);
      std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
      int elapsed_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
      std::cerr << " elapsed_time=" << elapsed_time << "\n";
      //
      os << eStatus << "\n";
      WriteMatrix(os, eMatCan_Tmat);
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
