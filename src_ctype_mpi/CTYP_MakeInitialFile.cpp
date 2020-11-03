#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PolytopeEquiStab.h"


int main(int argc, char* argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "MakeInitial_FileMatrix [FileIn] [FileOut]\n";
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
    std::cerr << "nbType=" << nbType << "\n";
    for (int iType=0; iType<nbType; iType++) {
      MyMatrix<Tmat> ePerfect_Tmat = ReadMatrix<Tmat>(is);
      //
      int eStatus=0;
      //
      std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
      MyMatrix<Tmat> eMatCan_Tmat = LinPolytopeIntegral_CanonicForm<Tmat>(ePerfect_Tmat);
      std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
      int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
      std::cerr << "iType=" << iType << " / " << nbType << " elapsed_seconds=" << elapsed_seconds << "\n";
      //
      os << eStatus << "\n";
      WriteMatrix(os, eMatCan_Tmat);
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
