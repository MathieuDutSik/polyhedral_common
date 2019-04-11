#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PerfectForm.h"



struct TypeIndex {
  int iProc;
  int idxMatrix;
  int iAdj;
};


int main(int argc, char* argv[])
{
  using T=mpq_class;
  using Tint=int;
  //
  try {
    if (argc != 3) {
      std::cerr << "MargePartialLogs [FileIn] [PrefixLog] [FileOut]\n";
      throw TerminalException{1};
    }
    std::string FileIn    = argv[1];
    std::string PrefixLog = argv[2];
    std::string FileOut   = argv[3];
    //
    std::string FileMatrix = BlDATA.ListStringValues.at("ListMatrixInput");
    //
    std::vector<MyMatrix<Tint>> ListNewMatrices;
    int iProc=0;
    while(true) {
      std::stringstream s;
      s << PrefixLog + iProc;
      string eFileLog(s.str());
      //
      std::vector<std::string> ListLines=ReadFullFile(eFileLog);
      
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
