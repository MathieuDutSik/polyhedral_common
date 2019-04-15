#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PerfectForm.h"


struct TypeIndex {
  int iProc;
  int idxMatrix;
};

template<typename T>
struct InformationMatrix {
  MyMatrix<T> eMat;
  int nbAdjacent;
  std::vector<int> ListIDaj;
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
    std::unordered_map<TypeIndex,InformationMatrix<Tint>> ListInfoMatrices;
    int iProc=0;
    while(true) {
      std::stringstream s;
      s << PrefixLog + iProc;
      string eFileLog(s.str());
      //
      std::vector<std::string> ListLines=ReadFullFile(eFileLog);
      std::vector<string> LStrParse;
      for (int iLine=0; iLine<ListLines.size(); iLine++) {
        std::string eLine = ListLines[iLine];
        // Looking for number of adjacent matrices
        LStrParse = STRING_ParseSingleLine(eLine, {"Number of Adjacent for idxMatrixF=", " nbAdjacent=", " END"});
        if (LStrParse.size() > 0) {
          int idxMatrixF=StringToInt(LStrParse[0]);
          int nbAdjacent=StringToInt(LStrParse[1]);
          TypeIndex eTyp{iProc, idxMatrixF};
          auto iter=ListInfoMatrices.find(eTyp);
          if (iter == ListInfoMatrices.end()) {
            std::cerr << "Failed to find entry\n";
            throw TerminalException{1};
          }
          iter.second.nbAdjacent = nbAdjacent;
        }
        //
        LStrParse = STRING_ParseSingleLine(eLine, {"Inserting New perfect form", "Obtained from ", "END"});
        if (LStrParse.size() > 0) {
          
        }
      }
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
