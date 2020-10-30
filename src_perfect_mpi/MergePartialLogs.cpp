#include "PerfectMPI_types.h"
#include "NumberTheory.h"
#include "Namelist.h"
#include "MatrixCanonicalForm.h"
#include "Temp_PerfectForm.h"
#include <unordered_map>

using TypeIndexRed = std::pair<int,int>;



namespace std {
  template <>
  struct hash<TypeIndexRed>
  {
    std::size_t operator()(const TypeIndexRed& k) const
    {
      using std::size_t;
      using std::hash;
      using std::string;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:

      return ((hash<int>()(k.first)
               ^ (hash<int>()(k.second) << 1)) >> 1);
    }
  };

}


template<typename T>
struct InformationMatrix {
  TypePerfectExch<T> ePerfect;
  int nbAdjacent;
  std::vector<int> ListIAdj;
};


int main(int argc, char* argv[])
{
  //  using T=mpq_class;
  using Tint=int;
  //
  try {
    if (argc != 3) {
      std::cerr << "MargePartialLogs [PrefixLog] [FileOut]\n";
      throw TerminalException{1};
    }
    std::string PrefixLog = argv[1];
    std::string FileOut   = argv[2];
    MyMatrix<Tint> ZerMat(0,0);
    InformationMatrix<Tint> eInfo { {0, ZerMat}, -1, {}};
    //
    std::unordered_map<TypeIndexRed,InformationMatrix<Tint>> ListInfoMatrices;
    int iProc=0;
    while(true) {
      std::stringstream s;
      s << PrefixLog << iProc;
      std::string eFileLog(s.str());
      std::cerr << "eFileLog=" << eFileLog << "\n";
      if (!(IsExistingFile(eFileLog))) {
        break;
      }
      std::cerr << "File exists, processing it\n";
      //
      std::vector<std::string> ListLines=ReadFullFile(eFileLog);
      std::vector<std::string> LStrParse;
      int nbLine=ListLines.size();
      for (int iLine=0; iLine<nbLine; iLine++) {
        std::cerr << "iLine=" << iLine << " / " << nbLine << "\n";
        std::string eLine = ListLines[iLine];
        // Looking for number of adjacent matrices
        LStrParse = STRING_ParseSingleLine(eLine, {"Number of Adjacent for idxMatrixF=", " nbAdjacent=", " END"});
        if (LStrParse.size() > 0) {
          int idxMatrixF=std::stoi(LStrParse[0]);
          int nbAdjacent=std::stoi(LStrParse[1]);
          TypeIndexRed eTyp{iProc, idxMatrixF};
          auto iter=ListInfoMatrices.find(eTyp);
          if (iter == ListInfoMatrices.end()) {
            ListInfoMatrices[eTyp] = eInfo;
          }
          ListInfoMatrices[eTyp].nbAdjacent = nbAdjacent;
        }
        //
        LStrParse = STRING_ParseSingleLine(eLine, {"Inserting New perfect form", " idxMatrixCurrent=", " Obtained from ", "END"});
        if (LStrParse.size() > 0) {
          TypePerfectExch<Tint> ePerfect = ParseStringToPerfectExch<Tint>(LStrParse[0]);
          int idxMatrixCurrent=std::stoi(LStrParse[1]);
          TypeIndex eIndex = ParseStringToTypeIndex(LStrParse[2]);
          TypeIndexRed eIndexRed1{eIndex.iProc, eIndex.idxMatrix};
          auto iter1=ListInfoMatrices.find(eIndexRed1);
          if (iter1 == ListInfoMatrices.end()) {
            ListInfoMatrices[eIndexRed1] = {};
          }
          ListInfoMatrices[eIndexRed1].ListIAdj.push_back(eIndex.iAdj);
          //
          TypeIndexRed eIndexRed2{iProc, idxMatrixCurrent};
          auto iter2=ListInfoMatrices.find(eIndexRed2);
          if (iter2 == ListInfoMatrices.end()) {
            ListInfoMatrices[eIndexRed2] = {};
          }
          ListInfoMatrices[eIndexRed2].ePerfect = ePerfect;
        }
        //
        LStrParse = STRING_ParseSingleLine(eLine, {"Reading existing matrix=", " idxMatrixCurrent=", "END"});
        if (LStrParse.size() > 0) {
          TypePerfectExch<Tint> ePerfect = ParseStringToPerfectExch<Tint>(LStrParse[0]);
          int idxMatrix=std::stoi(LStrParse[1]);
          TypeIndexRed eIndexRed{iProc, idxMatrix};
          auto iter=ListInfoMatrices.find(eIndexRed);
          if (iter == ListInfoMatrices.end())
            ListInfoMatrices[eIndexRed] = {};
          ListInfoMatrices[eIndexRed].ePerfect = ePerfect;
        }
        //
        LStrParse = STRING_ParseSingleLine(eLine, {"Processed entry=", "END"});
        if (LStrParse.size() > 0) {
          TypeIndex eIndex = ParseStringToTypeIndex(LStrParse[0]);
          TypeIndexRed eIndexRed{eIndex.iProc, eIndex.idxMatrix};
          auto iter=ListInfoMatrices.find(eIndexRed);
          if (iter == ListInfoMatrices.end())
            ListInfoMatrices[eIndexRed] = {};
          ListInfoMatrices[eIndexRed].ListIAdj.push_back(eIndex.iAdj);
        }
      }
      iProc++;
    }
    //
    // Now writing down the file
    //
    std::ofstream os(FileOut);
    os << ListInfoMatrices.size();
    int SumStatus=0;
    for (auto & ePair : ListInfoMatrices) {
      int eStatus=0;
      if ((ePair.second.nbAdjacent == int(ePair.second.ListIAdj.size())) && ePair.second.nbAdjacent > 0)
        eStatus=1;
      SumStatus += eStatus;
      os << eStatus << "\n";
      if (ePair.second.ePerfect.incd == 0) {
        std::cerr << "The incidence is zero. This is not allowed\n";
        throw TerminalException{1};
      }
      os << ePair.second.ePerfect.incd << "\n";
      WriteMatrix(os, ePair.second.ePerfect.eMat);
    }
    std::cerr << "SumStatus=" << SumStatus << " |ListInfoMatrices|=" << ListInfoMatrices.size() << "\n";
    std::cerr << "Normal termination of the program\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
