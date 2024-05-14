// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "CtypeMPI_types.h"
#include "NumberTheory.h"
#include <unordered_map>

using TypeIndexRed = std::pair<int, int>;

namespace std {
template <> struct hash<TypeIndexRed> {
  std::size_t operator()(const TypeIndexRed &k) const {
    using std::hash;
    using std::size_t;
    using std::string;

    // Compute individual hash values for first,
    // second and third and combine them using XOR
    // and bit shifting:

    return ((hash<int>()(k.first) ^ (hash<int>()(k.second) << 1)) >> 1);
  }
};

// clang-format off
}  // namespace std
// clang-format on

template <typename T> struct InformationMatrix {
  TypeCtypeExch<T> eCtype;
  int nbAdjacent;
  std::vector<int> ListIAdj;
};

int main(int argc, char *argv[]) {
  //  using T=mpq_class;
  using Tint = int;
  //
  try {
    if (argc != 3) {
      std::cerr << "CTYP_MergeLogs [PrefixLog] [FileOut]\n";
      throw TerminalException{1};
    }
    std::string PrefixLog = argv[1];
    std::string FileOut = argv[2];
    //
    std::unordered_map<TypeIndexRed, InformationMatrix<Tint>> ListInfoMatrices;
    int iProc = 0;
    while (true) {
      std::string eFileLog = PrefixLog + std::to_string(iProc);
      std::cerr << "eFileLog=" << eFileLog << "\n";
      if (!(IsExistingFile(eFileLog))) {
        break;
      }
      std::cerr << "File exists, processing it\n";
      //
      std::vector<std::string> ListLines = ReadFullFile(eFileLog);
      std::vector<std::string> LStrParse;
      int nbLine = ListLines.size();
      for (int iLine = 0; iLine < nbLine; iLine++) {
        std::cerr << "iLine=" << iLine << " / " << nbLine << "\n";
        std::string eLine = ListLines[iLine];
        // Looking for number of adjacent matrices
        LStrParse = STRING_ParseSingleLine(
            eLine,
            {"Number of Adjacent for idxMatrixF=", " nbAdjacent=", " END"});
        if (LStrParse.size() > 0) {
          int idxMatrixF = std::stoi(LStrParse[0]);
          int nbAdjacent = std::stoi(LStrParse[1]);
          TypeIndexRed eTyp{iProc, idxMatrixF};
          InformationMatrix<Tint> &eData = ListInfoMatrices[eTyp];
          eData.nbAdjacent = nbAdjacent;
        }
        //
        LStrParse = STRING_ParseSingleLine(
            eLine, {"Inserting New ctype",
                    " idxMatrixCurrent=", " Obtained from ", "END"});
        if (LStrParse.size() > 0) {
          TypeCtypeExch<Tint> eCtype =
              ParseStringToCtypeExch<Tint>(LStrParse[0]);
          int idxMatrixCurrent = std::stoi(LStrParse[1]);
          TypeIndex eIndex = ParseStringToTypeIndex(LStrParse[2]);
          TypeIndexRed eIndexRed1{eIndex.iProc, eIndex.idxMatrix};
          InformationMatrix<Tint> &eData1 = ListInfoMatrices[eIndexRed1];
          eData1.ListIAdj.push_back(eIndex.iAdj);
          //
          TypeIndexRed eIndexRed2{iProc, idxMatrixCurrent};
          InformationMatrix<Tint> &eData2 = ListInfoMatrices[eIndexRed2];
          eData2.eCtype = eCtype;
        }
        //
        LStrParse = STRING_ParseSingleLine(
            eLine, {"Reading existing matrix=", " idxMatrixCurrent=", "END"});
        if (LStrParse.size() > 0) {
          TypeCtypeExch<Tint> eCtype =
              ParseStringToCtypeExch<Tint>(LStrParse[0]);
          int idxMatrix = std::stoi(LStrParse[1]);
          TypeIndexRed eIndexRed{iProc, idxMatrix};
          InformationMatrix<Tint> &eData = ListInfoMatrices[eIndexRed];
          eData.eCtype = eCtype;
        }
        //
        LStrParse = STRING_ParseSingleLine(eLine, {"Processed entry=", "END"});
        if (LStrParse.size() > 0) {
          TypeIndex eIndex = ParseStringToTypeIndex(LStrParse[0]);
          TypeIndexRed eIndexRed{eIndex.iProc, eIndex.idxMatrix};
          InformationMatrix<Tint> &eData = ListInfoMatrices[eIndexRed];
          eData.ListIAdj.push_back(eIndex.iAdj);
        }
      }
      iProc++;
    }
    //
    // Now writing down the file
    //
    std::ofstream os(FileOut);
    os << ListInfoMatrices.size();
    int iCtype = 0;
    int nbDone = 0;
    int nbNotDone = 0;
    Tint MaxCoeff = 0;
    for (auto &ePair : ListInfoMatrices) {
      int eStatus = 0;
      if ((ePair.second.nbAdjacent ==
           static_cast<int>(ePair.second.ListIAdj.size())) &&
          ePair.second.nbAdjacent > 0)
        eStatus = 1;
      os << eStatus << "\n";
      WriteMatrix(os, ePair.second.eCtype.eMat);
      int nbRow = ePair.second.eCtype.eMat.rows();
      int nbCol = ePair.second.eCtype.eMat.cols();
      for (int iRow = 0; iRow < nbRow; iRow++) {
        for (int iCol = 0; iCol < nbCol; iCol++) {
          Tint val = T_abs(ePair.second.eCtype.eMat(iRow, iCol));
          MaxCoeff = std::max(MaxCoeff, val);
        }
      }
      std::cerr << "iCtype=" << iCtype
                << " nbAdjacent=" << ePair.second.nbAdjacent
                << " status=" << eStatus << "\n";
      if (eStatus == 1)
        nbDone++;
      else
        nbNotDone++;
      iCtype++;
    }
    int TotalForm = nbDone + nbNotDone;
    std::cerr << "MaxCoeff = " << MaxCoeff << "\n";
    std::cerr << "nbDone=" << nbDone << " nbNotDone=" << nbNotDone
              << " TotalForm=" << TotalForm << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
