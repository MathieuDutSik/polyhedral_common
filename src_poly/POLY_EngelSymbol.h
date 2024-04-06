/// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_ENGELSYMBOL_H_
#define SRC_POLY_POLY_ENGELSYMBOL_H_

struct EngelPolyhedralSubordination {
  int n;
  std::vector<CollectedResult<int>> TheSub;
  std::vector<vectface> ListListFace;
};

template <typename T>
EngelPolyhedralSubordination
ComputeEngelPolyhedralSubordination(MyMatrix<T> const &EXT,
                                    MyMatrix<T> const &FAC) {
  int nbFac = FAC.rows();
  int nbExt = EXT.rows();
  vectface FACset(nbExt);
  int n = FAC.cols();
  std::cerr << "nbFac=" << nbFac << " nbExt=" << nbExt << " n=" << n << "\n";
  MyVector<T> eFac;
  MyVector<T> eExt;
  for (int iFac = 0; iFac < nbFac; iFac++) {
    eFac = FAC.row(iFac);
    Face eFace(nbExt);
    //    std::cerr << "iFac=" << iFac << "\n";
    for (int iExt = 0; iExt < nbExt; iExt++) {
      eExt = EXT.row(iExt);
      T eScal = eFac.dot(eExt);
      int eValIns;
      if (eScal == 0)
        eValIns = 1;
      else
        eValIns = 0;
      eFace[iExt] = eValIns;
      //      std::cerr << "  iExt=" << iExt << " val=" << eValIns << "\n";
    }
    FACset.push_back(eFace);
  }
  std::vector<vectface> ListListFace;
  ListListFace.emplace_back(std::move(FACset));
  std::vector<CollectedResult<int>> TheSub;
  for (int eDim = 0; eDim < n - 1; eDim++) {
    int TheRank = n - 2 - eDim;
    std::cerr << "eDim=" << eDim << " TheRank=" << TheRank << "\n";
    std::vector<int> ListSizes;
    std::unordered_set<Face> NewListFace_set;
    std::cerr << "  siz=" << ListListFace[eDim].size() << "\n";
    for (auto &eFace : ListListFace[eDim]) {
      int nb = eFace.count();
      std::vector<int> eList(nb);
      boost::dynamic_bitset<>::size_type aRow = eFace.find_first();
      for (int i = 0; i < nb; i++) {
        eList[i] = static_cast<int>(aRow);
        aRow = eFace.find_next(aRow);
      }
      std::unordered_set<Face> ListSubFace;
      for (auto &fFace : FACset) {
        Face gFace(nbExt);
        std::vector<int> gList;
        int eIncd = 0;
        for (auto &eVal : eList) {
          if (fFace[eVal] == 1) {
            gList.push_back(eVal);
            gFace[eVal] = 1;
            eIncd++;
          }
        }
        bool IsFace;
        if (eIncd < TheRank) {
          IsFace = false;
        } else {
          MyMatrix<T> EXTmat = SelectRow(EXT, gList);
          int rank = RankMat(EXTmat);
          IsFace = rank == TheRank;
        }
        if (IsFace)
          ListSubFace.insert(gFace);
      }
      int eSize = ListSubFace.size();
      for (auto &rgFace : ListSubFace)
        NewListFace_set.insert(rgFace);
      ListSizes.push_back(eSize);
    }
    TheSub.push_back(Collected(ListSizes));
    vectface NewListFace_vect(nbExt);
    for (auto &rFace : NewListFace_set)
      NewListFace_vect.push_back(rFace);
    ListListFace.emplace_back(std::move(NewListFace_vect));
  }
  return {n, std::move(TheSub), std::move(ListListFace)};
}

template <typename T>
void ComputeFileFaceLatticeInfo(std::string const &eFile,
                                MyMatrix<T> const &EXT,
                                MyMatrix<T> const &FAC) {
  EngelPolyhedralSubordination eEngel =
      ComputeEngelPolyhedralSubordination(EXT, FAC);
  std::ofstream os(eFile);
  os << "return [";
  int len = eEngel.ListListFace.size();
  int nbExt = EXT.rows();
  for (int iDim = 0; iDim < len; iDim++) {
    if (iDim > 0)
      os << ",\n";
    int nbFace = eEngel.ListListFace[iDim].size();
    os << "[";
    for (int iFace = 0; iFace < nbFace; iFace++) {
      if (iFace > 0)
        os << ",";
      Face eFace = eEngel.ListListFace[iDim][iFace];
      bool IsFirst = true;
      os << "[";
      for (int iExt = 0; iExt < nbExt; iExt++) {
        if (eFace[iExt] == 1) {
          if (!IsFirst)
            os << ",";
          IsFirst = false;
          int eVal = iExt + 1;
          os << eVal;
        }
      }
      os << "]";
    }
    os << "]";
  }
  os << "];\n";
}

template <typename T>
void ComputeEngelPolyhedralSubordinationFile(std::string const &eFile,
                                             MyMatrix<T> const &EXT,
                                             MyMatrix<T> const &FAC) {
  EngelPolyhedralSubordination eEngel =
      ComputeEngelPolyhedralSubordination(EXT, FAC);
  std::ofstream os(eFile);
  os << "return [";
  int len = eEngel.TheSub.size();
  for (int i = 0; i < len; i++) {
    if (i > 0)
      os << ",\n";
    CollectedResult<int> eColl = eEngel.TheSub[i];
    int nbEnt = eColl.LVal.size();
    os << "[";
    for (int iEnt = 0; iEnt < nbEnt; iEnt++) {
      int eVal = eColl.LVal[iEnt];
      int eMult = eColl.LMult[iEnt];
      if (iEnt > 0)
        os << ",";
      os << "[" << eVal << "," << eMult << "]";
    }
    os << "]";
  }
  os << "];\n";
}

// clang-format off
#endif  // SRC_POLY_POLY_ENGELSYMBOL_H_
// clang-format on
