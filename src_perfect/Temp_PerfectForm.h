// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_TEMP_PERFECTFORM_H_
#define SRC_PERFECT_TEMP_PERFECTFORM_H_

// clang-format off>
#include "MatrixGroup.h"
#include "Parallel_Classes_Types.h"
#include "PolytopeEquiStab.h"
#include "Positivity.h"
#include "Tspace_General.h"
#include <map>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_PERFECT_FORM
#define DEBUG_FLIP
#endif

template <typename T, typename Tint> struct NakedPerfect {
  MyMatrix<T> eGram;
  MyMatrix<Tint> SHV;
  MyMatrix<Tint> SHVred;
  MyMatrix<T> PerfDomEXT;
  std::vector<std::vector<int>> ListBlock;
  std::vector<int> ListPos;
};

template <typename T, typename Tint>
MyMatrix<T> GetNakedPerfectConeClassical(MyMatrix<Tint> const &M) {
  int nbRow = M.rows();
  int n = M.cols();
  int dimSymm = n * (n + 1) / 2;
  MyMatrix<T> RetMat(nbRow, dimSymm);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<Tint> V = GetMatrixRow(M, iRow);
    MyMatrix<T> M(n, n);
    for (int u = 0; u < n; u++)
      for (int v = 0; v < n; v++)
        M(u, v) = V(u) * V(v);
    MyVector<T> Vm = SymmetricMatrixToVector(M);
    AssignMatrixRow(RetMat, iRow, Vm);
  }
  return RetMat;
}

template <typename T, typename Tint>
NakedPerfect<T, Tint> GetNakedPerfectCone(LinSpaceMatrix<T> const &LinSpa,
                                          MyMatrix<T> const &eGram,
                                          Tshortest<T, Tint> const &RecSHV,
                                          std::ostream& os) {
  int nbSHV = RecSHV.SHV.rows();
  std::vector<int> ListPos(nbSHV);
  int nbMat = LinSpa.ListMat.size();
  int n = eGram.rows();
  MyMatrix<T> RyshkovLoc(nbSHV, nbMat);
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<Tint> eVect = GetMatrixRow(RecSHV.SHV, iSHV);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      T eSum = EvaluationQuadForm<T, Tint>(LinSpa.ListMat[iMat], eVect);
      RyshkovLoc(iSHV, iMat) = eSum;
    }
  }
  //  os << "RyshkovLoc=\n";
  //  WriteMatrix(os, RyshkovLoc);
  int nbBlock = 0;
  std::vector<std::vector<int>> ListBlock;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> eVec1 = GetMatrixRow(RyshkovLoc, iSHV);
    bool IsMatch = false;
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      if (!IsMatch) {
        int jSHV = ListBlock[iBlock][0];
        MyVector<T> eVec2 = GetMatrixRow(RyshkovLoc, jSHV);
        bool test = IsVectorPositiveMultiple(eVec1, eVec2);
        if (test) {
          IsMatch = true;
          ListBlock[iBlock].push_back(iSHV);
          ListPos[iSHV] = iBlock;
        }
      }
    }
    if (!IsMatch) {
      ListBlock.push_back({iSHV});
      ListPos[iSHV] = nbBlock;
      nbBlock++;
    }
  }
#ifdef DEBUG_PERFECT_FORM
  os << "m=" << n << " nbBlock=" << nbBlock << " nbSHV=" << nbSHV
     << " nbMat=" << nbMat << "\n";
#endif
  MyMatrix<T> PerfDomEXT(nbBlock, nbMat);
  MyMatrix<Tint> SHVred(nbBlock, n);
  for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
    int iSHV = ListBlock[iBlock][0];
    MyVector<Tint> eVect = GetMatrixRow(RecSHV.SHV, iSHV);
    AssignMatrixRow(SHVred, iBlock, eVect);
    for (int iMat = 0; iMat < nbMat; iMat++)
      PerfDomEXT(iBlock, iMat) = RyshkovLoc(iSHV, iMat);
  }
  return {eGram, RecSHV.SHV, SHVred, PerfDomEXT, ListBlock, ListPos};
}


template<typename T, typename Tgroup>
struct RyshkovGRP {
  MyMatrix<T> PerfDomEXT;
  Tgroup GRPsub;
};




template <typename T, typename Tgroup>
RyshkovGRP<T,Tgroup> GetNakedPerfectCone_GRP(LinSpaceMatrix<T> const &LinSpa,
                                             MyMatrix<T> const &eGram,
                                             MyMatrix<T> const& SHV_T,
                                             Tgroup const& GRP,
                                             std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int nbSHV = SHV_T.rows();
  std::vector<int> ListPos(nbSHV);
  int nbMat = LinSpa.ListMat.size();
  int n = eGram.rows();
  MyMatrix<T> RyshkovLoc(nbSHV, nbMat);
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> eVect = GetMatrixRow(SHV_T, iSHV);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      T eSum = EvaluationQuadForm<T, T>(LinSpa.ListMat[iMat], eVect);
      RyshkovLoc(iSHV, iMat) = eSum;
    }
  }
  //  os << "RyshkovLoc=\n";
  //  WriteMatrix(os, RyshkovLoc);
  int nbBlock = 0;
  std::vector<std::vector<int>> ListBlock;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> eVec1 = GetMatrixRow(RyshkovLoc, iSHV);
    bool IsMatch = false;
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      if (!IsMatch) {
        int jSHV = ListBlock[iBlock][0];
        MyVector<T> eVec2 = GetMatrixRow(RyshkovLoc, jSHV);
        bool test = IsVectorPositiveMultiple(eVec1, eVec2);
        if (test) {
          IsMatch = true;
          ListBlock[iBlock].push_back(iSHV);
          ListPos[iSHV] = iBlock;
        }
      }
    }
    if (!IsMatch) {
      ListBlock.push_back({iSHV});
      ListPos[iSHV] = nbBlock;
      nbBlock++;
    }
  }
#ifdef DEBUG_PERFECT_FORM
  os << "m=" << n << " nbBlock=" << nbBlock << " nbSHV=" << nbSHV
     << " nbMat=" << nbMat << "\n";
#endif
  MyMatrix<T> PerfDomEXT(nbBlock, nbMat);
  MyMatrix<T> SHVred(nbBlock, n);
  for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
    int iSHV = ListBlock[iBlock][0];
    MyVector<T> eVect = GetMatrixRow(SHV_T, iSHV);
    AssignMatrixRow(SHVred, iBlock, eVect);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      PerfDomEXT(iBlock, iMat) = RyshkovLoc(iSHV, iMat);
    }
  }
  std::vector<Telt> l_gens;
  for (auto & eGen: GRP.GeneratorsOfGroup()) {
    std::vector<Tidx> eList(nbBlock);
    for (int iBlock=0; iBlock<nbBlock; iBlock++) {
      Tidx iSHV = ListBlock[iBlock][0];
      Tidx jSHV = eGen.at(iSHV);
      Tidx jBlock = ListPos[jSHV];
      eList[iBlock] = jBlock;
    }
    Telt eElt(eList);
    l_gens.push_back(eElt);
  }
  Tgroup GRPsub(l_gens, nbBlock);
  return {PerfDomEXT, GRPsub};
}



struct PerfEquivInfo {
  int iOrbit;
  MyMatrix<int> eMatEquiv;
  Face eInc;
};

template <typename T, typename Tint, typename Tgroup>
Tgroup MapLatticeGroupToConeGroup(NakedPerfect<T, Tint> const &eNaked,
                                  Tgroup const &GRPshv) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int nbBlock = eNaked.ListBlock.size();
  std::vector<Telt> ListGen;
  std::vector<Telt> LGen = GRPshv.GeneratorsOfGroup();
  for (auto &eGen : LGen) {
    std::vector<Tidx> v(nbBlock);
    for (int iBlock = 0; iBlock < nbBlock; iBlock++) {
      int iSHV = eNaked.ListBlock[iBlock][0];
      int jSHV = OnPoints(iSHV, eGen);
      int jBlock = eNaked.ListPos[jSHV];
      v[iBlock] = jBlock;
    }
    ListGen.push_back(Telt(v));
  }
  return Tgroup(ListGen, nbBlock);
}

struct QueryEquivInfo {
  bool result;
  int nbOrbit;
  PerfEquivInfo eEquiv;
};

template <typename Tint>
bool TestInclusionSHV(MyMatrix<Tint> const &TheSHVbig,
                      MyMatrix<Tint> const &TheSHVsma) {
  int nbRowSma = TheSHVsma.rows();
  int nbRowBig = TheSHVbig.rows();
  int n = TheSHVsma.cols();
  for (int iRowSma = 0; iRowSma < nbRowSma; iRowSma++) {
    bool WeMatch = false;
    for (int iRowBig = 0; iRowBig < nbRowBig; iRowBig++)
      if (!WeMatch) {
        int SumErr = 0;
        for (int i = 0; i < n; i++) {
          Tint eVal = TheSHVbig(iRowBig, i);
          Tint fVal = TheSHVsma(iRowSma, i);
          if (eVal != fVal)
            SumErr++;
        }
        if (SumErr == 0)
          WeMatch = true;
      }
    if (!WeMatch)
      return false;
  }
  return true;
}

template <typename T, typename Tint, typename Fadmissible, typename Fshortest>
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
Kernel_Flipping_Perfect(Fadmissible f_admissible,
                        Fshortest f_shortest,
                        MyMatrix<T> const &eMatIn, MyMatrix<T> const &eMatDir,
                        [[maybe_unused]] std::ostream & os) {
  std::vector<MyMatrix<T>> ListMat;
  std::vector<Tshortest<T, Tint>> ListShort;
  // Memoization procedure
  auto RetrieveShortestDesc =
      [&](MyMatrix<T> const &eMat) -> Tshortest<T, Tint> {
    int len = ListMat.size();
    for (int i = 0; i < len; i++) {
      if (ListMat[i] == eMat) {
#ifdef DEBUG_FLIP
        os << "PERF: Memoization success\n";
#endif
        return ListShort[i];
      }
    }
#ifdef DEBUG_FLIP
    os << "PERF: CALL ShortestFunction\n";
#endif
    Tshortest<T, Tint> RecSHV = f_shortest(eMat);
    ListMat.push_back(eMat);
    ListShort.push_back(RecSHV);
    return RecSHV;
  };
#ifdef DEBUG_FLIP
  os << "PERF: RetrieveShortestDesc, case 1 (eMatIn)\n";
#endif
  Tshortest<T, Tint> const RecSHVperf = RetrieveShortestDesc(eMatIn);
#ifdef DEBUG_FLIP
  os << "Kernel_Flipping_Perfect : SHVinformation=\n";
  int nbSHV = RecSHVperf.SHV.rows();
  int n = RecSHVperf.SHV.cols();
  int nbZero_sumMat = 0;
  std::vector<int> ListScal(nbSHV);
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<Tint> V = GetMatrixRow(RecSHVperf.SHV, iSHV);
    T sumMatIn = EvaluationQuadForm(eMatIn, V);
    T sumMatDir = EvaluationQuadForm(eMatDir, V);
    int eVal = 1;
    if (sumMatDir == 0) {
      nbZero_sumMat++;
      eVal = 0;
    }
    ListScal[iSHV] = eVal;
    os << "iSHV=" << iSHV << " V=";
    for (int i = 0; i < n; i++)
      os << V(i) << " ";
    os << "sumMatIn=" << sumMatIn << " sumMatDir=" << sumMatDir << "\n";
  }
  MyMatrix<Tint> SHVface(nbZero_sumMat, n);
  int idx = 0;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    if (ListScal[iSHV] == 0) {
      for (int i = 0; i < n; i++)
        SHVface(idx, i) = RecSHVperf.SHV(iSHV, i);
      idx++;
    }
  }
  os << "SHVface=\n";
  WriteMatrix(os, SHVface);
  os << "nbZero_sumMat=" << nbZero_sumMat << "\n";
  MyMatrix<T> ConeClassicalInput =
      GetNakedPerfectConeClassical<T>(RecSHVperf.SHV);
  os << "RankMat(ConeClassicalInput)=" << RankMat(ConeClassicalInput) << "\n";
  MyMatrix<T> ConeClassicalFace = GetNakedPerfectConeClassical<T>(SHVface);
  os << "RankMat(ConeClassicalFace)=" << RankMat(ConeClassicalFace) << "\n";
#endif
  T TheUpperBound = 1;
  T TheLowerBound = 0;
#ifdef DEBUG_FLIP
  int iterLoop = 0;
#endif
  while (true) {
    MyMatrix<T> Qupp = eMatIn + TheUpperBound * eMatDir;
#ifdef DEBUG_FLIP
    os << "PERF: CALL IsAdmissible (Qupp)\n";
#endif
    bool test = f_admissible(Qupp);
#ifdef DEBUG_FLIP
    iterLoop++;
    os << "iterLoop=" << iterLoop << " TheUpperBound=" << TheUpperBound
       << " TheLowerBound=" << TheLowerBound << " test=" << test << "\n";
#endif
    if (!test) {
      TheUpperBound = (TheUpperBound + TheLowerBound) / 2;
    } else {
#ifdef DEBUG_FLIP
      os << "PERF: Qupp=\n";
      WriteMatrix(os, Qupp);
      os << "PERF: RetrieveShortestDesc, case 2 (Qupp)\n";
#endif
      Tshortest<T, Tint> RecSHV = RetrieveShortestDesc(Qupp);
#ifdef DEBUG_FLIP
      os << "ITER: RecSHV.eMin=" << RecSHV.eMin << "\n";
      os << "ITER: RecSHV.SHV=\n";
      WriteMatrix(os, RecSHV.SHV);
#endif
      if (RecSHV.eMin == RecSHVperf.eMin) {
        T nLow = TheUpperBound;
        T nUpp = 2 * TheUpperBound;
        TheLowerBound = nLow;
        TheUpperBound = nUpp;
      } else {
        break;
      }
    }
  }
#ifdef DEBUG_FLIP
  os << "FIRST LOOP FINISHED TheUpperBound=" << TheUpperBound
     << " TheLowerBound=" << TheLowerBound << "\n";
#endif
  while (true) {
#ifdef DEBUG_FLIP
    os << "Now TheUpperBound=" << TheUpperBound
       << " TheLowerBound=" << TheLowerBound << "\n";
#endif
    MyMatrix<T> Qlow = eMatIn + TheLowerBound * eMatDir;
    MyMatrix<T> Qupp = eMatIn + TheUpperBound * eMatDir;
#ifdef DEBUG_FLIP
    os << "PERF: RetrieveShortestDesc, case 3 (Qlow)\n";
#endif
    Tshortest<T, Tint> RecSHVlow = RetrieveShortestDesc(Qlow);
#ifdef DEBUG_FLIP
    os << "PERF: RetrieveShortestDesc, case 4 (Qupp)\n";
#endif
    Tshortest<T, Tint> RecSHVupp = RetrieveShortestDesc(Qupp);
#ifdef DEBUG_FLIP
    os << "RecSHVupp.eMin=" << RecSHVupp.eMin
       << " RecSHVlow.eMin=" << RecSHVlow.eMin << "\n";
    MyMatrix<T> ConeClassicalLow =
        GetNakedPerfectConeClassical<T>(RecSHVlow.SHV);
    os << "RankMat(ConeClassicalLow)=" << RankMat(ConeClassicalLow) << "\n";
    MyMatrix<T> ConeClassicalUpp =
        GetNakedPerfectConeClassical<T>(RecSHVupp.SHV);
    os << "RankMat(ConeClassicalUpp)=" << RankMat(ConeClassicalUpp) << "\n";
#endif
    bool test1 = RecSHVupp.eMin == RecSHVperf.eMin;
    bool test2 = TestInclusionSHV(RecSHVperf.SHV, RecSHVlow.SHV);
    if (test1) {
#ifdef DEBUG_FLIP
      os << "Return Qupp\n";
#endif
      return {std::move(Qupp), std::move(RecSHVupp)};
    }
    if (!test2) {
#ifdef DEBUG_FLIP
      os << "Qperf=\n";
      WriteMatrix(os, eMatIn);
      os << "Qlow=\n";
      WriteMatrix(os, Qlow);
      //
      os << "RecSHVperf.SHV=\n";
      WriteMatrix(os, RecSHVperf.SHV);
      os << "RecSHVlow.SHV=\n";
      WriteMatrix(os, RecSHVlow.SHV);
      os << "Return Qlow\n";
#endif
      return {std::move(Qlow), std::move(RecSHVlow)};
    }
    T TheGamma = (TheLowerBound + TheUpperBound) / 2;
    MyMatrix<T> Qgamma = eMatIn + TheGamma * eMatDir;
#ifdef DEBUG_FLIP
    os << "Qgamma=\n";
    WriteMatrix(os, Qgamma);
#endif
#ifdef DEBUG_FLIP
    os << "PERF: RetrieveShortestDesc, case 5 (Qgamma)\n";
#endif
    Tshortest<T, Tint> RecSHVgamma = RetrieveShortestDesc(Qgamma);
#ifdef DEBUG_FLIP
    os << "|RecSHVgamma.SHV|=" << RecSHVgamma.SHV.rows() << "\n";
    WriteMatrix(os, RecSHVgamma.SHV);
    MyMatrix<T> ConeClassicalGamma =
        GetNakedPerfectConeClassical<T>(RecSHVgamma.SHV);
    os << "RankMat(ConeClassicalGamma)=" << RankMat(ConeClassicalGamma) << "\n";
#endif
    if (RecSHVgamma.eMin >= RecSHVperf.eMin) {
      TheLowerBound = TheGamma;
    } else {
#ifdef DEBUG_FLIP
      os << "Assigning TheUpperBound to TheGamma=" << TheGamma << "\n";
#endif
      TheUpperBound = TheGamma;
      int nbRowGamma = RecSHVgamma.SHV.rows();
      for (int iRowGamma = 0; iRowGamma < nbRowGamma; iRowGamma++) {
        MyVector<Tint> eVectShort = GetMatrixRow(RecSHVgamma.SHV, iRowGamma);
#ifdef DEBUG_FLIP
        os << "iRowGamma=" << iRowGamma << " / " << nbRowGamma
           << " eVectShort=";
        WriteVectorNoDim(os, eVectShort);
#endif
        T rVal = EvaluationQuadForm<T, Tint>(eMatDir, eVectShort);
        T qVal = EvaluationQuadForm<T, Tint>(eMatIn, eVectShort);
        if (rVal < 0) {
          T TheVal = (RecSHVperf.eMin - qVal) / rVal;
          if (TheVal < TheUpperBound) {
#ifdef DEBUG_FLIP
            os << "iRowGamma=" << iRowGamma
               << " Assigning TheUpperBound to TheVal=" << TheVal << "\n";
            os << "rVal=" << rVal << " qVal=" << qVal
               << " RecSHVperf.eMin=" << RecSHVperf.eMin << "\n";
#endif
            TheUpperBound = TheVal;
          }
        }
      }
    }
  }
}

template <typename T, typename Tint>
size_t ComputeInvariantPerfectTspace(size_t const &seed,
                                     MyMatrix<T> const &eGram,
                                     Tshortest<T, Tint> const &RecSHV,
                                     std::ostream &os) {
#ifdef DEBUG_PERFECT_TSPACE
  os << "PERF_TSPACE: ComputeInvariantPerfectTspace, begin\n";
#endif
  using Tidx_value = int16_t;
  int n = eGram.rows();
  MyMatrix<T> eG = eGram / RecSHV.eMin;
  size_t n_row = RecSHV.SHV.rows();
  int nbSHV = 2 * n_row;

  MyVector<T> V(n);
  int iSign, jSign;
  auto f1 = [&](size_t iRow) -> void {
    size_t iRowRed = iRow % n_row;
    iSign = 1;
    if (iRow >= n_row) {
      iSign = -1;
    }
    for (int i=0; i<n; i++) {
      T eSum(0);
      for (int j=0; j<n; j++)
        eSum += eG(i, j) * RecSHV.SHV(iRowRed, i);
      V(i) = eSum;
    }
  };
  auto f2 = [&](size_t jRow) -> T {
    size_t jRowRed = jRow % n_row;
    jSign = 1;
    if (jRow >= n_row) {
      jSign = -1;
    }
    T eSum(0);
    for (int i=0; i<n; i++)
      eSum += V(i) * RecSHV.SHV(jRowRed, i);
    T fSum = eSum * (iSign * jSign);
    return fSum;
  };
  WeightMatrix<true, T, Tidx_value> WMat(nbSHV, f1, f2, os);
#ifdef DEBUG_PERFECT_TSPACE
  os << "PERF_TSPACE: ComputeInvariantPerfectTspace, We have WMat\n";
#endif
  size_t hash = GetInvariantWeightMatrix(seed, WMat);
#ifdef DEBUG_PERFECT_TSPACE
  os << "PERF_TSPACE: ComputeInvariantPerfectTspace, We have hash\n";
#endif
  return hash;
}

template <typename T, typename Tint>
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
Flipping_Perfect(MyMatrix<T> const &eMatIn, MyMatrix<T> const &eMatDir,
                 std::ostream &os) {
  auto f_admissible=[&os](MyMatrix<T> const &eMat) -> bool {
    return IsPositiveDefinite<T>(eMat, os);
  };
  auto f_shortest=[&os](MyMatrix<T> const &eMat) -> Tshortest<T, Tint> {
    return T_ShortestVector<T, Tint>(eMat, os);
  };
  return Kernel_Flipping_Perfect<T,Tint,decltype(f_admissible),decltype(f_shortest)>(f_admissible, f_shortest, eMatIn, eMatDir, os);
}

template <typename T, typename Tint>
std::pair<MyMatrix<T>, Tshortest<T, Tint>> GetOnePerfectForm(LinSpaceMatrix<T> const &LinSpa,
                              std::ostream &os) {
  int nbMat = LinSpa.ListMat.size();
  MyMatrix<T> ThePerfMat = LinSpa.SuperMat;
  while (true) {
    Tshortest<T, Tint> RecSHV = T_ShortestVector<T, Tint>(ThePerfMat, os);
    int nbShort = RecSHV.SHV.rows();
    MyMatrix<T> ScalMat(nbShort, nbMat);
    for (int iShort = 0; iShort < nbShort; iShort++) {
      MyVector<Tint> eVectShort = RecSHV.SHV.row(iShort);
      for (int iMat = 0; iMat < nbMat; iMat++) {
        T eNorm = EvaluationQuadForm<T, Tint>(LinSpa.ListMat[iMat], eVectShort);
        ScalMat(iShort, iMat) = eNorm;
      }
    }
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(ScalMat);
    int TheRank = eSelect.TheRank;
    if (TheRank == nbMat) {
      return {ThePerfMat, RecSHV};
    }
    MyVector<T> eVect = eSelect.NSP.row(0);
    MyMatrix<T> DirMat = LINSPA_GetMatrixInTspace(LinSpa, eVect);
    MyMatrix<T> eMatRet =
        Flipping_Perfect<T, Tint>(ThePerfMat, DirMat, os).first;
    ThePerfMat = eMatRet;
  }
}

template <typename T, typename Tint> struct SimplePerfect {
  MyMatrix<T> Gram;
};

template <typename T, typename Tint>
std::istream &operator>>(std::istream &is, SimplePerfect<T, Tint> &obj) {
  MyMatrix<T> eG = ReadMatrix<T>(is);
  obj = {eG};
  return is;
}

template <typename T, typename Tint>
std::ostream &operator<<(std::ostream &os, SimplePerfect<T, Tint> const &obj) {
  WriteMatrix(os, obj.Gram);
  return os;
}

template <typename T> struct SimplePerfectInv {
  size_t hash;
};

template <typename T>
bool operator==(SimplePerfectInv<T> const &x, SimplePerfectInv<T> const &y) {
  return x.hash == y.hash;
}

template <typename T>
std::istream &operator>>(std::istream &is, SimplePerfectInv<T> &obj) {
  size_t hash;
  is >> hash;
  obj = {hash};
  return is;
}
template <typename T>
std::ostream &operator<<(std::ostream &os, SimplePerfectInv<T> const &obj) {
  os << obj.hash;
  return os;
}

template <typename T>
bool operator<(SimplePerfectInv<T> const &x, SimplePerfectInv<T> const &y) {
  return x.hash < y.hash;
}

template <typename T, typename Tint>
struct invariant_info<SimplePerfect<T, Tint>> {
  typedef SimplePerfectInv<T> invariant_type;
};

template <typename T, typename Tint> struct equiv_info<SimplePerfect<T, Tint>> {
  typedef MyMatrix<Tint> equiv_type;
};

template <typename T, typename Tint>
SimplePerfectInv<T> ComputeInvariantSimplePerfect(MyMatrix<T> const &eGram,
                                                  std::ostream &os) {
  Tshortest<T, Tint> RecSHV = T_ShortestVector<T, Tint>(eGram, os);
  size_t seed = 1234;
  size_t hash = ComputeInvariantPerfectTspace<T,Tint>(seed, eGram, RecSHV, os);
  return {hash};
}

template<typename T, typename Tint>
MyMatrix<T> conversion_and_duplication(MyMatrix<Tint> const& SHV) {
  int dim = SHV.cols();
  int nbSHV = SHV.rows();
  MyMatrix<T> SHV_T(2 * nbSHV, dim);
  for (int i_row=0; i_row<nbSHV; i_row++) {
    for (int i=0; i<dim; i++) {
      T val = UniversalScalarConversion<T,Tint>(SHV(i_row,i));
      SHV_T(2*i_row  , i) = val;
      SHV_T(2*i_row+1, i) = -val;
    }
  }
  return SHV_T;
}


template <typename T, typename Tint>
MyMatrix<T> get_shv_t(MyMatrix<T> const& eMat, Tshortest<T, Tint> const &RecSHV, std::ostream& os) {
  MyMatrix<T> SHVorig_T = UniversalMatrixConversion<T, Tint>(RecSHV.SHV);
  if (IsFullDimZbasis(RecSHV.SHV)) {
    return conversion_and_duplication<T,Tint>(RecSHV.SHV);
  }
  MyMatrix<Tint> SHV = ExtractInvariantVectorFamilyZbasis<T, Tint>(eMat, os);
  return UniversalMatrixConversion<T, Tint>(SHV);
}


template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
SimplePerfect_TestEquivalence(LinSpaceMatrix<T> const &LinSpa,
                              MyMatrix<T> const &eMat1,
                              MyMatrix<T> const &eMat2,
                              Tshortest<T, Tint> const &RecSHV1,
                              Tshortest<T, Tint> const &RecSHV2,
                              std::ostream &os) {
  MyMatrix<T> SHV1_T = get_shv_t(eMat1, RecSHV1, os);
  MyMatrix<T> SHV2_T = get_shv_t(eMat2, RecSHV2, os);
  std::optional<MyMatrix<T>> opt = LINSPA_TestEquivalenceGramMatrix_SHV<T,Tgroup>(LinSpa, eMat1, eMat2, SHV1_T, SHV2_T, os);
  if (!opt) {
    return {};
  }
  MyMatrix<T> const& M_T = *opt;
  MyMatrix<Tint> M = UniversalMatrixConversion<Tint, T>(M_T);
  return M;
}

template <typename T, typename Tint, typename Tgroup>
Tgroup SimplePerfect_Stabilizer(LinSpaceMatrix<T> const &LinSpa,
                                MyMatrix<T> const &eMat,
                                Tshortest<T, Tint> const &RecSHV,
                                std::ostream& os) {
  //
  // Functionality for checking quality of equivalences
  //
  MyMatrix<T> SHVorig_T = conversion_and_duplication<T,Tint>(RecSHV.SHV);
  MyMatrix<T> SHV_T = get_shv_t(eMat, RecSHV, os);
  Result_ComputeStabilizer_SHV<T,Tgroup> result = LINSPA_ComputeStabilizer_SHV<T,Tgroup>(LinSpa, eMat, SHV_T, os);
  if (SHVorig_T == SHV_T) {
    // This is the most likely scenario: The original
    // SHVorig is adequate
    return result.get_perm_group(SHV_T);
  } else {
    std::vector<MyMatrix<T>> l_matr = result.get_list_matrix(SHV_T, eMat, LinSpa);
    return get_perm_group_from_list_matrices<T,Tgroup>(l_matr, SHVorig_T);
  }
}

// Test if a quadratic form is eutactic
// A quadratic form is eutactic if the inverse of eGram can be expressed
// as a positive linear combination of the outer products v ⊗ v of shortest vectors
template <typename T>
std::optional<MyVector<T>> IsEutactic(MyMatrix<T> const &eGram, MyMatrix<T> const& SHV_T, std::string const& eutacticity, std::ostream& os) {
#ifdef DEBUG_PERFECT_FORM
  os << "PERFECT: IsEutactic, beginning\n";
#endif
  int n = eGram.rows();
  if (eGram.cols() != n) {
    std::cerr << "Error: eGram must be square\n";
    return {};
  }
  if (SHV_T.cols() != n) {
    std::cerr << "Error: SHV_T must have n columns\n";
    return {};
  }

  int nbSHV = SHV_T.rows();
  if (nbSHV == 0) {
#ifdef DEBUG_PERFECT_FORM
    os << "PERFECT: No shortest vectors provided\n";
#endif
    return {};
  }

  // Dimension of the space of symmetric matrices
  int dimSymm = n * (n + 1) / 2;

#ifdef DEBUG_PERFECT_FORM
  os << "PERFECT: n=" << n << " nbSHV=" << nbSHV << " dimSymm=" << dimSymm << "\n";
#endif

  // Compute the inverse of eGram
  MyMatrix<T> eGramInv;
  eGramInv = Inverse(eGram);

  // Convert inverse to vector form
  MyVector<T> eGramInvVec = SymmetricMatrixToVector(eGramInv);

#ifdef DEBUG_PERFECT_FORM
  os << "PERFECT: eGram inverse computed and vectorized\n";
#endif

  // Create matrix of outer products v ⊗ v for each shortest vector v
  MyMatrix<T> OuterProductMat(nbSHV, dimSymm);

  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> v = GetMatrixRow(SHV_T, iSHV);

    // Compute outer product v ⊗ v as a symmetric matrix
    MyMatrix<T> OuterProd(n, n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        OuterProd(i, j) = v(i) * v(j);
      }
    }

    // Convert symmetric matrix to vector form
    MyVector<T> OuterProdVec = SymmetricMatrixToVector(OuterProd);
    AssignMatrixRow(OuterProductMat, iSHV, OuterProdVec);
  }

#ifdef DEBUG_PERFECT_FORM
  os << "PERFECT: OuterProductMat constructed, size: " << OuterProductMat.rows() << "x" << OuterProductMat.cols() << "\n";
#endif

  // Solve the linear system: OuterProductMat^T * coeffs = eGramInvVec
  // We want to find coefficients such that eGramInv = sum(coeffs[i] * v_i ⊗ v_i)
  MyMatrix<T> OuterProductMatT = TransposedMat(OuterProductMat);

  // Check if the system is solvable
  MyMatrix<T> AugmentedMat = ConcatenateMatVec(OuterProductMatT, eGramInvVec);

#ifdef DEBUG_PERFECT_FORM
  int rankA = RankMat(OuterProductMatT);
  int rankAug = RankMat(AugmentedMat);
  os << "PERFECT: rankA=" << rankA << " rankAug=" << rankAug << "\n";
#endif

  // Try to solve the system to get coefficients
  if (eutacticity == "Eutactic") {
    return SolutionMatStrictlyPositive<T>(OuterProductMatT, eGramInvVec, os);
  }
  if (eutacticity == "WeaklyEutactic") {
    return SolutionMatNonnegative<T>(OuterProductMatT, eGramInvVec, os);
  }
  std::cerr << "The input is not as expected\n";
  throw TerminalException{1};
}

// clang-format off
#endif  // SRC_PERFECT_TEMP_PERFECTFORM_H_
// clang-format on
