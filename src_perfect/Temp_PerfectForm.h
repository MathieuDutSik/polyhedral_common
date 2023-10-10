// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_TEMP_PERFECTFORM_H_
#define SRC_PERFECT_TEMP_PERFECTFORM_H_

#include "MatrixGroup.h"
#include "Parallel_Classes_Types.h"
#include "Temp_Positivity.h"
#include "Temp_Tspace_General.h"
#include <map>
#include <string>
#include <utility>
#include <vector>

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
                                          Tshortest<T, Tint> const &RecSHV) {
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
  //  std::cerr << "RyshkovLoc=\n";
  //  WriteMatrix(std::cerr, RyshkovLoc);
  int nbBlock = 0;
  std::vector<std::vector<int>> ListBlock;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<T> eVec1 = GetMatrixRow(RyshkovLoc, iSHV);
    bool IsMatch = false;
    for (int iBlock = 0; iBlock < nbBlock; iBlock++)
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
    if (!IsMatch) {
      ListBlock.push_back({iSHV});
      ListPos[iSHV] = nbBlock;
      nbBlock++;
    }
  }
  std::cerr << "nbBlock=" << nbBlock << "\n";
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

struct PerfEquivInfo {
  int iOrbit;
  MyMatrix<int> eMatEquiv;
  Face eInc;
};

template <typename T, typename Tint, typename Tgroup> struct SinglePerfect {
  MyMatrix<T> eGram;
  MyMatrix<Tint> SHV;
  MyMatrix<T> PerfDomEXT;
  Tgroup PerfDomGRP;
  std::vector<std::vector<int>> ListBlock;
  std::vector<PerfEquivInfo> ListEquivInfo;
  int eStatus;
  int eCons;
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

template <typename T, typename Tint, typename Tgroup>
SinglePerfect<T, Tint, Tgroup> GetPerfectCone(LinSpaceMatrix<T> const &LinSpa,
                                              MyMatrix<T> const &eGram,
                                              Tshortest<T, int> const &RecSHV) {
  NakedPerfect<T, Tint> eNaked = GetNakedPerfectCone(LinSpa, eGram, RecSHV);
  Tgroup TheGRPshv = PERF_Automorphism(LinSpa, eGram, RecSHV.SHV);
  Tgroup PerfDomGRP = MapLatticeGroupToConeGroup(eNaked, TheGRPshv);
  return {
      eGram, RecSHV.SHV, eNaked.PerfDomEXT, PerfDomGRP, eNaked.ListBlock, {},
      0,     0};
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
PERF_TestEquivalence(LinSpaceMatrix<T> const &LinSpa, MyMatrix<T> const &ePerf1,
                     MyMatrix<T> const &ePerf2, MyMatrix<Tint> const &SHV1,
                     MyMatrix<Tint> const &SHV2) {
  using Telt = typename Tgroup::Telt;
  using Tidx_value = int16_t;
  MyMatrix<T> T_SHV1 = UniversalMatrixConversion<T, Tint>(SHV1);
  MyMatrix<T> T_SHV2 = UniversalMatrixConversion<T, Tint>(SHV2);
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat1 =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV1, ePerf1,
                                                     LinSpa.ListComm);
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat2 =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV2, ePerf2,
                                                     LinSpa.ListComm);
  std::optional<Telt> eResEquiv =
      GetEquivalenceAsymmetricMatrix<std::vector<T>, T, Telt>(WMat1, WMat2);
  if (!eResEquiv) {
    return {};
  }
  MyMatrix<T> M3 =
      RepresentVertexPermutation(T_SHV1, T_SHV2, eResEquiv.TheEquiv);
  if (IsIntegralMatrix(M3)) {
    MyMatrix<Tint> eMatEquiv = UniversalMatrixConversion<Tint, T>(M3);
    return eMatEquiv;
  }
  std::cerr << "Need to write some code here\n";
  throw TerminalException{1};
}

struct QueryEquivInfo {
  bool result;
  int nbOrbit;
  PerfEquivInfo eEquiv;
};

template <typename T, typename Tint> struct RecShort {
  std::function<bool(MyMatrix<T> const &)> IsAdmissible;
  std::function<Tshortest<T, Tint>(MyMatrix<T> const &)> ShortestFunction;
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

template <typename T, typename Tint>
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
Kernel_Flipping_Perfect(RecShort<T, Tint> const &eRecShort,
                        MyMatrix<T> const &eMatIn, MyMatrix<T> const &eMatDir) {
  std::vector<MyMatrix<T>> ListMat;
  std::vector<Tshortest<T, Tint>> ListShort;
  // Memoization procedure
  auto RetriveShortestDesc =
      [&](MyMatrix<T> const &eMat) -> Tshortest<T, Tint> {
    int len = ListMat.size();
    for (int i = 0; i < len; i++)
      if (ListMat[i] == eMat)
        return ListShort[i];
    Tshortest<T, Tint> RecSHV = eRecShort.ShortestFunction(eMat);
    ListMat.push_back(eMat);
    ListShort.push_back(RecSHV);
    return RecSHV;
  };
  Tshortest<T, Tint> const RecSHVperf = RetriveShortestDesc(eMatIn);
#ifdef DEBUG_FLIP
  std::cerr << "Kernel_Flipping_Perfect : SHVinformation=\n";
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
    std::cerr << "iSHV=" << iSHV << " V=";
    for (int i = 0; i < n; i++)
      std::cerr << V(i) << " ";
    std::cerr << "sumMatIn=" << sumMatIn << " sumMatDir=" << sumMatDir << "\n";
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
  std::cerr << "SHVface=\n";
  WriteMatrix(std::cerr, SHVface);
  std::cerr << "nbZero_sumMat=" << nbZero_sumMat << "\n";
  MyMatrix<T> ConeClassicalInput =
      GetNakedPerfectConeClassical<T>(RecSHVperf.SHV);
  std::cerr << "RankMat(ConeClassicalInput)=" << RankMat(ConeClassicalInput)
            << "\n";
  MyMatrix<T> ConeClassicalFace = GetNakedPerfectConeClassical<T>(SHVface);
  std::cerr << "RankMat(ConeClassicalFace)=" << RankMat(ConeClassicalFace)
            << "\n";
#endif
  T TheUpperBound = 1;
  T TheLowerBound = 0;
#ifdef DEBUG_FLIP
  int iterLoop = 0;
#endif
  while (true) {
    MyMatrix<T> Qupp = eMatIn + TheUpperBound * eMatDir;
    bool test = eRecShort.IsAdmissible(Qupp);
#ifdef DEBUG_FLIP
    iterLoop++;
    std::cerr << "iterLoop=" << iterLoop << " TheUpperBound=" << TheUpperBound
              << " TheLowerBound=" << TheLowerBound << " test=" << test << "\n";
#endif
    if (!test) {
      TheUpperBound = (TheUpperBound + TheLowerBound) / 2;
    } else {
      Tshortest<T, Tint> RecSHV = RetriveShortestDesc(Qupp);
#ifdef DEBUG_FLIP
      std::cerr << "ITER: RecSHV.eMin=" << RecSHV.eMin << "\n";
      std::cerr << "ITER: RecSHV.SHV=\n";
      WriteMatrix(std::cerr, RecSHV.SHV);
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
  std::cerr << "FIRST LOOP FINISHED TheUpperBound=" << TheUpperBound
            << " TheLowerBound=" << TheLowerBound << "\n";
#endif
  while (true) {
#ifdef DEBUG_FLIP
    std::cerr << "Now TheUpperBound=" << TheUpperBound
              << " TheLowerBound=" << TheLowerBound << "\n";
#endif
    MyMatrix<T> Qlow = eMatIn + TheLowerBound * eMatDir;
    MyMatrix<T> Qupp = eMatIn + TheUpperBound * eMatDir;
    Tshortest<T, Tint> RecSHVlow = RetriveShortestDesc(Qlow);
    Tshortest<T, Tint> RecSHVupp = RetriveShortestDesc(Qupp);
#ifdef DEBUG_FLIP
    std::cerr << "RecSHVupp.eMin=" << RecSHVupp.eMin
              << " RecSHVlow.eMin=" << RecSHVlow.eMin << "\n";
    MyMatrix<T> ConeClassicalLow =
        GetNakedPerfectConeClassical<T>(RecSHVlow.SHV);
    std::cerr << "RankMat(ConeClassicalLow)=" << RankMat(ConeClassicalLow)
              << "\n";
    MyMatrix<T> ConeClassicalUpp =
        GetNakedPerfectConeClassical<T>(RecSHVupp.SHV);
    std::cerr << "RankMat(ConeClassicalUpp)=" << RankMat(ConeClassicalUpp)
              << "\n";
#endif
    bool test1 = RecSHVupp.eMin == RecSHVperf.eMin;
    bool test2 = TestInclusionSHV(RecSHVperf.SHV, RecSHVlow.SHV);
    if (test1) {
#ifdef DEBUG_FLIP
      std::cerr << "Return Qupp\n";
#endif
      return {std::move(Qupp), std::move(RecSHVupp)};
    }
    if (!test2) {
#ifdef DEBUG_FLIP
      std::cerr << "Qperf=\n";
      WriteMatrix(std::cerr, eMatIn);
      std::cerr << "Qlow=\n";
      WriteMatrix(std::cerr, Qlow);
      //
      std::cerr << "RecSHVperf.SHV=\n";
      WriteMatrix(std::cerr, RecSHVperf.SHV);
      std::cerr << "RecSHVlow.SHV=\n";
      WriteMatrix(std::cerr, RecSHVlow.SHV);
      std::cerr << "Return Qlow\n";
#endif
      return {std::move(Qlow), std::move(RecSHVlow)};
    }
    T TheGamma = (TheLowerBound + TheUpperBound) / 2;
    MyMatrix<T> Qgamma = eMatIn + TheGamma * eMatDir;
#ifdef DEBUG_FLIP
    std::cerr << "Qgamma=\n";
    WriteMatrix(std::cerr, Qgamma);
#endif
    Tshortest<T, Tint> RecSHVgamma = RetriveShortestDesc(Qgamma);
#ifdef DEBUG_FLIP
    std::cerr << "|RecSHVgamma.SHV|=" << RecSHVgamma.SHV.rows() << "\n";
    WriteMatrix(std::cerr, RecSHVgamma.SHV);
    MyMatrix<T> ConeClassicalGamma =
        GetNakedPerfectConeClassical<T>(RecSHVgamma.SHV);
    std::cerr << "RankMat(ConeClassicalGamma)=" << RankMat(ConeClassicalGamma)
              << "\n";
#endif
    if (RecSHVgamma.eMin >= RecSHVperf.eMin) {
      TheLowerBound = TheGamma;
    } else {
#ifdef DEBUG_FLIP
      std::cerr << "Assigning TheUpperBound to TheGamma=" << TheGamma << "\n";
#endif
      TheUpperBound = TheGamma;
      int nbRowGamma = RecSHVgamma.SHV.rows();
      for (int iRowGamma = 0; iRowGamma < nbRowGamma; iRowGamma++) {
        MyVector<Tint> eVectShort = GetMatrixRow(RecSHVgamma.SHV, iRowGamma);
#ifdef DEBUG_FLIP
        std::cerr << "iRowGamma=" << iRowGamma << " / " << nbRowGamma
                  << " eVectShort=";
        WriteVectorNoDim(std::cerr, eVectShort);
#endif
        T rVal = EvaluationQuadForm<T, Tint>(eMatDir, eVectShort);
        T qVal = EvaluationQuadForm<T, Tint>(eMatIn, eVectShort);
        if (rVal < 0) {
          T TheVal = (RecSHVperf.eMin - qVal) / rVal;
          if (TheVal < TheUpperBound) {
#ifdef DEBUG_FLIP
            std::cerr << "iRowGamma=" << iRowGamma
                      << " Assigning TheUpperBound to TheVal=" << TheVal
                      << "\n";
            std::cerr << "rVal=" << rVal << " qVal=" << qVal
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
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
Flipping_Perfect(MyMatrix<T> const &eMatIn, MyMatrix<T> const &eMatDir) {
  std::function<bool(MyMatrix<T> const &)> IsAdmissible =
      [](MyMatrix<T> const &eMat) -> bool {
    return IsPositiveDefinite<T>(eMat);
  };
  std::function<Tshortest<T, Tint>(MyMatrix<T> const &)> ShortestFunction =
      [](MyMatrix<T> const &eMat) -> Tshortest<T, Tint> {
    return T_ShortestVector<T, Tint>(eMat);
  };
  RecShort<T, Tint> eRecShort{IsAdmissible, ShortestFunction};
  return Kernel_Flipping_Perfect(eRecShort, eMatIn, eMatDir);
}

template <typename T>
MyMatrix<T> GetOnePerfectForm(LinSpaceMatrix<T> const &LinSpa) {
  int nbMat = LinSpa.ListMat.size();
  MyMatrix<T> ThePerfMat = LinSpa.SuperMat;
  while (true) {
    Tshortest<T, int> RecSHV = T_ShortestVector<T, int>(ThePerfMat);
    int nbShort = RecSHV.SHV.rows();
    MyMatrix<T> ScalMat(nbShort, nbMat);
    for (int iShort = 0; iShort < nbShort; iShort++) {
      MyVector<int> eVectShort = RecSHV.SHV.row(iShort);
      for (int iMat = 0; iMat < nbMat; iMat++) {
        T eNorm = EvaluationQuadForm<T, int>(LinSpa.ListMat[iMat], eVectShort);
        ScalMat(iShort, iMat) = eNorm;
      }
    }
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(ScalMat);
    int TheRank = eSelect.TheRank;
    if (TheRank == nbMat)
      break;
    MyVector<T> eVect = eSelect.NSP.row(0);
    MyMatrix<T> DirMat = LINSPA_GetMatrixInTspace(LinSpa, eVect);
    MyMatrix<T> eMatRet = Flipping_Perfect<T, int>(ThePerfMat, DirMat).first;
    ThePerfMat = eMatRet;
  }
  return ThePerfMat;
}

template <typename T, typename Tint, typename Tgroup>
Tgroup PERF_Automorphism(LinSpaceMatrix<T> const &LinSpa,
                         MyMatrix<T> const &ePerf, MyMatrix<Tint> const &SHV) {
  using Tidx_value = int16_t;
  MyMatrix<T> T_SHV = UniversalMatrixConversion<T, Tint>(SHV);
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV, ePerf,
                                                     LinSpa.ListComm);
  return GetStabilizerAsymmetricMatrix<std::vector<T>, Tgroup, Tidx_value>(
      WMat);
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
SimplePerfectInv<T> ComputeInvariantSimplePerfect(MyMatrix<T> const &eGram) {
  using Tidx_value = int16_t;
  int n = eGram.rows();
  Tshortest<T, Tint> RecSHV = T_ShortestVector<T, Tint>(eGram);
  MyMatrix<T> eG = eGram / RecSHV.eMin;
  int nbSHV = RecSHV.SHV.size();
  MyVector<Tint> V1(n), V2(n);
  auto f = [&](size_t i, size_t j) -> T {
    for (int iCol = 0; iCol < n; iCol++) {
      V1(i) = RecSHV.SHV(i, iCol);
      V2(i) = RecSHV.SHV(j, iCol);
    }
    return ScalarProductQuadForm<T, Tint>(eG, V1, V2);
  };
  WeightMatrix<true, T, Tidx_value> WMat(nbSHV, f);
  size_t hash = GetInvariantWeightMatrix(WMat);
  return {hash};
}

template <typename T> struct DataLinSpa {
  LinSpaceMatrix<T> LinSpa;
  bool SavingPerfect;
  bool FullDataInMemory;
  std::string PrefixPerfect;
  std::string PrefixPolyhedral;
  bool ReturnAll;
  mpz_class UpperLimitMethod4;
  bool NeedCheckStabilization;
};

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
SimplePerfect_TestEquivalence(DataLinSpa<T> const &eData,
                              MyMatrix<T> const &Gram1,
                              MyMatrix<T> const &Gram2) {
  using Telt = typename Tgroup::Telt;
  using Tidx_value = int16_t;
  Tshortest<T, Tint> RecSHV1 = T_ShortestVector<T, Tint>(Gram1);
  Tshortest<T, Tint> RecSHV2 = T_ShortestVector<T, Tint>(Gram2);
  MyMatrix<T> T_SHV1 = UniversalMatrixConversion<T, Tint>(RecSHV1.SHV);
  MyMatrix<T> T_SHV2 = UniversalMatrixConversion<T, Tint>(RecSHV2.SHV);
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat1 =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV1, Gram1,
                                                     eData.LinSpa.ListComm);
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat2 =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV2, Gram2,
                                                     eData.LinSpa.ListComm);
  std::optional<Telt> eResEquiv =
      GetEquivalenceAsymmetricMatrix<std::vector<T>, Telt>(WMat1, WMat2);
  if (!eResEquiv) {
    return {};
  }
  MyMatrix<T> M3 = RepresentVertexPermutation(T_SHV1, T_SHV2, *eResEquiv);
  if (IsIntegralMatrix(M3)) {
    MyMatrix<Tint> eMatEquiv = UniversalMatrixConversion<Tint, T>(M3);
    return eMatEquiv;
  }
  std::vector<MyVector<T>> ListMatVect;
  for (auto &eMat : eData.LinSpa.ListMat) {
    MyVector<T> eVect = SymmetricMatrixToVector(eMat);
    ListMatVect.push_back(eVect);
  }
  MyMatrix<T> ListMatVectB = MatrixFromVectorFamily(ListMatVect);
  auto f_correct = [&](MyMatrix<T> const &M) -> bool {
    if (!IsIntegralMatrix(M))
      return false;
    if (eData.NeedCheckStabilization) {
      for (auto &eMat : eData.LinSpa.ListMat) {
        MyMatrix<T> eProd = M * eMat * TransposedMat(M);
        MyVector<T> eVect = SymmetricMatrixToVector(eProd);
        std::optional<MyVector<T>> opt = SolutionMat(ListMatVectB, eVect);
        if (!opt)
          return false;
      }
    }
    return true;
  };
  auto ConvertEquiv = [](std::optional<MyMatrix<T>> const &eEq)
      -> std::optional<MyMatrix<Tint>> {
    if (!eEq)
      return {};
    MyMatrix<Tint> eMat_I = UniversalMatrixConversion<Tint, T>(*eEq);
    return eMat_I;
  };
  Tgroup GRP1 = GetStabilizerAsymmetricMatrix<std::vector<T>, Tgroup>(WMat1);
  if (GRP1.size() < eData.UpperLimitMethod4) {
    return ConvertEquiv(LinPolytopeIntegral_Isomorphism_Method4(
        T_SHV1, T_SHV2, GRP1, *eResEquiv, f_correct));
  } else {
    std::optional<MyMatrix<T>> fResEquiv =
        LinPolytopeIntegral_Isomorphism_Method8(T_SHV1, T_SHV2, GRP1,
                                                *eResEquiv);
    if (!fResEquiv)
      return {};
    if (f_correct(*fResEquiv))
      return ConvertEquiv(fResEquiv);
    return ConvertEquiv(LinPolytopeIntegral_Isomorphism_Method4(
        T_SHV1, T_SHV2, GRP1, *eResEquiv, f_correct));
  }
}

template <typename T, typename Tint, typename Tgroup>
Tgroup SimplePerfect_Stabilizer(DataLinSpa<T> const &eData,
                                MyMatrix<T> const &Gram,
                                Tshortest<T, Tint> const &RecSHV) {
  using Telt = typename Tgroup::Telt;
  using Tidx_value = int16_t;
  //
  // Functionality for checking quality of equivalences
  //
  std::vector<MyVector<T>> ListMatVect;
  for (auto &eMat : eData.LinSpa.ListMat) {
    MyVector<T> eVect = SymmetricMatrixToVector(eMat);
    ListMatVect.push_back(eVect);
  }
  MyMatrix<T> ListMatVectB = MatrixFromVectorFamily(ListMatVect);
  MyMatrix<T> T_SHV = UniversalMatrixConversion<T, Tint>(RecSHV.SHV);
  auto f_correct = [&](MyMatrix<T> const &M) -> bool {
    if (!IsIntegralMatrix(M))
      return false;
    if (eData.NeedCheckStabilization) {
      for (auto &eMat : eData.LinSpa.ListMat) {
        MyMatrix<T> eProd = M * eMat * TransposedMat(M);
        MyVector<T> eVect = SymmetricMatrixToVector(eProd);
        std::optional<MyVector<T>> opt = SolutionMat(ListMatVectB, eVect);
        if (!opt)
          return false;
      }
    }
    return true;
  };
  auto IsCorrectGroup = [&](Tgroup const &g) -> bool {
    std::vector<Telt> LGen = g.GeneratorsOfGroup();
    for (auto &eGen : LGen) {
      MyMatrix<T> M = RepresentVertexPermutation(T_SHV, T_SHV, eGen);
      if (!f_correct(M))
        return false;
    }
    return true;
  };
  //
  // Now the computation itself
  //
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV, Gram,
                                                     eData.LinSpa.ListComm);
  Tgroup GRPshv1 = GetStabilizerAsymmetricMatrix<std::vector<T>, Tgroup>(WMat);
  if (IsCorrectGroup(GRPshv1))
    return GRPshv1;
  if (GRPshv1.size() < eData.UpperLimitMethod4) {
    return LinPolytopeIntegral_Stabilizer_Method4(T_SHV, GRPshv1, f_correct);
  } else {
    Tgroup GRPshv2 = LinPolytopeIntegral_Stabilizer_Method8(T_SHV, GRPshv1);
    if (IsCorrectGroup(GRPshv2))
      return GRPshv2;
    return LinPolytopeIntegral_Stabilizer_Method4(T_SHV, GRPshv2, f_correct);
  }
}

// clang-format off
#endif  // SRC_PERFECT_TEMP_PERFECTFORM_H_
// clang-format on
