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
#ifdef DEBUG_PERFECT_FORM
  std::cerr << "m=" << n << " nbBlock=" << nbBlock << " nbSHV=" << nbSHV
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
                                              Tshortest<T, Tint> const &RecSHV) {
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
                     MyMatrix<Tint> const &SHV2, std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  using Tidx_value = int16_t;
  MyMatrix<T> T_SHV1 = UniversalMatrixConversion<T, Tint>(SHV1);
  MyMatrix<T> T_SHV2 = UniversalMatrixConversion<T, Tint>(SHV2);
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat1 =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV1, ePerf1,
                                                     LinSpa.ListComm, os);
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat2 =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV2, ePerf2,
                                                     LinSpa.ListComm, os);
  std::optional<Telt> opt =
    GetEquivalenceAsymmetricMatrix<std::vector<T>, Telt, Tidx_value>(WMat1, WMat2, os);
  if (!opt) {
    return {};
  }
  Telt const& TheEquiv = *opt;
  MyMatrix<T> M3 = RepresentVertexPermutation(T_SHV1, T_SHV2, TheEquiv);
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
    Tshortest<T, Tint> RecSHV = eRecShort.ShortestFunction(eMat);
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
    bool test = eRecShort.IsAdmissible(Qupp);
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
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
Flipping_Perfect(MyMatrix<T> const &eMatIn, MyMatrix<T> const &eMatDir,
                 std::ostream &os) {
  std::function<bool(MyMatrix<T> const &)> IsAdmissible =
      [&os](MyMatrix<T> const &eMat) -> bool {
        return IsPositiveDefinite<T>(eMat, os);
  };
  std::function<Tshortest<T, Tint>(MyMatrix<T> const &)> ShortestFunction =
      [&os](MyMatrix<T> const &eMat) -> Tshortest<T, Tint> {
    return T_ShortestVector<T, Tint>(eMat, os);
  };
  RecShort<T, Tint> eRecShort{IsAdmissible, ShortestFunction};
  return Kernel_Flipping_Perfect(eRecShort, eMatIn, eMatDir, os);
}

template <typename T, typename Tint>
MyMatrix<T> GetOnePerfectForm(LinSpaceMatrix<T> const &LinSpa,
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
    if (TheRank == nbMat)
      break;
    MyVector<T> eVect = eSelect.NSP.row(0);
    MyMatrix<T> DirMat = LINSPA_GetMatrixInTspace(LinSpa, eVect);
    MyMatrix<T> eMatRet =
        Flipping_Perfect<T, Tint>(ThePerfMat, DirMat, os).first;
    ThePerfMat = eMatRet;
  }
  return ThePerfMat;
}

template <typename T, typename Tint, typename Tgroup>
Tgroup PERF_Automorphism(LinSpaceMatrix<T> const &LinSpa,
                         MyMatrix<T> const &ePerf, MyMatrix<Tint> const &SHV,
                         std::ostream& os) {
  using Tidx_value = int16_t;
  MyMatrix<T> T_SHV = UniversalMatrixConversion<T, Tint>(SHV);
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV, ePerf,
                                                     LinSpa.ListComm, os);
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
SimplePerfectInv<T> ComputeInvariantSimplePerfect(MyMatrix<T> const &eGram,
                                                  std::ostream &os) {
  using Tidx_value = int16_t;
  int n = eGram.rows();
  Tshortest<T, Tint> RecSHV = T_ShortestVector<T, Tint>(eGram, os);
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

template <typename T, typename Tgroup, typename Fcorrect>
std::optional<MyMatrix<T>> LinPolytopeIntegral_Isomorphism_Method4(
    MyMatrix<T> const &EXT1_T, MyMatrix<T> const &EXT2_T, Tgroup const &GRP1,
    typename Tgroup::Telt const &ePerm, Fcorrect f_correct) {
  using Telt = typename Tgroup::Telt;
  for (auto &fPerm : GRP1) {
    Telt eEquivCand = fPerm * ePerm;
    MyMatrix<T> eBigMat = FindTransformation(EXT1_T, EXT2_T, eEquivCand);
    if (f_correct(eBigMat))
      return eBigMat;
  }
  return {};
}

template <typename T, typename Tgroup, typename Fcorrect>
Tgroup LinPolytopeIntegral_Stabilizer_Method4(MyMatrix<T> const &EXT_T,
                                              Tgroup const &GRPisom,
                                              Fcorrect f_correct) {
  static_assert(
      is_ring_field<T>::value,
      "Requires T to be a field in LinPolytopeIntegral_Stabilizer_Method4");
  using Telt = typename Tgroup::Telt;
  int nbVert = EXT_T.rows();
  std::vector<Telt> generatorList;
  Tgroup GRPret(nbVert);
  auto fInsert = [&](Telt const &ePerm) -> void {
    bool test = GRPret.isin(ePerm);
    if (!test) {
      generatorList.push_back(ePerm);
      GRPret = Tgroup(generatorList, nbVert);
    }
  };
  for (auto &ePerm : GRPisom) {
    MyMatrix<T> eBigMat = FindTransformation(EXT_T, EXT_T, ePerm);
    if (f_correct(eBigMat))
      fInsert(ePerm);
  }
  return GRPret;
}

template <typename T, typename Tint, typename Tgroup>
std::optional<MyMatrix<Tint>>
SimplePerfect_TestEquivalence(LinSpaceMatrix<T> const &eData,
                              MyMatrix<T> const &Gram1,
                              MyMatrix<T> const &Gram2, std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx_value = int16_t;
  Tshortest<T, Tint> RecSHV1 = T_ShortestVector<T, Tint>(Gram1, os);
  Tshortest<T, Tint> RecSHV2 = T_ShortestVector<T, Tint>(Gram2, os);
  MyMatrix<T> T_SHV1 = UniversalMatrixConversion<T, Tint>(RecSHV1.SHV);
  MyMatrix<T> T_SHV2 = UniversalMatrixConversion<T, Tint>(RecSHV2.SHV);
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat1 =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV1, Gram1,
                                                     eData.LinSpa.ListComm, os);
  WeightMatrix<false, std::vector<T>, Tidx_value> WMat2 =
      GetWeightMatrix_ListComm<false, T, Tidx_value>(T_SHV2, Gram2,
                                                     eData.LinSpa.ListComm, os);
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
Tgroup SimplePerfect_Stabilizer(LinSpaceMatrix<T> const &eData,
                                MyMatrix<T> const &Gram,
                                Tshortest<T, Tint> const &RecSHV,
                                std::ostream& os) {
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
                                                     eData.LinSpa.ListComm, os);
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
