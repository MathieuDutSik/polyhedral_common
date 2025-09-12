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
#define DEBUG_PERFECT_TSPACE_HASH
#define DEBUG_INITIAL_PERFECT
#define DEBUG_PERFECT_REPR
#define DEBUG_FLIP
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_FLIP
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

// Use a set, better asymptotics, but copy needed
template <typename Tint>
bool TestInclusionSHV_set(MyMatrix<Tint> const &TheSHVbig,
                          MyMatrix<Tint> const &TheSHVsma) {
  std::unordered_set<MyVector<Tint>> set_big;
  int nbRowBig = TheSHVbig.rows();
  for (int iRowBig = 0; iRowBig < nbRowBig; iRowBig++) {
    MyVector<Tint> V = GetMatrixRow(TheSHVbig, iRowBig);
    set_big.insert(V);
  }
  int nbRowSma = TheSHVsma.rows();
  for (int iRowSma = 0; iRowSma < nbRowSma; iRowSma++) {
    MyVector<Tint> V = GetMatrixRow(TheSHVsma, iRowSma);
    if (set_big.count(V) == 0) {
      return false;
    }
  }
  return true;
}

// Use a quadratic algorithm. But not copy needed
template <typename Tint>
bool TestInclusionSHV_quad(MyMatrix<Tint> const &TheSHVbig,
                           MyMatrix<Tint> const &TheSHVsma) {
  int nbRowSma = TheSHVsma.rows();
  int nbRowBig = TheSHVbig.rows();
  int n = TheSHVsma.cols();
  auto is_equal_vector=[&](int iRowSma, int iRowBig) -> bool {
    for (int i = 0; i < n; i++) {
      Tint const& eVal = TheSHVbig(iRowBig, i);
      Tint const& fVal = TheSHVsma(iRowSma, i);
      if (eVal != fVal) {
        return false;
      }
    }
    return true;
  };
  auto is_small_vect_in_big=[&](int iRowSma) -> bool {
    for (int iRowBig = 0; iRowBig < nbRowBig; iRowBig++) {
      if (is_equal_vector(iRowSma, iRowBig)) {
        return true;
      }
    }
    return false;
  };
  for (int iRowSma = 0; iRowSma < nbRowSma; iRowSma++) {
    if (!is_small_vect_in_big(iRowSma)) {
      return false;
    }
  }
  return true;
}

template <typename Tint>
bool TestInclusionSHV(MyMatrix<Tint> const &TheSHVbig,
                      MyMatrix<Tint> const &TheSHVsma) {
  if (TheSHVbig.rows() > 100) {
    return TestInclusionSHV_set(TheSHVbig, TheSHVsma);
  } else {
    return TestInclusionSHV_quad(TheSHVbig, TheSHVsma);
  }
}

// Memoization procedure. The computation with f_shortest is expensive
// so keeping track of what has been computed is a good idea.
// We are using the std::list for storing the data since it makes the
// const& references stable when the ListShort gets extended.
// This is a problem of std::vector that if you have a reference to an entry
// and the std::vector gets extended then the reference may get invalidated
// when we pass a threshold like 1, 2, 4, 8, 16. When we pass the threshold
// the underlying vector gets deallocated and a newer larger gets
// allocated.
template<typename Fcomp, typename Tin, typename Tout>
struct Memoization {
  Fcomp f_comp;
  std::vector<Tin> ListI;
  std::list<Tout> ListO;
  Memoization(Fcomp _f_comp) : f_comp(_f_comp) {
  }
  size_t get_index(Tin const& input) {
    size_t len = ListI.size();
    for (size_t i = 0; i < len; i++) {
      if (ListI[i] == input) {
        return i;
      }
    }
    Tout output = f_comp(input);
    ListI.push_back(input);
    ListO.push_back(output);
    return len;
  }
  Tout const& comp(Tin const& input) {
    size_t index = get_index(input);
    auto iter = ListO.cbegin();
    for (size_t u=0; u<index; u++) {
      iter++;
    }
    return *iter;
  }
};






/*
  This is the code for given a form eMatIn, a direction of change eMatDir,
  to compute a coefficient lambda>0 such that eMatIn + lambda eMatDir has
  more shortest vector that eMatIn + epsilon eMatDir for epsilon infinitely
  small.
  ---
  Two functions need to be provided on input:
  + The f_admissible function that tells you whether you are in the admissible
    domain or not.
  + The f_shortest that returns the shortest vectors.
  ---
  There are two cases in which this routine is used:
  + The perfect forms enumeration
  + The copositive simplex algorithm
  See "A simplex algorithm for rational cp-factorization" arXiv:1807.01382
  ---
  Those flipping algorithms tend to be more intricate than one would expect.
  See the following for some example description:
  + Algorithm 2 (page 8) in "A simplex algorithm for rational cp-factorization" arXiv:1807.01382
  + Algorithm 2 (page 10) in "Enumerating perfect forms" arXiv:0901.1587
 */
template <typename T, typename Tint, typename Fadmissible, typename Fshortest>
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
Kernel_Flipping_Perfect(Fadmissible f_admissible,
                        Fshortest f_shortest,
                        MyMatrix<T> const &eMatIn, MyMatrix<T> const &eMatDir,
                        [[maybe_unused]] std::ostream & os) {
#ifdef SANITY_CHECK_FLIP
  if (f_admissible(eMatDir)) {
    std::cerr << "PERF: If the direction of change belong to the space, then we are not going to find new vectors\n";
    throw TerminalException{1};
  }
#endif
  Memoization<Fadmissible,MyMatrix<T>,bool> memo_admissible(f_admissible);
  Memoization<Fshortest,MyMatrix<T>,Tshortest<T, Tint>> memo_shortest(f_shortest);
#ifdef DEBUG_FLIP
  os << "PERF: Kernel_Flipping_Perfect, case 1 (eMatIn)\n";
#endif
  Tshortest<T, Tint> const& RecSHV_in = memo_shortest.comp(eMatIn);
  //
  // Initial checks of the input
  //
#ifdef DEBUG_FLIP
  os << "PERF: Kernel_Flipping_Perfect : SHVinformation=\n";
  int nbSHV = RecSHV_in.SHV.rows();
  int n = RecSHV_in.SHV.cols();
  int nbZero_sumMat = 0;
  std::vector<int> ListScal(nbSHV);
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<Tint> V = GetMatrixRow(RecSHV_in.SHV, iSHV);
    T sumMatIn = EvaluationQuadForm(eMatIn, V);
    T sumMatDir = EvaluationQuadForm(eMatDir, V);
    int eVal = 1;
    if (sumMatDir == 0) {
      nbZero_sumMat++;
      eVal = 0;
    }
    ListScal[iSHV] = eVal;
    os << "PERF: iSHV=" << iSHV << " V=";
    for (int i = 0; i < n; i++)
      os << V(i) << " ";
    os << "sumMatIn=" << sumMatIn << " sumMatDir=" << sumMatDir << "\n";
  }
  MyMatrix<Tint> SHVface(nbZero_sumMat, n);
  int idx = 0;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    if (ListScal[iSHV] == 0) {
      for (int i = 0; i < n; i++)
        SHVface(idx, i) = RecSHV_in.SHV(iSHV, i);
      idx++;
    }
  }
  os << "PERF: SHVface=\n";
  WriteMatrix(os, SHVface);
  os << "PERF: nbZero_sumMat=" << nbZero_sumMat << "\n";
  MyMatrix<T> ConeClassicalInput =
      GetNakedPerfectConeClassical<T>(RecSHV_in.SHV);
  os << "PERF: RankMat(ConeClassicalInput)=" << RankMat(ConeClassicalInput) << "\n";
  MyMatrix<T> ConeClassicalFace = GetNakedPerfectConeClassical<T>(SHVface);
  os << "PERF: RankMat(ConeClassicalFace)=" << RankMat(ConeClassicalFace) << "\n";
#endif
  //
  // Running the algorithm
  // First loop where we find an interval where the solution may exist.
  //
  T TheUpperBound = 1;
  T TheLowerBound = 0;
#ifdef DEBUG_FLIP
  int iterLoop = 0;
#endif
  while (true) {
    MyMatrix<T> Q_upp = eMatIn + TheUpperBound * eMatDir;
#ifdef DEBUG_FLIP
    os << "PERF: CALL IsAdmissible (Q_upp)\n";
#endif
    bool test = memo_admissible.comp(Q_upp);
#ifdef DEBUG_FLIP
    iterLoop++;
    os << "PERF: iterLoop=" << iterLoop << " TheUpperBound=" << TheUpperBound
       << " TheLowerBound=" << TheLowerBound << " test=" << test << "\n";
#endif
    if (!test) {
      // We are outside, so reduce the range.
      TheUpperBound = (TheUpperBound + TheLowerBound) / 2;
    } else {
#ifdef DEBUG_FLIP
      os << "PERF: Q_upp=\n";
      WriteMatrix(os, Q_upp);
      os << "PERF: Kernel_Flipping_Perfect, case 2 (Q_upp)\n";
#endif
      Tshortest<T, Tint> const& RecSHV_upp = memo_shortest.comp(Q_upp);
#ifdef DEBUG_FLIP
      os << "ITER: RecSHV_upp.min=" << RecSHV_upp.min << "\n";
      os << "ITER: RecSHV_upp.SHV=\n";
      WriteMatrix(os, RecSHV_upp.SHV);
#endif
      if (RecSHV_upp.min == RecSHV_in.min) {
        T nLow = TheUpperBound;
        T nUpp = 2 * TheUpperBound;
        TheLowerBound = nLow;
        TheUpperBound = nUpp;
      } else {
        break;
      }
    }
  }
#ifdef SANITY_CHECK_FLIP
  {
    MyMatrix<T> Q_low = eMatIn + TheLowerBound * eMatDir;
    MyMatrix<T> Q_upp = eMatIn + TheUpperBound * eMatDir;
    if (!memo_admissible.comp(Q_low)) {
      std::cerr << "PERF: We should have Qlow being admissible\n";
      throw TerminalException{1};
    }
    if (!memo_admissible.comp(Q_upp)) {
      std::cerr << "PERF: We should have Qlow being admissible\n";
      throw TerminalException{1};
    }
    Tshortest<T, Tint> const& RecSHV_upp = memo_shortest.comp(Q_upp);
    if (RecSHV_upp.min >= RecSHV_in.min) {
      std::cerr << "PERF: We should have RecSHV_upp.min < RecSHV_in.min\n";
      throw TerminalException{1};
    }
  }
#endif
  //
  // Second while loop
  //
  // Now we must have Q_low / Q_upp admissible
  // and RecSHV_upp.min < RecSHV_in.min
#ifdef DEBUG_FLIP
  os << "PERF: FIRST LOOP FINISHED TheUpperBound=" << TheUpperBound
     << " TheLowerBound=" << TheLowerBound << "\n";
#endif
  while (true) {
#ifdef DEBUG_FLIP
    os << "PERF: Now TheLowerBound=" << TheLowerBound
       << " TheUpperBound=" << TheUpperBound << "\n";
#endif
    MyMatrix<T> Q_low = eMatIn + TheLowerBound * eMatDir;
    MyMatrix<T> Q_upp = eMatIn + TheUpperBound * eMatDir;
#ifdef DEBUG_FLIP
    os << "PERF: Kernel_Flipping_Perfect, case 3 (Q_low / Q_upp)\n";
#endif
    Tshortest<T, Tint> const& RecSHV_low = memo_shortest.comp(Q_low);
    Tshortest<T, Tint> const& RecSHV_upp = memo_shortest.comp(Q_upp);
#ifdef DEBUG_FLIP
    os << "PERF: Kernel_Flipping_Perfect, case 4 SHV_low.min=" << RecSHV_low.min << " SHV_upp.min=" << RecSHV_upp.min << "\n";
    MyMatrix<T> ConeClassicalLow =
        GetNakedPerfectConeClassical<T>(RecSHV_low.SHV);
    os << "PERF: RankMat(ConeClassicalLow)=" << RankMat(ConeClassicalLow) << "\n";
    MyMatrix<T> ConeClassicalUpp =
        GetNakedPerfectConeClassical<T>(RecSHV_upp.SHV);
    os << "PERF: RankMat(ConeClassicalUpp)=" << RankMat(ConeClassicalUpp) << "\n";
#endif
    bool test1 = RecSHV_upp.min == RecSHV_in.min;
#ifdef SANITY_CHECK_FLIP
    if (RecSHV_upp.min > RecSHV_in.min) {
      std::cerr << "PERF: The shortest vectors should be of the same norm or below the input\n";
      throw TerminalException{1};
    }
#endif
    bool test2 = TestInclusionSHV(RecSHV_in.SHV, RecSHV_low.SHV);
#ifdef DEBUG_FLIP
    os << "PERF: RecSHV_upp.min = " << RecSHV_upp.min << "\n";
    os << "PERF: RecSHV_in.min  = " << RecSHV_in.min << "\n";
    os << "PERF: test1=" << test1 << " test2=" << test2 << "\n";
#endif
    if (test1) {
#ifdef DEBUG_FLIP
      os << "PERF: Return Q_upp\n";
#endif
      return {std::move(Q_upp), std::move(RecSHV_upp)};
    }
    if (!test2 && test1) {
#ifdef DEBUG_FLIP
      os << "PERF: Qperf=\n";
      WriteMatrix(os, eMatIn);
      os << "PERF: Q_low=\n";
      WriteMatrix(os, Q_low);
      //
      os << "PERF: RecSHV_in.SHV=\n";
      WriteMatrix(os, RecSHV_in.SHV);
      os << "PERF: RecSHV_low.SHV=\n";
      WriteMatrix(os, RecSHV_low.SHV);
      os << "PERF: Return Q_low\n";
#endif
      return {std::move(Q_low), std::move(RecSHV_low)};
    }
    T TheGamma = (TheLowerBound + TheUpperBound) / 2;
    MyMatrix<T> Q_gamma = eMatIn + TheGamma * eMatDir;
#ifdef DEBUG_FLIP
    os << "PERF: Q_gamma=\n";
    WriteMatrix(os, Q_gamma);
#endif
#ifdef DEBUG_FLIP
    os << "PERF: Kernel_Flipping_Perfect, case 5 (Q_gamma)\n";
#endif
    Tshortest<T, Tint> const& RecSHV_gamma = memo_shortest.comp(Q_gamma);
#ifdef DEBUG_FLIP
    os << "|RecSHV_gamma.SHV|=" << RecSHV_gamma.SHV.rows() << "\n";
    WriteMatrix(os, RecSHV_gamma.SHV);
    MyMatrix<T> ConeClassicalGamma =
        GetNakedPerfectConeClassical<T>(RecSHV_gamma.SHV);
    os << "PERF: RankMat(ConeClassicalGamma)=" << RankMat(ConeClassicalGamma) << "\n";
#endif
    if (RecSHV_gamma.min >= RecSHV_in.min) {
      TheLowerBound = TheGamma;
    } else {
#ifdef DEBUG_FLIP
      os << "PERF: Assigning TheUpperBound to TheGamma=" << TheGamma << "\n";
#endif
      TheUpperBound = TheGamma;
      int nbRowGamma = RecSHV_gamma.SHV.rows();
      for (int iRowGamma = 0; iRowGamma < nbRowGamma; iRowGamma++) {
        MyVector<Tint> V = GetMatrixRow(RecSHV_gamma.SHV, iRowGamma);
        T rVal = EvaluationQuadForm<T, Tint>(eMatDir, V);
        T qVal = EvaluationQuadForm<T, Tint>(eMatIn, V);
#ifdef DEBUG_FLIP
        os << "PERF: iRowGamma=" << iRowGamma << " / " << nbRowGamma
           << " V=";
        for (int i=0; i<eMatIn.rows(); i++) {
          os << " " << V(i);
        }
        os << " rVal=" << rVal << " qVal=" << qVal << "\n";
#endif
        if (rVal < 0) {
          T TheVal = (RecSHV_in.min - qVal) / rVal;
          if (TheVal < TheUpperBound) {
#ifdef DEBUG_FLIP
            os << "PERF: iRowGamma=" << iRowGamma
               << " Assigning TheUpperBound to TheVal=" << TheVal << "\n";
            os << "PERF: rVal=" << rVal << " qVal=" << qVal
               << " RecSHV_in.min=" << RecSHV_in.min << "\n";
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
  using Tidx_value = int16_t;
  int n = eGram.rows();
  MyMatrix<T> eG = eGram / RecSHV.min;
  size_t n_row = RecSHV.SHV.rows();
  int nbSHV = 2 * n_row;
#ifdef DEBUG_PERFECT_TSPACE_HASH
  os << "PERF_TSPACE: ComputeInvariantPerfectTspace, nbSHV=" << nbSHV << "\n";
#endif

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
        eSum += eG(i, j) * RecSHV.SHV(iRowRed, j);
      V(i) = eSum;
    }
  };
#ifdef DEBUG_PERFECT_TSPACE_HASH
  std::map<T,size_t> map_val;
  std::map<T,size_t> map_abs_val;
#endif
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
#ifdef DEBUG_PERFECT_TSPACE_HASH
    map_val[fSum] += 1;
    T fSum_abs = T_abs(fSum);
    map_abs_val[fSum_abs] += 1;
#endif
    return fSum;
  };
  WeightMatrix<true, T, Tidx_value> WMat(nbSHV, f1, f2, os);
#ifdef DEBUG_PERFECT_TSPACE_HASH
  int iter = 0;
  for (auto & kv : map_val) {
    os << "PERF_TSPACE: iter=" << iter << " val=" << kv.first << " mult=" << kv.second << "\n";
    iter += 1;
  }
  int iter_abs = 0;
  for (auto & kv : map_abs_val) {
    os << "PERF_TSPACE: iter_abs=" << iter_abs << " val=" << kv.first << " mult=" << kv.second << "\n";
    iter_abs += 1;
  }
#endif
#ifdef DEBUG_PERFECT_TSPACE_HASH
  os << "PERF_TSPACE: ComputeInvariantPerfectTspace, We have WMat\n";
#endif
  size_t hash = GetInvariantWeightMatrix(seed, WMat, os);
#ifdef DEBUG_PERFECT_TSPACE_HASH
  os << "PERF_TSPACE: ComputeInvariantPerfectTspace, We have hash, hash=" << hash << "\n";
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
MyMatrix<T> get_scal_mat(LinSpaceMatrix<T> const &LinSpa, Tshortest<T, Tint> const& RecSHV) {
  int nbMat = LinSpa.ListMat.size();
  int nbShort = RecSHV.SHV.rows();
  MyMatrix<T> ScalMat(nbShort, nbMat);
  for (int iShort = 0; iShort < nbShort; iShort++) {
    MyVector<Tint> eVectShort = RecSHV.SHV.row(iShort);
    for (int iMat = 0; iMat < nbMat; iMat++) {
      T eNorm = EvaluationQuadForm<T, Tint>(LinSpa.ListMat[iMat], eVectShort);
      ScalMat(iShort, iMat) = eNorm;
    }
  }
  return ScalMat;
}

template <typename T, typename Tint>
bool is_perfect_in_space(LinSpaceMatrix<T> const &LinSpa, Tshortest<T, Tint> const& RecSHV) {
  int nbMat = LinSpa.ListMat.size();
  MyMatrix<T> ScalMat = get_scal_mat<T,Tint>(LinSpa, RecSHV);
  return RankMat(ScalMat) == nbMat;
}


/*
  Finds a perfect form by starting from the super matrix of the space
  and then doing some moves if that is not perfect.
  ---
  The direction of the movements are obtained by computing the
  shortest vectors, getting the kernel.
  The vector in the kernel gives direction of change.
  However, if that direction is positive definite then we will not obtain
  a new form. This scenario is explained in "Enumerating Perfect Forms" arXiv:0901.1587
  ---
  During the running of this test, the minimal norm of the matrix remains the same,
  but the number of vectors attaining it is increasing.
 */
template <typename T, typename Tint>
std::pair<MyMatrix<T>, Tshortest<T, Tint>> GetOnePerfectForm(LinSpaceMatrix<T> const &LinSpa,
                              std::ostream &os) {
  int nbMat = LinSpa.ListMat.size();
  MyMatrix<T> ThePerfMat = LinSpa.SuperMat;
  Tshortest<T, Tint> RecSHV = T_ShortestVector<T, Tint>(ThePerfMat, os);
#ifdef DEBUG_INITIAL_PERFECT
  int iter = 0;
#endif
  while (true) {
    MyMatrix<T> ScalMat = get_scal_mat<T,Tint>(LinSpa, RecSHV);
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(ScalMat);
    int TheRank = eSelect.TheRank;
#ifdef DEBUG_INITIAL_PERFECT
    os << "PERF: GetOnePerfectForm, iter=" << iter << " min=" << RecSHV.min << " |SHV|=" << RecSHV.SHV.rows() << "\n";
#endif
    if (TheRank == nbMat) {
#ifdef DEBUG_INITIAL_PERFECT
      os << "PERF: GetOnePerfectForm, returning at iter=" << iter << "\n";
#endif
      return {std::move(ThePerfMat), std::move(RecSHV)};
    }
    MyVector<T> V = eSelect.NSP.row(0);
    auto iife_get_dir=[&]() -> MyMatrix<T> {
      MyMatrix<T> M = LINSPA_GetMatrixInTspace(LinSpa, V);
      if (IsPositiveDefinite<T>(M, os)) {
        // For a positive definite matrix, we need to take the opposite
        // because we need a direction outside of the cone.
        return -M;
      }
      return M;
    };
    MyMatrix<T> DirMat = iife_get_dir();
#ifdef DEBUG_INITIAL_PERFECT
    os << "PERF: GetOnePerfectForm, iter=" << iter << " DirMat=\n";
    WriteMatrix(os, DirMat);
#endif
    auto pair = Flipping_Perfect<T, Tint>(ThePerfMat, DirMat, os);
    ThePerfMat = pair.first;
    RecSHV = pair.second;
#ifdef DEBUG_INITIAL_PERFECT
    iter += 1;
#endif
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
#ifdef DEBUG_PERFECT_REPR
  os << "PERF: TestEquivalence, det1=" << DeterminantMat(eMat1) << " det2=" << DeterminantMat(eMat2) << " opt.has_value()=" << opt.has_value() << "\n";
  os << "PERF: SHV1_T=\n";
  WriteMatrix(os, SHV1_T);
  os << "PERF: SHV2_T=\n";
  WriteMatrix(os, SHV2_T);
#endif
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
