// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_LINEARPROGRAMMING_H_
#define SRC_POLY_POLY_LINEARPROGRAMMING_H_

#include "Basic_file.h"
#include "Basic_string.h"
#include "MAT_MatrixInt.h"
#include "POLY_PolytopeFct.h"
#include "POLY_cddlib.h"
#include <string>
#include <vector>

template <typename T> struct LpSolution {
  std::string method;
  bool PrimalDefined = false;
  bool DualDefined = false;
  //
  MyVector<T> DualSolution;
  //
  T OptimalValue;
  //
  MyVector<T> DirectSolution;
  MyVector<T> DirectSolutionExt;
  Face eFace;
  int rankDirectSol = -1;
  std::string Answer;
};


template <typename T>
void PrintLpSolution(LpSolution<T> const& eSol, std::ostream & os) {
  os << "method=" << eSol.method << "\n";
  os << "PrimalDefined=" << eSol.PrimalDefined << "\n";
  os << "DualDefined=" << eSol.DualDefined << "\n";
  os << "DualSolution=" << StringVector(eSol.DualSolution) << "\n";
  os << "OptimalValue=" << eSol.OptimalValue << "\n";
  os << "DirectSolution=" << StringVector(eSol.DirectSolution) << "\n";
  os << "DirectSolutionExt=" << StringVector(eSol.DirectSolutionExt) << "\n";
  os << "rankDirectSol=" << eSol.rankDirectSol << "\n";
  os << "Answer=" << eSol.Answer << "\n";
}

template <typename T>
void WriteInputFileCdd(std::string const &FileName, MyMatrix<T> const &ListIneq,
                       MyVector<T> const &ToBeMinimized) {
  int nbRow = ListIneq.rows();
  int nbCol = ListIneq.cols();
  std::ofstream os(FileName);
  os << "H-representation\n";
  os << "begin\n";
  os << " " << nbRow << " " << nbCol << " integer\n";
  for (int iRow = 0; iRow < nbRow; iRow++) {
    for (int iCol = 0; iCol < nbCol; iCol++)
      os << " " << ListIneq(iRow, iCol);
    os << "\n";
  }
  os << "end\n";
  os << "minimize\n";
  WriteVectorNoDim(os, ToBeMinimized);
}

template <typename T>
LpSolution<T> CDD_LinearProgramming_External(MyMatrix<T> const &InequalitySet,
                                             MyVector<T> const &ToBeMinimized) {
  std::cerr << "Begin CDD_LinearProgramming_External\n";
  std::string eStr = random_string(20);
  std::string FileIne = "/tmp/LP_" + eStr + ".ine";
  std::string FileLps = "/tmp/LP_" + eStr + ".lps";
  std::string FileErr = "/tmp/LP_" + eStr + ".error";
  std::string FileDdl = "/tmp/LP_" + eStr + ".ddl";
  std::string FileLog = "/tmp/LP_" + eStr + ".log";
  std::string FileCpp = "/tmp/LP_" + eStr + ".cpp";
  auto CleanFile = [&]() -> void {
    RemoveFileIfExist(FileIne);
    RemoveFileIfExist(FileLps);
    RemoveFileIfExist(FileErr);
    RemoveFileIfExist(FileDdl);
    RemoveFileIfExist(FileLog);
    RemoveFileIfExist(FileCpp);
  };
  CleanFile();
  //
  WriteInputFileCdd(FileIne, InequalitySet, ToBeMinimized);
  //
  std::string FileTestlp2 = "testlp2_gmp";
  std::string eComm1 =
      FileTestlp2 + " " + FileIne + " 2> " + FileErr + " > " + FileLog;
  int iret1 = system(eComm1.c_str());
  if (iret1 != 0) {
    std::cerr << "iret1=" << iret1 << "\n";
    std::cerr << "Call to " << FileTestlp2 << " failed\n";
    std::cerr << "eComm1=" << eComm1 << "\n";
    throw TerminalException{1};
  }
  std::string FilelpcddcleanerCpp = "lpcddcleanerCpp";
  std::string eComm2 = FilelpcddcleanerCpp + " < " + FileLog + " > " + FileCpp;
  int iret2 = system(eComm2.c_str());
  if (iret2 != 0) {
    std::cerr << "iret2=" << iret2 << "\n";
    std::cerr << "Call to " << FilelpcddcleanerCpp << "\n";
    std::cerr << "eComm2=" << eComm2 << "\n";
    throw TerminalException{1};
  }
  //
  std::ifstream is(FileCpp);
  LpSolution<T> eSol;
  eSol.method = "cdd";
  bool PrimalDefined;
  is >> PrimalDefined;
  eSol.PrimalDefined = PrimalDefined;
  //
  bool DualDefined;
  is >> DualDefined;
  eSol.DualDefined = DualDefined;
  //
  MyVector<T> DualSolution = ReadVector<T>(is);
  eSol.DualSolution = DualSolution;
  //
  T OptimalValue;
  is >> OptimalValue;
  eSol.OptimalValue = OptimalValue;
  //
  MyVector<T> DirectSolution = ReadVector<T>(is);
  eSol.DirectSolution = DirectSolution;
  //
  MyVector<T> DirectSolutionExt = ReadVector<T>(is);
  eSol.DirectSolutionExt = DirectSolutionExt;
  //
  std::string eAnswer;
  is >> eAnswer;
  eSol.Answer = eAnswer;
  //
  CleanFile();
  return eSol;
}

template <typename T>
LpSolution<T> CDD_LinearProgramming(MyMatrix<T> const &TheEXT,
                                    MyVector<T> const &eVect) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  cdd::dd_ErrorType error = cdd::dd_NoError;
  cdd::dd_matrixdata<T> *M;
  cdd::dd_LPSolverType solver =
      cdd::dd_DualSimplex; /* either DualSimplex or CrissCross */
  cdd::dd_lpdata<T>
      *lp; /* pointer to LP data structure that is not visible by user. */
  cdd::dd_colrange j;
  int d_input;
  int nbRow, nbCol, idx;
  std::vector<int> DualSolutionPos;
  std::vector<T> DualSolutionVal;
  std::vector<int> DirectSolutionPos;
  std::vector<T> DirectSolutionVal;
  M = cdd::MyMatrix_PolyFile2Matrix(TheEXT);
  M->representation = cdd::dd_Inequality;
  d_input = TheEXT.cols();
  for (j = 1; j <= d_input; j++)
    M->rowvec[j - 1] = eVect(j - 1);
  lp = cdd::dd_Matrix2LP(M, &error);
  lp->objective = cdd::dd_LPmin;
  nbRow = TheEXT.rows();
  nbCol = TheEXT.cols();
  dd_LPSolve(lp, solver, &error);
  bool PrimalDefined = false;
  bool DualDefined = false;
  std::string eAnswer;
  switch (lp->LPS) {
  case cdd::dd_Optimal:
    // Correspondence in cddlp.c: We have dual_solution, primal_solution and
    // optimal_value
    PrimalDefined = true;
    DualDefined = true;
    eAnswer = "dd_Optimal";
    break;
  case cdd::dd_Inconsistent:
    // Corrrespondence in cddlp.c: We have dual_direction
    PrimalDefined = false;
    DualDefined = true;
    eAnswer = "dd_Inconsistent";
    break;
  case cdd::dd_DualInconsistent:
    // Correspondence in cddlp.c: Nothing actually.
    PrimalDefined = true;
    DualDefined = false;
    eAnswer = "dd_DualInconsistent";
    break;
  case cdd::dd_StrucDualInconsistent:
    // Correspondence in cddlp.c: We have primal_direction
    PrimalDefined = true;
    DualDefined = false;
    eAnswer = "dd_StructDualInconsistent";
    break;
  case cdd::dd_LPSundecided:
    // Programming of this case is missing. Please work from cdd code
    std::cerr
        << "Programming of the case dd_LPSundecided is missing. Please work\n";
    throw TerminalException{1};
  case cdd::dd_StrucInconsistent:
    // Programming of this case is missing. Please work from cdd code
    std::cerr << "Programming of the case dd_StructInconsistent is missing. "
                 "Please work\n";
    throw TerminalException{1};
  case cdd::dd_Unbounded:
    // Programming of this case is missing. Please work from cdd code
    std::cerr
        << "Programming of the case dd_Unbounded is missing. Please work\n";
    throw TerminalException{1};
  case cdd::dd_DualUnbounded:
    // Programming of this case is missing. Please work from cdd code
    std::cerr
        << "Programming of the case dd_Unbounded is missing. Please work\n";
    throw TerminalException{1};
  }
  LpSolution<T> eSol;
  eSol.method = "cdd";
  eSol.Answer = eAnswer;
  if (PrimalDefined) {
    MyVector<T> eVectDirSol(nbCol - 1);
    MyVector<T> eVectDirSolExt(nbCol);
    eVectDirSolExt(0) = 1;
    for (j = 1; j < lp->d; j++) {
      eVectDirSol(j - 1) = lp->sol[j];
      eVectDirSolExt(j) = lp->sol[j];
    }
    eSol.PrimalDefined = true;
    eSol.DirectSolution = eVectDirSol;
    eSol.DirectSolutionExt = eVectDirSolExt;
    Face eFace(nbRow);
    if (DualDefined) {
      // The comparison with values makes sense only of the dual program is
      // defined. Otherwise, what we may get is actually a primal_direction and
      // it would just not make sense with negative values for eSum.
      for (int iRow = 0; iRow < nbRow; iRow++) {
        T eSum = 0;
        for (int iCol = 0; iCol < nbCol; iCol++)
          eSum += eVectDirSolExt(iCol) * TheEXT(iRow, iCol);
        if (eSum < 0) {
          std::cerr << "CDD_LinearProgramming Error iRow=" << iRow
                    << " eSum=" << eSum << "\n";
          std::cerr << "DualDefined=" << DualDefined
                    << " PrimalDefined=" << PrimalDefined << "\n";
          std::cerr << "eVectDirSolExt =";
          for (int iCol = 0; iCol < nbCol; iCol++)
            std::cerr << " " << eVectDirSolExt(iCol);
          std::cerr << "\n";
          std::cerr << "TheEXT=\n";
          WriteMatrix(std::cerr, TheEXT);
          std::cerr << "eVect=\n";
          WriteVectorNoDim(std::cerr, eVect);
          std::cerr << "Obtained vertex solution is not valid\n";
          std::cerr << "Please debug. Before calling TerminalEnding\n";
          throw TerminalException{1};
          // TerminalEnding();
        }
        if (eSum == 0)
          eFace[iRow] = 1;
      }
      int nbIncd = eFace.count();
      if (nbIncd == 0) {
        std::cerr << "We have nbIncd=" << nbIncd
                  << " while we should have 0 < nbIncd <= with nbRow=" << nbRow
                  << "\n";
        throw TerminalException{1};
      }
    }
    eSol.eFace = eFace;
    MyMatrix<T> eMatRed = SelectRow(TheEXT, eFace);
    int rnk = RankMat(eMatRed);
    eSol.rankDirectSol = rnk;
  }
  MyVector<T> eVectDualSolution = ZeroVector<T>(nbRow);
  if (DualDefined) {
    for (j = 1; j < lp->d; j++) {
      idx = lp->nbindex[j + 1];
      if (idx > 0)
        eVectDualSolution(idx - 1) = lp->dsol[j];
    }
    eSol.DualDefined = true;
    eSol.DualSolution = eVectDualSolution;
  }
  if (PrimalDefined && DualDefined)
    eSol.OptimalValue = lp->optvalue;
  dd_FreeMatrix(M);
  dd_FreeLPData(lp);
  return eSol;
}

template <typename T>
std::optional<MyMatrix<T>> AffinizeSubspace(MyMatrix<T> const &NSP) {
  int TheDim = NSP.rows();
  int nbCol = NSP.cols();
  if (TheDim == 0)
    return {};
  int FirstVertex = -1;
  for (int iRow = 0; iRow < TheDim; iRow++) {
    if (NSP(iRow, 0) != 0)
      FirstVertex = iRow;
  }
  if (FirstVertex == -1)
    return {};
  T ePivot = NSP(FirstVertex, 0);
  MyMatrix<T> NSPred(TheDim, nbCol);
  MyVector<T> eVertex(nbCol);
  for (int i = 0; i < nbCol; i++)
    NSPred(0, i) = NSP(FirstVertex, i) / ePivot;
  int pos = 0;
  for (int iRow = 0; iRow < TheDim; iRow++)
    if (iRow != FirstVertex) {
      pos++;
      for (int i = 0; i < nbCol; i++)
        NSPred(pos, i) = NSP(iRow, i) - NSP(iRow, 0) * NSPred(0, i);
    }
  return NSPred;
}

template <typename T>
std::optional<MyMatrix<T>>
ComputeStandardAffineSubspace(MyMatrix<T> const &ListEqua) {
  MyMatrix<T> NSP = NullspaceTrMat(ListEqua);
  return AffinizeSubspace(NSP);
}

template <typename T>
bool TestRealizabilityInequalitiesEqualities(MyMatrix<T> const &ListIneq,
                                             MyMatrix<T> const &ListEqua,
                                             MyVector<T> const &ToBeMinimized) {
  std::optional<MyMatrix<T>> EquivArr = ComputeStandardAffineSubspace(ListEqua);
  if (!EquivArr)
    return false;
  MyMatrix<T> NSPred = EquivArr.TheEquiv;
  MyMatrix<T> ListIneqRed = ListIneq * (NSPred.transpose());
  MyVector<T> ToBeMinimizedRed = NSPred * ToBeMinimized;
  LpSolution<T> eSol = CDD_LinearProgramming(ListIneqRed, ToBeMinimizedRed);
  if (eSol.PrimalDefined && eSol.DualDefined)
    return true;
  return false;
}

template <typename T>
LpSolution<T> GLPK_LinearProgramming(MyMatrix<T> const &ListIneq,
                                     MyVector<T> const &ToBeMinimized) {
  int dimTot = ToBeMinimized.size();
  MyMatrix<T> ListEqua(0, dimTot);
  GLPKoption eGLPKoption;
  LpSolutionSimple<double> eResSimple =
      GLPK_LinearProgramming_Kernel_Dense_LIBRARY(ListEqua, ListIneq,
                                                  ToBeMinimized, eGLPKoption);
  if (!eResSimple.PrimalDefined) {
    return CDD_LinearProgramming(ListIneq, ToBeMinimized);
  }
  int nbIneq = ListIneq.rows();
  std::cerr << "nbIneq=" << nbIneq << "\n";
  std::cerr << "|eResSimple.RowStatus|=" << eResSimple.RowStatus.size() << "\n";
  std::cerr << "|eResSimple.ColumnStatus|=" << eResSimple.ColumnStatus.size()
            << "\n";
  std::vector<int> ListRowSelect;
  for (int iIneq = 0; iIneq < nbIneq; iIneq++)
    if (eResSimple.RowStatus(iIneq) == 3)
      ListRowSelect.push_back(iIneq);
  MyMatrix<T> ListIneqSel = SelectRow(ListIneq, ListRowSelect);
  MyMatrix<T> NSP = NullspaceTrMat(ListIneqSel);
  int dimNSP = NSP.rows();
  if (dimNSP == 0) {
    return CDD_LinearProgramming(ListIneq, ToBeMinimized);
  }
  MyVector<T> TheVert;
  int nbCol = NSP.cols();
  if (dimNSP == 1) {
    MyVector<T> eNSP = GetMatrixRow(NSP, 0);
    TheVert = eNSP / eNSP(0);
  } else {
    MyVector<T> eFirstPoint;
    int iNSPselect = -1;
    for (int iNSP = 0; iNSP < dimNSP; iNSP++) {
      if (iNSPselect == -1) {
        T eVAL = NSP(iNSP, 0);
        if (eVAL != 0) {
          eFirstPoint = GetMatrixRow(NSP, iNSP) / eVAL;
          iNSPselect = iNSP;
        }
      }
    }
    if (iNSPselect == -1)
      return CDD_LinearProgramming(ListIneq, ToBeMinimized);
    MyMatrix<T> ColumnSpace(dimNSP, nbCol);
    AssignMatrixRow(ColumnSpace, 0, eFirstPoint);
    int pos = 0;
    for (int iNSP = 0; iNSP < dimNSP; iNSP++) {
      if (iNSP != iNSPselect) {
        pos++;
        MyVector<T> eNSP = GetMatrixRow(NSP, iNSP);
        MyVector<T> eVec = eNSP - eNSP(0) * eFirstPoint;
        AssignMatrixRow(ColumnSpace, pos, eVec);
      }
    }
    MyMatrix<T> SEC_ListIneq = ListIneq * TransposedMat(ColumnSpace);
    MyVector<T> SEC_ToBeMinimized = ColumnSpace * ToBeMinimized;
    LpSolution<T> TheLP = CDD_LinearProgramming(ListIneq, ToBeMinimized);
    if (TheLP.PrimalDefined && TheLP.DualDefined)
      TheVert = TransposedMat(ColumnSpace) * TheLP.DirectSolutionExt;
    else
      return CDD_LinearProgramming(ListIneq, ToBeMinimized);
  }
  T optimal = ScalarProduct(ToBeMinimized, TheVert);
  MyVector<T> TheVertRed(nbCol - 1);
  for (int i = 0; i < nbCol - 1; i++)
    TheVertRed(i) = TheVert(i + 1);
  int nbRow = ListIneq.rows();
  Face eFace(nbRow);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> eRow = GetMatrixRow(ListIneq, iRow);
    T scal = ScalarProduct(eRow, TheVert);
    if (scal < 0)
      return CDD_LinearProgramming(ListIneq, ToBeMinimized);
    if (scal == 0)
      eFace[iRow] = 1;
  }
  //
  LpSolution<T> eRes;
  eRes.method = "glpk";
  eRes.PrimalDefined = true;
  eRes.DualDefined = true;
  //  MyVector<T> DualSolution;
  eRes.OptimalValue = optimal;
  //
  eRes.DirectSolution = TheVertRed;
  eRes.DirectSolutionExt = TheVert;
  eRes.eFace = eFace;
  eRes.Answer = "dd_Optimal";
  return eRes;
}

template <typename T>
LpSolution<T> CDD_LinearProgramming_BugSearch(MyMatrix<T> const &TheEXT,
                                              MyVector<T> const &eVect) {
  LpSolution<T> eSol1 = CDD_LinearProgramming(TheEXT, eVect);
  LpSolution<T> eSol2 = CDD_LinearProgramming_External(TheEXT, eVect);
  if (eSol1.PrimalDefined != eSol2.PrimalDefined ||
      eSol1.DualDefined != eSol2.DualDefined) {
    WriteInputFileCdd("bugSearch.ine", TheEXT, eVect);
    std::cerr << "We find the bug we were after\n";
    throw TerminalException{1};
  }
  return eSol1;
}

template <typename T> MyMatrix<T> Polytopization(MyMatrix<T> const &EXT) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  T eVal, prov;
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyMatrix<T> nMat(nbRow, nbCol + 1);
  MyMatrix<T> eBasis(nbCol, nbCol);
  MyVector<T> eVect(nbCol + 1);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    eVal = -1;
    nMat(iRow, 0) = eVal;
    for (int iCol = 0; iCol < nbCol; iCol++)
      nMat(iRow, iCol + 1) = EXT(iRow, iCol);
  }
  for (int iCol = 0; iCol <= nbCol; iCol++) {
    T eSum = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      eSum += nMat(iRow, iCol);
    eVect(iCol) = eSum;
  }
  //
  LpSolution<T> eSol = CDD_LinearProgramming(nMat, eVect);
  if (!eSol.PrimalDefined) {
    std::cerr << "The optimization resulted in a result whose primal is not "
                 "defined\n";
    std::cerr << "This likely means that the set of inequalities on input does "
                 "not define\n";
    std::cerr << "A full dimensional polytope, i.e. that the inequalitity "
                 "together imply some equalities\n";
    throw TerminalException{1};
  }
  MyVector<T> SolDir = eSol.DirectSolution;
  ZeroAssignation(eBasis);
  for (int iCol = 0; iCol < nbCol; iCol++)
    eBasis(iCol, 0) = SolDir(iCol);
  int iColSelect = -1;
  for (int iCol = 0; iCol < nbCol; iCol++)
    if (iColSelect == -1) {
      if (eBasis(iCol, 0) != 0)
        iColSelect = iCol;
    }
  if (iColSelect == -1) {
    std::cerr << "Apparently, we did not find the column\n";
    std::cerr << "Please solve error in Polytopization\n";
    throw TerminalException{1};
  }
  int iRowWrite = 0;
  eVal = 1;
  for (int iCol = 0; iCol < nbCol; iCol++)
    if (iCol != iColSelect) {
      iRowWrite++;
      eBasis(iCol, iRowWrite) = eVal;
    }
  MyMatrix<T> EXTret = EXT * eBasis;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T eMult = 1 / EXTret(iRow, 0);
    for (int iCol = 0; iCol < nbCol; iCol++)
      EXTret(iRow, iCol) = EXTret(iRow, iCol) * eMult;
  }
  return EXTret;
}

template <typename T> struct SolVertex {
  Face eInc;
  MyVector<T> eVert;
};

template <typename T>
MyMatrix<T> SetIsobarycenter(MyMatrix<T> const& EXT) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyVector<T> eVect = MyVector<T>(nbCol);
  for (int iCol = 0; iCol < nbCol; iCol++) {
    T eSum = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      eSum += EXT(iRow, iCol);
    eSum = eSum / nbRow;
    eVect(iCol) = eSum;
  }
  MyMatrix<T> nMat(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    nMat(iRow, 0) = EXT(iRow, 0);
    for (int iCol = 1; iCol < nbCol; iCol++)
      nMat(iRow, iCol) = EXT(iRow, iCol) - eVect(iCol);
  }
  return nMat;
}

template <typename T>
vectface Kernel_FindVertices(MyMatrix<T> const &EXT, size_t const &nb) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyVector<T> eVect = MyVector<T>(nbCol);
  MyVector<T> TheVert = MyVector<T>(nbCol);
  vectface ListFace(EXT.rows());
  while (true) {
    for (int iCol = 0; iCol < nbCol; iCol++) {
      int a = random();
      int b = random();
      T eVal = a - b;
      eVect(iCol) = eVal;
    }
    LpSolution<T> eSol = CDD_LinearProgramming(EXT, eVect);
    MyVector<T> SolDir = eSol.DirectSolution;
    TheVert(0) = 1;
    for (int iCol = 0; iCol < nbCol - 1; iCol++)
      TheVert(iCol + 1) = SolDir(iCol);
    Face eInc(nbRow);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      T eSum = 0;
      for (int iCol = 0; iCol < nbCol; iCol++)
        eSum += EXT(iRow, iCol) * TheVert(iCol);
      if (eSum == 0)
        eInc[iRow] = 1;
    }
    MyMatrix<T> RnkMat = SelectRow(EXT, eInc);
    int TheRank = RankMat(RnkMat);
    if (TheRank == nbCol - 1) {
      ListFace.push_back(eInc);
      if (ListFace.size() == nb)
        return ListFace;
    }
  }
}

template <typename T>
vectface FindVertices(MyMatrix<T> const &EXT, int const &nb) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  MyMatrix<T> EXT_B = ColumnReduction(EXT);
  MyMatrix<T> EXT_C = Polytopization(EXT_B);
  MyMatrix<T> EXT_D = SetIsobarycenter(EXT_C);
  return Kernel_FindVertices(EXT_D, nb);
}

template <typename T> Face FindOneInitialVertex(MyMatrix<T> const &TheEXT) {
  return FindVertices(TheEXT, 1)[0];
}

template <typename T>
MyMatrix<T> LP_GetExpressionForLP(MyMatrix<T> const &EXT) {
  MyVector<T> OneVert = EXT.row(0);
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyMatrix<T> ListDiff(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> eDiff = GetMatrixRow(EXT, iRow) - OneVert;
    AssignMatrixRow(ListDiff, iRow, eDiff);
  }
  MyMatrix<T> TheBasis = RowReduction(ListDiff);
  int TheDim = TheBasis.rows();
  MyMatrix<T> NewListCoord(nbRow, TheDim + 1);
  MyVector<T> eEntOne(1);
  eEntOne(0) = 1;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> eVect1 = ListDiff.row(iRow);
    std::optional<MyVector<T>> opt = SolutionMat(TheBasis, eVect1);
    const MyVector<T> &eVect2 = *opt;
    MyVector<T> eVect3 = Concatenation(eEntOne, eVect2);
    NewListCoord.row(iRow) = eVect3;
  }
  MyVector<T> eIso = Isobarycenter(NewListCoord);
  eIso(0) = 0;
  MyMatrix<T> RetMat(nbRow, TheDim + 1);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> V = GetMatrixRow(NewListCoord, iRow) - eIso;
    AssignMatrixRow(RetMat, iRow, V);
  }
  return RetMat;
}

template <typename T>
Face FindViolatedFace(MyMatrix<T> const &EXT, MyVector<T> const &eVect) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  std::optional<MyVector<T>> opt = SolutionMat(EXT, eVect);
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  int nbColVect = eVect.size();
  if (nbCol != nbColVect) {
    std::cerr << "We have nbCol=" << nbCol << " but nbColVect=" << nbColVect
              << "\n";
    throw TerminalException{1};
  }
  if (!opt) {
    std::cerr << "Error in the input of FindVolatedFace\n";
    std::cerr << "The vector eVect should belong to the space spanned by EXT\n";
    throw TerminalException{1};
  }
  MyMatrix<T> EXT2(nbRow + 1, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++)
    EXT2.row(iRow) = EXT.row(iRow);
  MyVector<T> uVect = -eVect;
  AssignMatrixRow(EXT2, nbRow, uVect);
  MyMatrix<T> EXTpoly = Polytopization<T>(EXT2);
  MyMatrix<T> EXT_lp = LP_GetExpressionForLP(EXTpoly);

  int TheDim = EXT_lp.cols();
  MyVector<T> TheCritFacet = EXT_lp.row(nbRow);
  auto GetOneFacetDefiningVect = [&]() -> MyVector<T> {
    MyVector<T> eMinimize(TheDim);
    while (true) {
      // quite enough for us
      int siz = 3;
      int sizTot = 1 + 2 * siz;
      for (int iCol = 0; iCol < TheDim; iCol++) {
        int a = random() % sizTot;
        T eVal = a - siz;
        eMinimize(iCol) = eVal;
      }
      LpSolution<T> eSol = CDD_LinearProgramming(EXT_lp, eMinimize);
      if (eSol.PrimalDefined && eSol.DualDefined)
        if (eSol.rankDirectSol == TheDim - 1)
          return eMinimize;
    }
  };
  int nbIter = 20;
  while (true) {
    MyVector<T> eMinimize = GetOneFacetDefiningVect();
    for (int iter = 0; iter < nbIter; iter++) {
      LpSolution<T> eSol = CDD_LinearProgramming(EXT_lp, eMinimize);
      if (eSol.PrimalDefined && eSol.DualDefined &&
          eSol.rankDirectSol == TheDim - 1) {
        Face eFace = eSol.eFace;
        if (eFace[nbRow] == 0) {
          SelectionRowCol<T> eSelect = TMat_SelectRowCol(EXT);
          std::vector<int> ListColSelect = eSelect.ListColSelect;
          MyMatrix<T> EXTred = SelectColumn(EXT, ListColSelect);
          MyVector<T> eVectRed = SelectColumnVector(eVect, ListColSelect);
          MyVector<T> eFacet = FindFacetInequality(EXTred, eFace);
          T eScal = eFacet.dot(eVectRed);
          if (eScal >= 0) {
            std::cerr << "Error, we should have negative scalar product\n";
            std::cerr << "In order for the facet to separate\n";
            throw TerminalException{1};
          }
          Face eFaceRet(nbRow);
          for (int iRow = 0; iRow < nbRow; iRow++)
            eFaceRet[iRow] = eFace[iRow];
          return eFaceRet;
        }
      }
      eMinimize -= TheCritFacet;
    }
    nbIter++;
  }
}

template <typename T>
Face FindViolatedFaceFast(MyMatrix<T> const &EXT, MyVector<T> const &eVect) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  MyMatrix<T> EXT_ext(nbRow, nbCol+1);
  MyVector<T> ToMinimize(nbCol+1);
  for (int iRow=0; iRow<nbRow; iRow++) {
    EXT_ext(iRow,0);
    for (int iCol=0; iCol<nbCol; iCol++)
      EXT_ext(iRow,iCol+1) = EXT(iRow,iCol);
  }
  ToMinimize(0) = 0;
  for (int iCol=0; iCol<nbCol; iCol++)
    ToMinimize(iCol+1) = eVect(iCol);
  LpSolution<T> eSol = CDD_LinearProgramming(EXT_ext, ToMinimize);
  Face eFace(nbRow);
  for (int i_row=0; i_row<nbRow; i_row++) {
    T sum = 0;
    for (int i_col=0; i_col<nbCol; i_col++) {
      sum += EXT(i_row,i_col) * eSol.DirectSolution(i_col);
    }
    if (sum == 0)
      eFace[i_row] = 1;
  }
#ifdef DEBUG
  MyVector<T> eFAC = FindFacetInequalityCheck(EXT, eFace);
  T scal = eVect.dot(eFAC);
  if (scal >= 0) {
    std::cerr << "Faied to find a correct solution\n";
    throw TerminalException{1};
  }
#endif
  return eFace;
}


template <typename T> struct PosRelRes {
  bool eTestExist;
  MyVector<T> InternalVector;
  MyVector<T> TheRelat;
};

template <typename T>
PosRelRes<T>
SearchPositiveRelationSimple_DualMethod(MyMatrix<T> const &ListVect) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  int nbRow = ListVect.rows();
  int nbCol = ListVect.cols();
  MyMatrix<T> ListVectExt(nbRow, nbCol + 1);
  MyVector<T> eMinimize(nbCol + 1);
  eMinimize(0) = 0;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    ListVectExt(iRow, 0) = -1;
    for (int iCol = 0; iCol < nbCol; iCol++) {
      ListVectExt(iRow, iCol + 1) = ListVect(iRow, iCol);
      eMinimize(iCol + 1) += ListVect(iRow, iCol);
    }
  }
  LpSolution<T> eSol = CDD_LinearProgramming(ListVectExt, eMinimize);
  //  LpSolution<T> eSol=CDD_LinearProgramming_External(ListVectExt, eMinimize);
  PosRelRes<T> eResult;
  bool IsDone = false;
  if (eSol.PrimalDefined && eSol.DualDefined) {
    IsDone = true;
    eResult.eTestExist = false;
    eResult.InternalVector = eSol.DirectSolution;
    for (int iRow = 0; iRow < nbRow; iRow++) {
      T eScal = 0;
      for (int iCol = 0; iCol < nbCol; iCol++)
        eScal += eSol.DirectSolution(iCol) * ListVect(iRow, iCol);
      if (eScal <= 0) {
        std::cerr << "Error in SearchPositiveRelationSimple_DualMethod 1\n";
        std::cerr << "We have eScal=" << eScal << "\n";
        throw TerminalException{1};
      }
    }
  }
  if (!eSol.PrimalDefined && eSol.DualDefined) {
    IsDone = true;
    eResult.eTestExist = true;
    eResult.TheRelat = eSol.DualSolution;
    for (int iRow = 0; iRow < nbRow; iRow++)
      if (eSol.DualSolution(iRow) < 0) {
        std::cerr << "Error in SearchPositiveRelationSimple_DualMethod 2\n";
        std::cerr << "iRow=" << iRow << "\n";
        std::cerr << "We have DualSol=" << eSol.DualSolution(iRow) << " < 0\n";
        throw TerminalException{1};
      }
    for (int iCol = 0; iCol < nbCol; iCol++) {
      T eSum = 0;
      for (int iRow = 0; iRow < nbRow; iRow++)
        eSum += eSol.DualSolution(iRow) * ListVect(iRow, iCol);
      if (eSum != 0) {
        std::cerr << "Error in SearchPositiveRelationSimple_DualMethod 2\n";
        std::cerr << "iCol=" << iCol << "\n";
        std::cerr << "We have eSum=" << eSum << "\n";
        throw TerminalException{1};
      }
    }
  }
  if (!IsDone) {
    std::cerr << "Error. No value assigned\n";
    throw TerminalException{1};
  }
  return eResult;
}

struct Constraint {
  std::vector<int> ListStrictlyPositive;
  std::vector<int> ListPositive;
  std::vector<std::vector<int>> ListSetStrictPositive;
};

template <typename T>
PosRelRes<T> SearchPositiveRelation(MyMatrix<T> const &ListVect,
                                    Constraint const &eConstraint) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  MyMatrix<T> NSP = NullspaceMat(ListVect);
  int nbVect = ListVect.rows();
  int nbRelation = NSP.rows();
  std::vector<MyVector<T>> ListInequalities;
  for (auto &eVect : eConstraint.ListStrictlyPositive) {
    MyVector<T> eV(nbRelation + 1);
    eV(0) = -1;
    for (int iRel = 0; iRel < nbRelation; iRel++)
      eV(iRel + 1) = NSP(iRel, eVect);
    ListInequalities.push_back(eV);
  }
  for (auto &eVect : eConstraint.ListPositive) {
    MyVector<T> eV(nbRelation + 1);
    eV(0) = 0;
    for (int iRel = 0; iRel < nbRelation; iRel++)
      eV(iRel + 1) = NSP(iRel, eVect);
    ListInequalities.push_back(eV);
  }
  for (auto &eSet : eConstraint.ListSetStrictPositive) {
    MyVector<T> eV = ZeroVector<T>(nbRelation + 1);
    eV(0) = -1;
    for (auto &eVect : eSet) {
      for (int iRel = 0; iRel < nbRelation; iRel++)
        eV(iRel + 1) += NSP(iRel, eVect);
    }
    ListInequalities.push_back(eV);
  }
  MyVector<T> ToBeMinimized(nbRelation + 1);
  ToBeMinimized(0) = 0;
  for (int iRel = 0; iRel < nbRelation; iRel++) {
    T eSum = 0;
    for (int iVect = 0; iVect < nbVect; iVect++)
      eSum += NSP(iRel, iVect);
    ToBeMinimized(iRel + 1) = eSum;
  }
  MyMatrix<T> MatInequalities = MatrixFromVectorFamily(ListInequalities);
  LpSolution<T> eSol = CDD_LinearProgramming(MatInequalities, ToBeMinimized);
  PosRelRes<T> eResult;
  if (eSol.PrimalDefined && eSol.DualDefined) {
    MyVector<T> DirSol = eSol.DirectSolution;
    eResult.eTestExist = true;
    MyVector<T> eVectRel(nbRelation);
    for (int iRel = 0; iRel < nbRelation; iRel++)
      eVectRel(iRel) = DirSol(iRel);
    MyVector<T> TheRelat = ProductVectorMatrix(eVectRel, NSP);
    for (int iVect = 0; iVect < nbVect; iVect++) {
      if (TheRelat(iVect) < 0) {
        std::cerr << "We have a clear and present error\n";
        throw TerminalException{1};
      }
    }
    eResult.TheRelat = TheRelat;
  } else {
    eResult.eTestExist = false;
  }
  return eResult;
}

template <typename T>
PosRelRes<T> SearchPositiveRelationSimple(MyMatrix<T> const &ListVect) {
  int nbVect = ListVect.rows();
  std::vector<int> ListStrictlyPositive;
  std::vector<int> ListPositive(nbVect);
  for (int iVect = 0; iVect < nbVect; iVect++)
    ListPositive[iVect] = iVect;
  std::vector<std::vector<int>> ListSetStrictPositive = {ListPositive};
  Constraint eConstraint{ListStrictlyPositive, ListPositive,
                         ListSetStrictPositive};
  return SearchPositiveRelation(ListVect, eConstraint);
}

template <typename T>
std::optional<MyVector<T>>
SolutionMatNonnegative_Version1(MyMatrix<T> const &ListVect,
                                MyVector<T> const &eVect) {
  int nbVect = ListVect.rows();
  int nbCol = ListVect.cols();
  MyMatrix<T> InputListVect(nbVect + 1, nbCol);
  for (int iVect = 0; iVect < nbVect; iVect++)
    InputListVect.row(iVect) = ListVect.row(iVect);
  for (int iCol = 0; iCol < nbCol; iCol++)
    InputListVect(nbVect, iCol) = -eVect(iCol);
  //
  std::vector<int> ListStrictlyPositive;
  std::vector<int> ListPositive(nbVect + 1);
  for (int iVect = 0; iVect <= nbVect; iVect++)
    ListPositive[iVect] = iVect;
  std::vector<int> Part1(nbVect);
  for (int iVect = 0; iVect < nbVect; iVect++)
    Part1[iVect] = iVect;
  std::vector<int> Part2 = {nbVect};
  std::vector<std::vector<int>> ListSetStrictPositive = {Part1, Part2};
  Constraint eConstraint{ListStrictlyPositive, ListPositive,
                         ListSetStrictPositive};
  //
  PosRelRes<T> PRR = SearchPositiveRelation(InputListVect, eConstraint);
  if (!PRR.eTestExist)
    return {};
  MyVector<T> TheSol(nbVect);
  for (int iVect = 0; iVect < nbVect; iVect++)
    TheSol(iVect) = PRR.TheRelat(iVect) / PRR.TheRelat(nbVect);
  return TheSol;
}

template <typename T>
std::optional<MyVector<T>>
SolutionMatNonnegative_LP(MyMatrix<T> const &ListVect,
                          MyVector<T> const &eVect) {
  int nbVect = ListVect.rows();
  int nbCol = ListVect.cols();
  MyMatrix<T> ListIneq(nbVect, nbCol + 1);
  MyVector<T> eIneq(nbCol + 1);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    ListIneq(iVect, 0) = 0;
    for (int iCol = 0; iCol < nbCol; iCol++)
      ListIneq(iVect, 1 + iCol) = ListVect(iVect, iCol);
  }
  eIneq(0) = 0;
  for (int iCol = 0; iCol < nbCol; iCol++)
    eIneq(1 + iCol) = eVect(iCol);
  //
  LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, eIneq);
  if (!eSol.DualDefined) {
    return {};
  }
  MyVector<T> TheRet = -eSol.DualSolution;
  return TheRet;
}

template <typename T> struct SolutionMatNonnegativeComplete {
  std::optional<MyVector<T>> ExtremeRay;
  std::optional<MyVector<T>> SolNonnegative;
};

template <typename T>
SolutionMatNonnegativeComplete<T>
GetSolutionMatNonnegativeComplete(MyMatrix<T> const &ListVect,
                                  MyVector<T> const &eVect) {
  int nbVect = ListVect.rows();
  int nbCol = ListVect.cols();
  MyMatrix<T> ListIneq(nbVect, nbCol + 1);
  MyVector<T> eIneq(nbCol + 1);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    ListIneq(iVect, 0) = 0;
    for (int iCol = 0; iCol < nbCol; iCol++)
      ListIneq(iVect, 1 + iCol) = ListVect(iVect, iCol);
  }
  eIneq(0) = 0;
  for (int iCol = 0; iCol < nbCol; iCol++)
    eIneq(1 + iCol) = eVect(iCol);
  //
  LpSolution<T> eSol = CDD_LinearProgramming(ListIneq, eIneq);
  auto GetSolNonnegative = [&]() -> std::optional<MyVector<T>> {
    if (!eSol.DualDefined) {
      return {};
    }
    MyVector<T> TheRet = -eSol.DualSolution;
    return TheRet;
  };
  std::optional<MyVector<T>> SolNonnegative = GetSolNonnegative();
  auto GetExtremeRay = [&]() -> std::optional<MyVector<T>> {
    if (eSol.PrimalDefined) {
      MyVector<T> V = eSol.DirectSolution;
      MyVector<T> LScal = ListVect * V;
      for (int iVect = 0; iVect < nbVect; iVect++) {
        if (LScal(iVect) < 0) {
          std::cerr << "LScal=" << StringVectorGAP(LScal) << "\n";
          std::cerr << "Negative value at iVect=" << iVect
                    << " LScal(iVect)=" << LScal(iVect) << "\n";
          throw TerminalException{1};
        }
      }
      T eScal = V.dot(eVect);
      if (eScal > 0) {
        std::cerr << "eScal=" << eScal << "\n";
        std::cerr << "The direction is not a counter example\n";
        throw TerminalException{1};
      }
      if (eScal == 0 && !SolNonnegative) {
        std::cerr << "eScal = 0 and no Non-negative solution were found\n";
        std::cerr << "Possibly a bug here\n";
        throw TerminalException{1};
      }
      return V;
    }
    return {};
  };
  std::optional<MyVector<T>> ExtremeRay = GetExtremeRay();
  return {ExtremeRay, SolNonnegative};
}

template <typename T>
std::optional<MyVector<T>>
SolutionMatNonnegative_Check(MyMatrix<T> const &ListVect,
                             MyVector<T> const &eVect) {
  std::optional<MyVector<T>> opt1 =
      SolutionMatNonnegative_Version1(ListVect, eVect);
  std::optional<MyVector<T>> opt2 = SolutionMatNonnegative_LP(ListVect, eVect);
  if (opt1 && !opt2) {
    std::cerr << "opt1 is defined but not opt2, incoherent\n";
    throw TerminalException{1};
  }
  if (!opt1 && opt2) {
    std::cerr << "opt1 is not defined but opt2 is, incoherent\n";
    throw TerminalException{1};
  }
  if (!opt1 && !opt2)
    return {};
  MyVector<T> const &V = *opt1;
  int nbRow = V.size();
  if (nbRow != ListVect.rows()) {
    std::cerr << "Dimension incoherency\n";
    throw TerminalException{1};
  }
  for (int iRow = 0; iRow < nbRow; iRow++) {
    if (V(iRow) < 0) {
      std::cerr << "V should be non-negative\n";
      throw TerminalException{1};
    }
  }
  MyVector<T> prod = ListVect.transpose() * V;
  if (prod != eVect) {
    std::cerr << "prod does not matcheVect\n";
    throw TerminalException{1};
  }
  return V;
}

template <typename T>
std::optional<MyVector<T>>
SolutionMatNonnegative_FailSafe(MyMatrix<T> const &ListVect,
                                MyVector<T> const &eVect) {
  auto local_terminate = [&]() -> void {
    std::string FILE_FAC = "DEBUG_Nonnegative_FAC";
    WriteMatrixFile(FILE_FAC, ListVect);
    std::string FILE_INEQ = "DEBUG_Nonnegative_INEQ";
    WriteVectorFile(FILE_INEQ, eVect);
    throw TerminalException{1};
  };
  std::optional<MyVector<T>> opt_LP =
      SolutionMatNonnegative_LP(ListVect, eVect);
  if (opt_LP) {
    MyVector<T> const &V = *opt_LP;
    int nbRow = V.size();
    if (nbRow != ListVect.rows()) {
      std::cerr << "Dimension incoherency\n";
      local_terminate();
    }
    bool HasError = false;
    for (int iRow = 0; iRow < nbRow; iRow++) {
      T val = V(iRow);
      if (val < 0) {
        double val_d = UniversalScalarConversion<double, T>(val);
        std::cerr << "iRow=" << iRow << " val=" << val << " val_d=" << val_d
                  << "\n";
        HasError = true;
      }
    }
    MyVector<T> prod = ListVect.transpose() * V;
    if (prod != eVect) {
      std::cerr << "prod does not matcheVect\n";
      std::cerr << "|prod|=" << prod.rows() << " / " << prod.cols() << "\n";
      std::cerr << "|eVect|=" << eVect.rows() << " / " << eVect.cols() << "\n";
      for (int iRow = 0; iRow < prod.rows(); iRow++) {
        double prod_d = UniversalScalarConversion<double, T>(prod(iRow, 0));
        double eVect_d = UniversalScalarConversion<double, T>(eVect(iRow, 0));
        std::cerr << "iRow=" << iRow << " prod_d=" << prod_d
                  << " eVect_d=" << eVect_d << "\n";
      }
      HasError = true;
    }
    if (HasError) {
      local_terminate();
    }
    return V;
  }
  std::optional<MyVector<T>> opt_V1 =
      SolutionMatNonnegative_Version1(ListVect, eVect);
  if (opt_V1) {
    std::cerr << "opt_V1 is defined but not opt_LP, incoherent\n";
    local_terminate();
  }
  return {};
}

template <typename T>
std::optional<MyVector<T>> SolutionMatNonnegative(MyMatrix<T> const &ListVect,
                                                  MyVector<T> const &eVect) {
  //  return SolutionMatNonnegative_Check(ListVect, eVect);
  return SolutionMatNonnegative_FailSafe(ListVect, eVect);
}

template <typename T> Face ComputeSkeletonClarkson(MyMatrix<T> const &FACinp) {
  MyMatrix<T> FAC = ColumnReduction(FACinp);
  int n_fac = FAC.rows();
  int dim = FAC.cols();
  Face f_adj(n_fac * n_fac);
  for (int i_fac = 0; i_fac < n_fac; i_fac++) {
    MyMatrix<T> Equa(dim, 1);
    for (int i = 0; i < dim; i++)
      Equa(i, 0) = FAC(i_fac, i);
    MyMatrix<T> NSP = NullspaceMat(Equa);
    std::vector<int> ListIdx;
    for (int j_fac = 0; j_fac < i_fac; j_fac++) {
      if (f_adj[i_fac + j_fac * n_fac] == 1) {
        ListIdx.push_back(j_fac);
      }
    }
    for (int j_fac = i_fac + 1; j_fac < n_fac; j_fac++)
      ListIdx.push_back(j_fac);
    int n_ineq = ListIdx.size();
    MyMatrix<T> FAC_local(n_ineq, dim);
    for (int i_ineq = 0; i_ineq < n_ineq; i_ineq++) {
      int idx = ListIdx[i_ineq];
      MyVector<T> eFAC = GetMatrixRow(FAC, idx);
      MyVector<T> eFACred = NSP * eFAC;
      FAC_local(i_ineq, 0) = 0;
      for (int i = 0; i < dim - 1; i++)
        FAC_local(i_ineq, i + 1) = eFACred(i);
    }
    std::vector<int> ListIrred = cdd::RedundancyReductionClarkson(FAC_local);
    for (auto &eIrred : ListIrred) {
      int pos = ListIdx[eIrred];
      f_adj[pos + i_fac * n_fac] = 1;
    }
  }
  return f_adj;
}

template <typename T>
LpSolution<T> GLPK_LinearProgramming_Secure(MyMatrix<T> const &ListIneq,
                                            MyVector<T> const &ToBeMinimized) {
  LpSolution<T> TheLP = GLPK_LinearProgramming(ListIneq, ToBeMinimized);
  if (TheLP.method == "cdd")
    return TheLP;
  std::vector<int> ListRowSelect = FaceToVector(TheLP.eFace);
  MyMatrix<T> ListIneq_Incd = SelectRow(ListIneq, ListRowSelect);
  MyVector<T> eVectTest = ToBeMinimized;
  eVectTest(0) -= TheLP.OptimalValue;
  std::optional<MyVector<T>> opt =
      SolutionMatNonnegative(ListIneq_Incd, eVectTest);
  if (!opt)
    return CDD_LinearProgramming(ListIneq, ToBeMinimized);
  return TheLP;
}

template <typename T>
MyMatrix<T> KernelLinearDeterminedByInequalities(MyMatrix<T> const &FAC) {
  int nbCol = FAC.cols();
  int nbRow = FAC.rows();
  PosRelRes<T> eRes = SearchPositiveRelationSimple(FAC);
  std::cerr << "eRes.eTestExist=" << eRes.eTestExist << "\n";
  if (!eRes.eTestExist) {
    return IdentityMat<T>(nbCol);
  } else {
    std::vector<int> ListIdx;
    for (int iRow = 0; iRow < nbRow; iRow++)
      if (eRes.TheRelat(iRow) > 0)
        ListIdx.push_back(iRow);
    MyMatrix<T> FACred = SelectRow(FAC, ListIdx);
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(FACred);
    if (eSelect.NSP.rows() == 0)
      return MyMatrix<T>(0, nbCol);
    MyMatrix<T> FACproj = FAC * eSelect.NSP.transpose();
    MyMatrix<T> FACprojCor = SelectNonZeroRows(FACproj);
    MyMatrix<T> TheSpann;
    if (FACprojCor.rows() == 0)
      TheSpann = IdentityMat<T>(eSelect.NSP.rows());
    else
      TheSpann = KernelLinearDeterminedByInequalities(FACprojCor);
    if (TheSpann.rows() == 0)
      return MyMatrix<T>(0, nbCol);
    else
      return TheSpann * eSelect.NSP;
  }
}

template <typename T> bool IsFullDimensional(MyMatrix<T> const &FAC) {
  PosRelRes<T> eRes = SearchPositiveRelationSimple(FAC);
  if (!eRes.eTestExist) {
    return true;
  } else {
    return false;
  }
}

template <typename T>
MyMatrix<T> LinearDeterminedByInequalities(MyMatrix<T> const &FAC) {
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  return KernelLinearDeterminedByInequalities(SelectNonZeroRows(FAC));
}

/* Finds an interior point in the cone determined by the inequalities.
   No group used here nor any equalities.
 */
template <typename T>
MyVector<T> GetSpaceInteriorPoint_Basic(MyMatrix<T> const &FAC) {
  int n_rows = FAC.rows();
  int n_cols = FAC.cols();
  MyMatrix<T> ListInequalities = ZeroMatrix<T>(n_rows, n_cols + 1);
  MyVector<T> ToBeMinimized = ZeroVector<T>(n_cols + 1);
  for (int i_row = 0; i_row < n_rows; i_row++) {
    ListInequalities(i_row, 0) = -1;
    for (int i_col = 0; i_col < n_cols; i_col++)
      ListInequalities(i_row, i_col + 1) = FAC(i_row, i_col);
    for (int i_col = 0; i_col <= n_cols; i_col++)
      ToBeMinimized(i_col) += ListInequalities(i_row, i_col);
  }
  LpSolution<T> eSol = CDD_LinearProgramming(ListInequalities, ToBeMinimized);
  if (!eSol.PrimalDefined || !eSol.DualDefined) {
    std::cerr << "Failed to find an interior point by linear programming\n";
    std::cerr << "Maybe the cone is actually not full dimensional\n";
    throw TerminalException{1};
  }
  MyVector<T> eVect = eSol.DirectSolution;
#ifdef DEBUG_LINEAR_PROGRAM
  MyVector<T> ListScal = FAC * eVect;
  for (int i_row = 0; i_row < n_rows; i_row++) {
    T eScal = ListScal(i_row);
    if (eScal <= 0) {
      std::cerr << "We have eScal=" << eScal << " which should be positive\n";
      throw TerminalException{1};
    }
  }
#endif
  return eVect;
}

template <typename T>
MyMatrix<T> GetSpaceInteriorPoint(MyMatrix<T> const &FAC,
                                  MyMatrix<T> const &Equa) {
  MyMatrix<T> NSP = NullspaceMat(TransposedMat(Equa));
  MyMatrix<T> FACred = FAC * NSP.transpose();
  MyVector<T> eVectInt = GetSpaceInteriorPoint_Basic(FACred);
  MyVector<T> TheSol = NSP.transpose() * eVectInt;
#ifdef DEBUG_LINEAR_PROGRAM
  MyVector<T> ListScal_FAC = FAC * TheSol;
  for (int i_row = 0; i_row < FAC.rows(); i_row++) {
    T eScal = ListScal_FAC(i_row);
    if (eScal <= 0) {
      std::cerr << "We have eScal=" << eScal << " which should be positive\n";
      throw TerminalException{1};
    }
  }
  MyVector<T> ListScal_Equa = Equa * TheSol;
  for (int i_row = 0; i_row < Equa.rows(); i_row++) {
    T eScal = ListScal_Equa(i_row);
    if (eScal != 0) {
      std::cerr << "We have eScal=" << eScal << " which should be zero\n";
      throw TerminalException{1};
    }
  }
#endif
  return TheSol;
}

template <typename T>
MyMatrix<T> GetSpaceInteriorPointFace(MyMatrix<T> const &FAC, Face const &f) {
  int n_row = FAC.rows();
  int n = FAC.cols();
  int n_f = f.size();
  if (n_row != n_f) {
    std::cerr << "FAC and f have incoherent lengths\n";
    throw TerminalException{1};
  }
  int cnt = f.count();
  MyMatrix<T> FACred(n_row - cnt, n);
  MyMatrix<T> Equa(cnt, n);
  int pos_fac = 0;
  int pos_equa = 0;
  for (int i_row = 0; i_row < n_row; i_row++) {
    if (f[i_row] == 0) {
      FACred.row(pos_fac) = FAC.row(i_row);
      pos_fac++;
    } else {
      Equa.row(pos_equa) = FAC.row(i_row);
      pos_equa++;
    }
  }
  return GetSpaceInteriorPoint(FACred, Equa);
}

template <typename T> struct EmbeddedPolytope {
  MyMatrix<T> LinSpace;
  MyMatrix<T> ListIneq;
};

template <typename T>
EmbeddedPolytope<T> ComputeEmbeddedPolytope(MyMatrix<T> const &ListIneq,
                                            MyMatrix<T> const &ListEqua) {
  MyMatrix<T> NSP = NullspaceTrMat(ListEqua);
  MyMatrix<T> ListIneqRed = ListIneq * (NSP.transpose());
  MyMatrix<T> LinSpa = LinearDeterminedByInequalities(ListIneqRed);
  MyMatrix<T> FinalSpace = (LinSpa.transpose()) * NSP;
  std::optional<MyMatrix<T>> eRes = AffinizeSubspace(FinalSpace);
  if (!eRes) {
    std::cerr << "Call to ComputeEmbeddedPolytope failed\n";
    std::cerr << "Because the affinization operation failed\n";
    throw TerminalException{1};
  }
  MyMatrix<T> FinalSpace_Aff = *eRes;
  MyMatrix<T> ListIneqFinal = ListIneq * (FinalSpace_Aff.transpose());
  return {FinalSpace_Aff, ListIneqFinal};
}

// clang-format off
#endif  // SRC_POLY_POLY_LINEARPROGRAMMING_H_
// clang-format on
