#ifndef SRC_POLY_POLY_LINEARPROGRAMMING_GLPK_H_
#define SRC_POLY_POLY_LINEARPROGRAMMING_GLPK_H_

#include "POLY_PolytopeFct.h"
#include <glpk.h>
#include <string>
#include <vector>

/*
  GLPK does not return the unfeasibility proof as opposed to CDD.
  We may have to write an auxilliary program for that.
  It also uses real for discriminating the solution.
  ----
  equation xA >= Acst
           xB = Bcst
  optimize x MVect
 */
template <typename T>
LpSolutionSimple<T> GLPK_LinearProgramming_Kernel_Sparse_PROC(
    MySparseMatrix<T> const &Aspmat, MyVector<T> const &ListAconst,
    MySparseMatrix<T> const &Bspmat, MyVector<T> const &ListBconst,
    MyVector<T> const &ToBeMinimized, GLPKoption const &eGLPKoption) {
  int nbVar = Aspmat.cols();
  int nbIneq = Aspmat.rows();
  int nbEqua = Bspmat.rows();
  int nbConstIneq = ListAconst.rows();
  if (nbIneq != nbConstIneq) {
    std::cerr << "We have inconsistency in input\n";
    std::cerr << "nbConstIneq = " << nbConstIneq << "\n";
    std::cerr << "     nbIneq = " << nbIneq << "\n";
    throw TerminalException{1};
  }
  int nbConstEqua = ListBconst.rows();
  if (nbEqua != nbConstEqua) {
    std::cerr << "We have inconsistency in input\n";
    std::cerr << "nbConstEqua = " << nbConstEqua << "\n";
    std::cerr << "     nbEqua = " << nbEqua << "\n";
    throw TerminalException{1};
  }
  std::string eRand = random_string(20);
  std::string ePrefix = "/tmp/LP_" + eRand + "/";
  CreateDirectory(ePrefix);
  //
  std::string FileMath = ePrefix + "GLP.mod";
  std::string FileOut = ePrefix + "GLP.out";
  std::string FileErr1 = ePrefix + "GLP.err1";
  std::string FileErr2 = ePrefix + "GLP.err2";
  std::function<void(void)> CleanAtLeaving = [&](void) -> void {
    /*
    RemoveFileIfExist(FileMath);
    RemoveFileIfExist(FileOut);
    RemoveFileIfExist(FileErr1);
    RemoveFileIfExist(FileErr2);*/
  };
  //
  // Definition of needed data types and routines
  //
  std::function<void(std::ostream &, T)> PrintValue =
      [&](std::ostream &os, T const &eVal) -> void {
    if (eGLPKoption.UseDouble) {
      double eVal_d = UniversalScalarConversion<double, T>(eVal);
      os << eVal_d;
    } else {
      os << eVal;
    }
  };
  struct PairCV {
    int iCol;
    T eVal;
  };
  std::function<void(std::ostream &, std::vector<PairCV>)>
      PrintConstraintVector =
          [&](std::ostream &os, std::vector<PairCV> const &LPair) -> void {
    int len = LPair.size();
    bool IsFirst = true;
    for (int i = 0; i < len; i++) {
      int iVar = 1 + LPair[i].iCol;
      T eVal = LPair[i].eVal;
      if (eVal != 0) {
        if (IsFirst) {
          if (eVal == 1) {
            os << " x" << iVar;
          } else {
            if (eVal == -1) {
              os << " - x" << iVar;
            } else {
              if (eVal < 0) {
                os << " - ";
                PrintValue(os, -eVal);
                os << " * x" << iVar;
              } else {
                os << " ";
                PrintValue(os, eVal);
                os << " * x" << iVar;
              }
            }
          }
        } else {
          if (eVal == 1) {
            os << " + x" << iVar;
          } else {
            if (eVal == -1) {
              os << " - x" << iVar;
            } else {
              if (eVal < 0) {
                os << " - ";
                PrintValue(os, -eVal);
                os << " * x" << iVar;
              } else {
                os << " + ";
                PrintValue(os, eVal);
                os << " * x" << iVar;
              }
            }
          }
        }
        IsFirst = false;
      }
    }
  };
  //
  // Now printing the model description
  //
  std::ofstream OUTfs;
  OUTfs.open(FileMath);
  OUTfs << "# GLP_IntegerLinearProgramming\n";
  for (int iVar = 1; iVar <= nbVar; iVar++) {
    OUTfs << "var x" << iVar << ";\n";
  }
  OUTfs << "minimize obj:";
  std::vector<PairCV> ToBeMinVect(nbVar);
  for (int iVar = 0; iVar < nbVar; iVar++) {
    PairCV ePair{iVar, ToBeMinimized[iVar]};
    ToBeMinVect[iVar] = ePair;
  }
  PrintConstraintVector(OUTfs, ToBeMinVect);
  OUTfs << ";\n";
  //
  std::function<void(std::vector<PairCV> &, T &)> EquaCanonicalize =
      [&](std::vector<PairCV> &ListPair, T &eVal) -> void {
    int lenA = ListPair.size();
    int len = lenA + 1;
    MyVector<T> eVect(len);
    eVect(0) = eVal;
    for (int i = 0; i < lenA; i++)
      eVect(i + 1) = ListPair[i].eVal;
    MyVector<T> eVectNew = RemoveFractionVector(eVect);
    //
    eVal = eVectNew(0);
    for (int i = 0; i < lenA; i++)
      ListPair[i].eVal = eVectNew(i + 1);
  };
  int iConst = 0;
  std::vector<std::vector<PairCV>> LLPair(nbIneq);
  for (int k = 0; k < Aspmat.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(Aspmat, k); it; ++it) {
      T eVal = it.value();
      int iRow = it.row();
      int iCol = it.col();
      PairCV ePair{iCol, eVal};
      LLPair[iRow].push_back(ePair);
    }
  std::function<LpSolutionSimple<T>(void)> StandardUnfeas =
      [&](void) -> LpSolutionSimple<T> {
    CleanAtLeaving();
    int nbRow = -1;
    int nbCol = -1;
    MyVector<int> ColumnStatus(1);
    MyVector<int> RowStatus(1);
    return {false, 0, nbRow, nbCol, {}, {}, ColumnStatus, RowStatus};
  };
  for (int iIneq = 0; iIneq < nbIneq; iIneq++) {
    iConst++;
    T SumAbs = 0;
    for (auto &ePair : LLPair[iIneq])
      SumAbs += T_abs(ePair.eVal);
    if (SumAbs == 0 && ListAconst[iIneq] > 0) {
      return StandardUnfeas();
    } else {
      if (SumAbs > 0) {
        T eRhs = ListAconst[iIneq];
        EquaCanonicalize(LLPair[iIneq], eRhs);
        //
        OUTfs << "s.t. c" << iConst << ":";
        PrintConstraintVector(OUTfs, LLPair[iIneq]);
        OUTfs << " >= ";
        PrintValue(OUTfs, eRhs);
        OUTfs << ";\n";
      }
    }
  }
  std::vector<std::vector<PairCV>> LLPairEqua(nbEqua);
  for (int k = 0; k < Bspmat.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(Bspmat, k); it; ++it) {
      T eVal = it.value();
      int iRow = it.row();
      int iCol = it.col();
      PairCV ePair{iCol, eVal};
      LLPairEqua[iRow].push_back(ePair);
    }
  for (int iEqua = 0; iEqua < nbEqua; iEqua++) {
    iConst++;
    T SumAbs = 0;
    for (auto &ePair : LLPairEqua[iEqua])
      SumAbs += T_abs(ePair.eVal);
    if (SumAbs == 0 && ListBconst[iEqua] != 0) {
      return StandardUnfeas();
    } else {
      if (SumAbs > 0) {
        T eRhs = ListBconst[iEqua];
        EquaCanonicalize(LLPairEqua[iEqua], eRhs);
        //
        OUTfs << "s.t. c" << iConst << ":";
        PrintConstraintVector(OUTfs, LLPairEqua[iEqua]);
        OUTfs << " = ";
        PrintValue(OUTfs, eRhs);
        OUTfs << ";\n";
      }
    }
  }
  OUTfs << "solve;\n";
  OUTfs << "display\n";
  for (int iVar = 1; iVar <= nbVar; iVar++) {
    if (iVar > 1)
      OUTfs << ",";
    OUTfs << " x" << iVar;
  }
  OUTfs << ";\n";
  OUTfs << "end;\n";
  OUTfs.close();
  //
  // Now running the GLPK program
  //
  std::string eCommand = "glpsol";
  if (eGLPKoption.UseExact) {
    eCommand += " --exact";
  }
  if (eGLPKoption.UseXcheck) {
    eCommand += " --xcheck";
  }
  eCommand += " --output " + FileOut + " --math " + FileMath;
  eCommand += " > " + FileErr1 + " 2> " + FileErr2;
  std::cerr << "eCommand=" << eCommand << "\n";
  int iret = system(eCommand.c_str());
  if (iret == -1) {
    printf("Oh dear, something went wrong with glpsol! %s\n", strerror(errno));
    throw TerminalException{1};
  }
  if (!IsExistingFile(FileErr1) || !IsExistingFile(FileErr2) ||
      !IsExistingFile(FileOut)) {
    std::cerr << "Not all the files of the computation have been created\n";
    std::cerr << "FileErr1 = " << FileErr1 << "\n";
    std::cerr << "FileErr2 = " << FileErr2 << "\n";
    std::cerr << "FileOut = " << FileOut << "\n";
    throw TerminalException{1};
  }
  //
  // Now reading the data output
  //
  std::vector<std::string> RESUL;
  std::ifstream INfs(FileOut);
  std::string line;
  while (getline(INfs, line))
    RESUL.push_back(line);
  INfs.close();
  int nbLine = RESUL.size();
  //  std::cerr << "nbLine=" << nbLine << "\n";
  //
  bool IsUndefined = false;
  int nbRow = -1;
  int nbCol = -1;
  for (int iLine = 0; iLine < nbLine; iLine++) {
    std::vector<std::string> LSplit = STRING_Split(RESUL[iLine], " ");
    int len = LSplit.size();
    if (len > 1) {
      if (LSplit[0] == "Status:" && LSplit[1] == "UNDEFINED")
        IsUndefined = true;
      if (LSplit[0] == "Rows:") {
        nbRow = atoi(LSplit[1].c_str());
      }
      if (LSplit[0] == "Columns:") {
        nbCol = atoi(LSplit[1].c_str());
      }
    }
  }
  if (nbRow == -1 || nbCol == -1) {
    std::cerr << "We have nbRow=" << nbRow << " and nbCol=" << nbCol << "\n";
    std::cerr << "This is inconsistent\n";
    throw TerminalException{1};
  }
  MyVector<int> RowStatus(nbRow);
  MyVector<int> ColumnStatus(nbCol);
  int iLineCritColumn = -1;
  int iLineCritRow = -1;
  //  std::cerr << "Before determination of iLineCritS\n";
  for (int iLine = 0; iLine < nbLine; iLine++) {
    std::vector<std::string> Ucol = STRING_Split(RESUL[iLine], "Column name");
    if (Ucol.size() > 1 && iLineCritColumn == -1)
      iLineCritColumn = iLine;
    std::vector<std::string> Urow = STRING_Split(RESUL[iLine], "Row name");
    if (Urow.size() > 1 && iLineCritRow == -1)
      iLineCritRow = iLine;
  }
  if (iLineCritColumn == -1 || iLineCritRow == -1) {
    std::cerr << "We did not find iLineCritColumn or iLineCritRow\n";
    throw TerminalException{1};
  }
  //  std::cerr << "Before nbCol iteration\n";
  for (int iCol = 0; iCol < nbCol; iCol++) {
    std::vector<std::string> U =
        STRING_Split(RESUL[iLineCritColumn + 2 + iCol], " ");
    std::string eStat = U[2];
    int eVal = 0;
    if (eStat == "B")
      eVal = 1;
    if (eStat == "NF")
      eVal = 2;
    if (eStat == "NL")
      eVal = 3;
    if (eStat == "NS")
      eVal = 4;
    if (eVal == 0) {
      std::cerr << "eStat=" << eStat << "\n";
      std::cerr << "Failed to assign eVal for columns\n";
      throw TerminalException{1};
    }
    ColumnStatus[iCol] = eVal;
  }
  //  std::cerr << "Before nbRow iteration\n";
  for (int iRow = 0; iRow < nbRow; iRow++) {
    std::vector<std::string> U =
        STRING_Split(RESUL[iLineCritRow + 2 + iRow], " ");
    std::string eStat = U[2];
    int eVal = 0;
    if (eStat == "B")
      eVal = 1;
    if (eStat == "NF")
      eVal = 2;
    if (eStat == "NL")
      eVal = 3;
    if (eStat == "NS")
      eVal = 4;
    if (eVal == 0) {
      std::cerr << "Failed to assign eVal for rows\n";
      throw TerminalException{1};
    }
    RowStatus[iRow] = eVal;
  }
  if (IsUndefined) {
    CleanAtLeaving();
    return {false, 0, nbRow, nbCol, {}, {}, ColumnStatus, RowStatus};
  }
  MyVector<T> DirectSolution(nbVar);
  MyVector<T> DirectSolutionExt(nbVar + 1);
  DirectSolutionExt(0) = 1;
  //  std::cerr << "Before nbRow iteration\n";
  for (int iVar = 1; iVar <= nbVar; iVar++) {
    //    std::cerr << "iVar=" << iVar << "\n";
    std::string eVar = "x" + IntToString(iVar);
    std::string strBreak = RESUL[iLineCritColumn + 2 + iVar - 1];
    std::vector<std::string> U = STRING_Split(strBreak, eVar);
    if (U.size() == 1) {
      std::cerr << "strBreak=" << strBreak << "\n";
      std::cerr << "Likely inconsistency in the code\n";
      throw TerminalException{1};
    }
    std::string eLine = RESUL[iLineCritColumn + 2 + iVar - 1];
    //    std::cerr << "  eLine=" << eLine << "\n";
    std::vector<std::string> U2 = STRING_Split(eLine, " ");
    std::string eValStr = U2[U2.size() - 2];
    //    std::cerr << "  |U2|=" << U2.size() << "\n";
    //    std::cerr << "  eValStr=" << eValStr << "\n";
    double eVal = atof(eValStr.c_str());
    //    std::cerr << "  eVal=" << eVal << "\n";
    DirectSolution(iVar - 1) = eVal;
  }
  T TheOptimal = 0;
  for (int iVar = 0; iVar < nbVar; iVar++)
    TheOptimal += DirectSolution(iVar) * ToBeMinimized(iVar);
  CleanAtLeaving();
  return {true,           TheOptimal,        nbRow,        nbCol,
          DirectSolution, DirectSolutionExt, ColumnStatus, RowStatus};
}

LpSolutionSimple<double> GLPK_LinearProgramming_Kernel_Sparse_LIBRARY(
    MySparseMatrix<double> const &Aspmat, MyVector<double> const &ListAconst,
    MySparseMatrix<double> const &Bspmat,
    [[maybe_unused]] MyVector<double> const &ListBconst,
    MyVector<double> const &ToBeMinimized,
    [[maybe_unused]] GLPKoption const &eGLPKoption) {
  glp_prob *prob = NULL;
  //  glp_erase_prob(prob);
  prob = glp_create_prob();
  glp_bfcp *bfcp = NULL;
  glp_set_bfcp(prob, bfcp);
  int nbVar = Aspmat.cols();
  int nbIneq = Aspmat.rows();
  int nbRow = nbIneq;
  int nbCol = nbVar;
  MyVector<int> RowStatus(nbIneq);
  MyVector<int> ColumnStatus(nbVar);
  int nbEqua = Bspmat.rows();
  if (nbEqua > 0) {
    std::cerr << "Right now that code does not handle equations\n";
    throw TerminalException{1};
  }
  struct PairCV {
    int iCol;
    double eVal;
  };
  int m = nbIneq + 1;
  glp_add_rows(prob, m);
  for (int i = 1; i <= m; i++) {
    double lb, ub = 0;
    if (i == 1)
      lb = 0;
    else
      lb = ListAconst(i - 2);
    glp_set_row_bnds(prob, i, GLP_LO, lb, ub);
  }
  glp_add_cols(prob, nbVar);
  for (int j = 1; j <= nbVar; j++) {
    glp_set_col_kind(prob, j, GLP_CV);
    double lb = 0, ub = 0;
    glp_set_col_bnds(prob, j, GLP_FR, lb, ub);
  }
  std::vector<std::vector<PairCV>> LLPairIneq(nbIneq);
  for (int k = 0; k < Aspmat.outerSize(); ++k)
    for (typename MySparseMatrix<double>::InnerIterator it(Aspmat, k); it;
         ++it) {
      double eVal = it.value();
      int iRow = it.row();
      int iCol = it.col();
      PairCV ePair{iCol, eVal};
      LLPairIneq[iRow].push_back(ePair);
    }
  for (int i = 1; i <= m; i++) {
    std::vector<int> ind{0};
    std::vector<double> val{double(0)};
    if (i > 1) {
      for (auto &ePair : LLPairIneq[i - 2]) {
        ind.push_back(ePair.iCol + 1);
        val.push_back(ePair.eVal);
      }
    }
    int len = ind.size() - 1;
    glp_set_mat_row(prob, i, len, ind.data(), val.data());
  }
  glp_set_obj_dir(prob, GLP_MIN);
  glp_set_obj_coef(prob, 0, double(0));
  for (int i = 1; i <= nbVar; i++) {
    double eVal = ToBeMinimized(i - 1);
    glp_set_obj_coef(prob, i, eVal);
  }
  glp_smcp eSmcp;
  glp_init_smcp(&eSmcp);
  eSmcp.presolve = GLP_ON;
  eSmcp.r_test = GLP_RT_STD;
  glp_simplex(prob, &eSmcp);
  if (glp_get_status(prob) != GLP_OPT) {
    return {false, double(0), nbRow, nbCol, {}, {}, ColumnStatus, RowStatus};
  }

  std::cerr << "DirectSolution, step 1\n";
  MyVector<double> DirectSolution(nbVar);
  std::cerr << "DirectSolution, step 2\n";
  MyVector<double> DirectSolutionExt(nbVar + 1);
  std::cerr << "DirectSolution, step 3\n";
  DirectSolutionExt(0) = 1;
  for (int i = 1; i <= nbVar; i++) {
    double eVal = glp_get_col_prim(prob, i);
    DirectSolution(i - 1) = eVal;
    DirectSolutionExt(i) = eVal;
  }
  auto ConversionValue = [&](int val) -> int {
    if (val == GLP_BS)
      return 1;
    if (val == GLP_NF)
      return 2;
    if (val == GLP_NL)
      return 3;
    if (val == GLP_NS)
      return 4;
    if (val == GLP_NU)
      return 5;
    return -1;
  };
  for (int iIneq = 0; iIneq < nbIneq; iIneq++) {
    int val = glp_get_row_stat(prob, iIneq + 2);
    RowStatus(iIneq) = ConversionValue(val);
  }
  for (int iCol = 0; iCol < nbVar; iCol++) {
    int val = glp_get_col_stat(prob, iCol + 1);
    ColumnStatus(iCol) = ConversionValue(val);
  }
  double TheOptimal = glp_get_obj_val(prob);

  glp_erase_prob(prob);
  return {true,           TheOptimal,        nbRow,        nbCol,
          DirectSolution, DirectSolutionExt, ColumnStatus, RowStatus};
}

template <typename T> struct LinProgSparseDecomp {
  MySparseMatrix<T> Aspmat;
  MySparseMatrix<T> Bspmat;
  MyVector<T> ListAconst;
  MyVector<T> ListBconst;
};

template <typename T>
LinProgSparseDecomp<T> GetSparseDecomposition(MyMatrix<T> const &ListEqua,
                                              MyMatrix<T> const &ListIneq) {
  auto fDecompose = [](MyMatrix<T> const &eMat, MySparseMatrix<T> &eSpMat,
                       MyVector<T> &eV) -> void {
    int nbRow = eMat.rows();
    int nbCol = eMat.cols();
    eV = MyVector<T>(nbRow);
    for (int iRow = 0; iRow < nbRow; iRow++)
      eV(iRow) = -eMat(iRow, 0);
    //
    int nnz = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 1; iCol < nbCol; iCol++)
        if (eMat(iRow, iCol) != 0)
          nnz++;
    using T2 = Eigen::Triplet<T>;
    std::vector<T2> tripletList(nnz);
    int iZ = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 1; iCol < nbCol; iCol++) {
        T eVal = eMat(iRow, iCol);
        if (eVal != 0) {
          tripletList[iZ] = T2(iRow, iCol - 1, eVal);
          iZ++;
        }
      }
    eSpMat = MySparseMatrix<T>(nbRow, nbCol - 1);
    eSpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  };
  MySparseMatrix<T> Aspmat, Bspmat;
  MyVector<T> ListAconst, ListBconst;
  fDecompose(ListEqua, Bspmat, ListBconst);
  fDecompose(ListIneq, Aspmat, ListAconst);
  return {Aspmat, Bspmat, ListAconst, ListBconst};
}

template <typename T>
LinProgSparseDecomp<double>
GetSparseDecompositionDouble(MyMatrix<T> const &ListEqua,
                             MyMatrix<T> const &ListIneq) {
  auto fDecompose = [](MyMatrix<T> const &eMat, MySparseMatrix<double> &eSpMat,
                       MyVector<double> &eV) -> void {
    int nbRow = eMat.rows();
    int nbCol = eMat.cols();
    eV = MyVector<double>(nbRow);
    for (int iRow = 0; iRow < nbRow; iRow++)
      eV(iRow) = UniversalScalarConversion<double, T>(-eMat(iRow, 0));
    //
    int nnz = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 1; iCol < nbCol; iCol++)
        if (eMat(iRow, iCol) != 0)
          nnz++;
    using T2 = Eigen::Triplet<double>;
    std::vector<T2> tripletList(nnz);
    int iZ = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 1; iCol < nbCol; iCol++) {
        T eVal = eMat(iRow, iCol);
        if (eVal != 0) {
          double eVal_d = UniversalScalarConversion<double, T>(eVal);
          tripletList[iZ] = T2(iRow, iCol - 1, eVal_d);
          iZ++;
        }
      }
    eSpMat = MySparseMatrix<double>(nbRow, nbCol - 1);
    eSpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  };
  MySparseMatrix<double> Aspmat, Bspmat;
  MyVector<double> ListAconst, ListBconst;
  fDecompose(ListEqua, Bspmat, ListBconst);
  fDecompose(ListIneq, Aspmat, ListAconst);
  return {Aspmat, Bspmat, ListAconst, ListBconst};
}

template <typename T>
LpSolutionSimple<T> GLPK_LinearProgramming_Kernel_Dense_PROC(
    MyMatrix<T> const &ListEqua, MyMatrix<T> const &ListIneq,
    MyVector<T> const &ToBeMinimized, GLPKoption const &eGLPKoption) {
  LinProgSparseDecomp<T> RecSpDecomp =
      GetSparseDecomposition(ListEqua, ListIneq);
  return GLPK_LinearProgramming_Kernel_Sparse_PROC(
      RecSpDecomp.Aspmat, RecSpDecomp.ListAconst, RecSpDecomp.Bspmat,
      RecSpDecomp.ListBconst, ToBeMinimized, eGLPKoption);
}

template <typename T>
LpSolutionSimple<double> GLPK_LinearProgramming_Kernel_Dense_LIBRARY(
    MyMatrix<T> const &ListEqua, MyMatrix<T> const &ListIneq,
    MyVector<T> const &ToBeMinimized, GLPKoption const &eGLPKoption) {
  LinProgSparseDecomp<double> RecSpDecomp =
      GetSparseDecompositionDouble(ListEqua, ListIneq);
  MyVector<double> ToBeMinimized_d =
      UniversalVectorConversion<double, T>(ToBeMinimized);
  //  return GLPK_LinearProgramming_Kernel_Sparse_LIBRARY(RecSpDecomp.Aspmat,
  //  RecSpDecomp.ListAconst, RecSpDecomp.Bspmat, RecSpDecomp.ListBconst,
  //  ToBeMinimized_d, eGLPKoption);
  return GLPK_LinearProgramming_Kernel_Sparse_PROC(
      RecSpDecomp.Aspmat, RecSpDecomp.ListAconst, RecSpDecomp.Bspmat,
      RecSpDecomp.ListBconst, ToBeMinimized_d, eGLPKoption);
}

#endif //  SRC_POLY_POLY_LINEARPROGRAMMING_GLPK_H_
