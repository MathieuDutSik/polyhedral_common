#ifndef TEMP_SHORT_VECTOR_Undefinite_H
#define TEMP_SHORT_VECTOR_Undefinite_H

#include "MAT_Matrix.h"
#include "MAT_MatrixFLT.h"
#include "Temp_Positivity.h"

template <typename Tint, typename T, typename Tfloat>
std::optional<MyVector<Tint>>
GetShortVector_unlimited_float_kernel(MyMatrix<T> const &M, T const &CritNorm,
                                      bool const &StrictIneq,
                                      bool const &NeedNonZero) {
  int n = M.rows();
  MyMatrix<Tfloat> M_f(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      Tfloat eVal_f = UniversalScalarConversion<Tfloat, T>(M(i, j));
      M_f(i, j) = eVal_f;
    }
  MyMatrix<Tfloat> ListEigVect(n, n);
  MyVector<Tfloat> ListEigVal(n);
  jacobi_double(M_f, ListEigVal, ListEigVect);
  int nbNeg = 0;
  for (int i = 0; i < n; i++) {
    Tfloat sum = 0;
    for (int j = 0; j < n; j++)
      sum += T_abs(ListEigVect(i, j));
    for (int j = 0; j < n; j++)
      ListEigVect(i, j) /= sum;
    if (ListEigVal(i) < 0)
      nbNeg++;
  }
  /*
  std::cerr << "ListEigVal=\n";
  WriteVector(std::cerr, ListEigVal);
  std::cerr << "ListEigVect=\n";
  WriteMatrix(std::cerr, ListEigVect);*/

  if (nbNeg == 0) {
    return {};
  }
  int eMult = 1;
  while (true) {
    //    std::cerr << "eMult=" << eMult << "\n";
    for (int i = 0; i < n; i++) {
      MyVector<Tint> eVect(n);
      bool IsZero = true;
      /*      std::cerr << "  i=" << i << "    ListEigVect(i,:)=";
      for (int j=0; j<n; j++) {
        std::cerr << " " << ListEigVect(i,j);
      }
      std::cerr << "\n";*/
      for (int j = 0; j < n; j++) {
        Tfloat eVal_f = eMult * ListEigVect(i, j);
        Tint eVal = UniversalNearestScalarInteger<Tint, Tfloat>(eVal_f);
        //	std::cerr << "j=" << j << "  eVal_f=" << eVal_f << "  eVal=" <<
        //eVal << "\n";
        if (eVal != 0)
          IsZero = false;
        eVect(j) = eVal;
      }
      /*
      std::cerr << "\n";
      std::cerr << "    eVect=";
      for (int j=0; j<n; j++)
        std::cerr << " " << eVect(j);
        std::cerr << "\n";*/
      if (!IsZero || !NeedNonZero) {
        T norm = EvaluationQuadForm(M, eVect);
        /*
        std::cerr << "    norm=" << norm << " CritNorm=" << CritNorm << "\n";
        std::cerr << "    StrictIneq=" << StrictIneq << "\n";
        if (norm < CritNorm) {
          std::cerr << "    We have norm < CritNorm\n";
        }
        else {
          std::cerr << "    We DO NOT have norm < CritNorm\n";
        }
        if (norm <= CritNorm) {
          std::cerr << "    We have norm <= CritNorm\n";
        }
        else {
          std::cerr << "    We DO NOT have norm <= CritNorm\n";
          }*/
        if ((!StrictIneq && norm <= CritNorm) || norm < CritNorm) {
          return eVect;
        }
      }
    }
    eMult++;
  }
}

template <typename Tint, typename T>
MyVector<Tint> GetShortVector_unlimited_float(MyMatrix<T> const &M,
                                              T const &CritNorm,
                                              bool const &StrictIneq,
                                              bool const &NeedNonZero) {
  //  std::cerr << "Before call to GetShortVector_unlimited_float_kernel\n";
  std::optional<MyVector<Tint>> opt =
      GetShortVector_unlimited_float_kernel<Tint, T, double>(
          M, CritNorm, StrictIneq, NeedNonZero);
  //  std::cerr << " After call to GetShortVector_unlimited_float_kernel\n";
  if (opt)
    return *opt;
#ifdef USE_MPREAL
  int theprec = 20;
  while (true) {
    std::cerr << "theprec=" << theprec << "\n";
    mpfr::mpreal::set_default_prec(theprec);
    std::optional<MyVector<Tint>> opt =
        GetShortVector_unlimited_float_kernel<Tint, T, mpfr::mpreal>(
            M, CritNorm, StrictIneq, NeedNonZero);
    if (optt)
      return *opt;
    theprec += 8;
  }
#else
  std::cerr << "mpreal might be needed after all\n";
  throw TerminalException{1};
#endif
}

#endif
