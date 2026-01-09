// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_TSPACE_LISTMATTRANSFORM_H_
#define SRC_LATT_TSPACE_LISTMATTRANSFORM_H_

template <typename T>
MyMatrix<T> GetMatrixFromBasis(std::vector<MyMatrix<T>> const &ListMat,
                               MyVector<T> const &eVect) {
  int n = ListMat[0].rows();
  MyMatrix<T> RetMat = ZeroMatrix<T>(n, n);
  int nbMat = ListMat.size();
  int siz = eVect.size();
  if (siz != nbMat) {
    std::cerr << "TSPACE: Error in GetMatrixFromBasis\n";
    std::cerr << "TSPACE:  siz=" << siz << "\n";
    std::cerr << "TSPACE: nbMat=" << nbMat << "\n";
    std::cerr << "TSPACE: But they should be both equal\n";
    throw TerminalException{1};
  }
  for (int iMat = 0; iMat < nbMat; iMat++) {
    RetMat += eVect(iMat) * ListMat[iMat];
  }
  return RetMat;
}

//
// We search for the set of matrices satisfying g M g^T = M for all g in ListGen
//
template <typename T>
std::vector<MyMatrix<T>>
BasisInvariantForm(int const &n, std::vector<MyMatrix<T>> const &ListGen,
                   [[maybe_unused]] std::ostream &os) {
  std::vector<std::vector<int>> ListCoeff;
  MyMatrix<int> ListCoeffRev(n, n);
  int pos = 0;
  for (int iLin = 0; iLin < n; iLin++) {
    for (int iCol = 0; iCol <= iLin; iCol++) {
      ListCoeff.push_back({iLin, iCol});
      ListCoeffRev(iLin, iCol) = pos;
      ListCoeffRev(iCol, iLin) = pos;
      pos++;
    }
  }
  int nbCoeff = pos;
  auto FuncPos = [&](int const &i, int const &j) -> int {
    return ListCoeffRev(i, j);
  };
  std::vector<MyVector<T>> ListEquations;
  for (auto &eGen : ListGen) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j <= i; j++) {
        MyVector<T> TheEquation = ZeroVector<T>(nbCoeff);
        TheEquation(FuncPos(i, j)) += 1;
        for (int k = 0; k < n; k++) {
          for (int l = 0; l < n; l++) {
            int pos = FuncPos(k, l);
            TheEquation(pos) += -eGen(i, k) * eGen(j, l);
          }
        }
        ListEquations.push_back(TheEquation);
      }
    }
  }
  MyMatrix<T> MatEquations = MatrixFromVectorFamily(ListEquations);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: Before call to NullspaceTrMat nbGen=" << ListGen.size()
     << " MatEquations(rows/cols)=" << MatEquations.rows() << " / "
     << MatEquations.cols() << "\n";
#endif
  MyMatrix<T> NSP = NullspaceTrMat(MatEquations);
#ifdef DEBUG_TSPACE_FUNCTIONS
  os << "TSPACE: After call to NullspaceTrMat |NSP|=" << NSP.rows() << " / "
     << NSP.cols() << "\n";
#endif
  int dimSpa = NSP.rows();
  std::vector<MyMatrix<T>> TheBasis(dimSpa);
  for (int iDim = 0; iDim < dimSpa; iDim++) {
    MyVector<T> eRow = GetMatrixRow(NSP, iDim);
    MyMatrix<T> eMat(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        eMat(i, j) = eRow(FuncPos(i, j));
    TheBasis[iDim] = eMat;
  }
#ifdef DEBUG_TSPACE_FUNCTIONS
  for (auto &eBasis : TheBasis) {
    for (auto &eGen : ListGen) {
      MyMatrix<T> eProd = eGen * eBasis * eGen.transpose();
      if (eProd != eBasis) {
        std::cerr
            << "TSPACE: We have eProd <> eBasis, so a linear algebra bug\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  return TheBasis;
}

template <typename T, typename Tint>
std::vector<MyMatrix<T>> get_spec_basis_shv_equal(std::vector<MyMatrix<T>> const& ListMat, std::vector<MyVector<Tint>> const &ListVect) {
  int nbVect = ListVect.size();
  int nbMat = ListMat.size();
  MyMatrix<T> MatrixValuesDiff(nbVect - 1, nbMat);
  for (int iMat = 0; iMat < nbMat; iMat++) {
    MyMatrix<T> const& eMat = ListMat[iMat];
    MyVector<Tint> const& V1 = ListVect[0];
    T val1 = EvaluationQuadForm(eMat, V1);
    for (int iVect = 0; iVect < nbVect - 1; iVect++) {
      MyVector<Tint> const& V2 = ListVect[iVect+1];
      T val2 = EvaluationQuadForm(eMat, V2);
      MatrixValuesDiff(iVect, iMat) = val2 - val1;
    }
  }
  MyMatrix<T> NSP = NullspaceTrMat(MatrixValuesDiff);
  int dimSpa = NSP.rows();
  std::vector<MyMatrix<T>> TheBasis(dimSpa);
  for (int iDim = 0; iDim < dimSpa; iDim++) {
    MyVector<T> eRow = GetMatrixRow(NSP, iDim);
    TheBasis[iDim] = GetMatrixFromBasis(ListMat, eRow);
  }
  return TheBasis;
}

template <typename T, typename Tint>
std::vector<MyMatrix<T>> get_spec_basis_shv_zero(std::vector<MyMatrix<T>> const& ListMat, std::vector<MyVector<Tint>> const &ListVect) {
  int nbVect = ListVect.size();
  int nbMat = ListMat.size();
  MyMatrix<T> MatrixValues(nbVect, nbMat);
  for (int iMat = 0; iMat < nbMat; iMat++) {
    MyMatrix<T> const& eMat = ListMat[iMat];
    for (int iVect = 0; iVect < nbVect; iVect++) {
      MyVector<Tint> const& V = ListVect[iVect];
      T val = EvaluationQuadForm(eMat, V);
      MatrixValues(iVect, iMat) = val;
    }
  }
  MyMatrix<T> NSP = NullspaceTrMat(MatrixValues);
  int dimSpa = NSP.rows();
  std::vector<MyMatrix<T>> TheBasis(dimSpa);
  for (int iDim = 0; iDim < dimSpa; iDim++) {
    MyVector<T> eRow = GetMatrixRow(NSP, iDim);
    TheBasis[iDim] = GetMatrixFromBasis(ListMat, eRow);
  }
  return TheBasis;
}


template <typename T, typename Tint>
std::vector<MyMatrix<T>> get_spec_basis_shv_kernel(std::vector<MyMatrix<T>> const& ListMat, std::vector<MyVector<Tint>> const &ListVect) {
  int nbVect = ListVect.size();
  int nbMat = ListMat.size();
  MyMatrix<T> MatrixScalars(nbVect * nbVect, nbMat);
  for (int iMat = 0; iMat < nbMat; iMat++) {
    MyMatrix<T> const& eMat = ListMat[iMat];
    int pos = 0;
    for (int iVect = 0; iVect < nbVect; iVect++) {
      MyVector<Tint> const& V0 = ListVect[iVect];
      for (int jVect = 0; jVect < nbVect; jVect++) {
        MyVector<Tint> const& V1 = ListVect[jVect];
        T val = ScalarProductQuadForm(eMat, V0, V1);
        MatrixScalars(pos, iMat) = val;
        pos += 1;
      }
    }
  }
  MyMatrix<T> NSP = NullspaceTrMat(MatrixScalars);
  int dimSpa = NSP.rows();
  std::vector<MyMatrix<T>> TheBasis(dimSpa);
  for (int iDim = 0; iDim < dimSpa; iDim++) {
    MyVector<T> eRow = GetMatrixRow(NSP, iDim);
    TheBasis[iDim] = GetMatrixFromBasis(ListMat, eRow);
  }
  return TheBasis;
}

// clang-format off
#endif  // SRC_LATT_TSPACE_LISTMATTRANSFORM_H_
// clang-format on
