// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_SHORT_SHORT_ENUMERATION_H_
#define SRC_SHORT_SHORT_ENUMERATION_H_

// Canonicalize the vector under the action of the symmetric group
// and inversion (group size 2^n n!)
//
MyVector<int> CyclicCanonicalization_SymN(MyVector<int> const &fCand,
                                          int const &d) {
  int n = fCand.size();
  std::vector<int> RetCand;
  for (int i = 0; i < n; i++) {
    int aVal = fCand(i);
    int bVal;
    if (aVal < 0) {
      bVal = (-aVal) % d;
    } else {
      bVal = aVal % d;
    }
    int cVal;
    if (2 * bVal > d) {
      cVal = d - bVal;
    } else {
      cVal = bVal;
    }
    RetCand.push_back(cVal);
  }
  int len = fCand.size();
  std::vector<int> eList = SortingPerm<int, int>(RetCand);
  MyVector<int> nCand(len);
  for (int i = 0; i < len; i++) {
    int j = eList[i];
    nCand(i) = RetCand[j];
  }
  return nCand;
}

// Canicalize the vector under symn action + the multiplication by
// non-zero element of Fd (field of the prime d)
MyVector<int> CyclicCanonicalization_SymN_fact(MyVector<int> const &V,
                                               int const &d) {
  MyVector<int> Vcan = CyclicCanonicalization_SymN(V, d);
  for (int mult = 2; mult <= d - 1; mult++) {
    MyVector<int> Vm = mult * V;
    MyVector<int> VmCan = CyclicCanonicalization_SymN(Vm, d);
    if (VmCan < Vcan)
      Vcan = VmCan;
  }
  return Vcan;
}


std::vector<std::vector<int>> SHORT_GetCandidateCyclic_Optimized(int const &n,
                                                                 int const &d) {
  int res = d % 2;
  int MaxVal = (d - res) / 2;
  auto Canonicalization =
      [&d](std::vector<int> const &fCand) -> std::vector<int> {
    std::vector<int> RetCand;
    for (auto &eVal : fCand) {
      int nVal;
      if (2 * eVal > d) {
        nVal = d - eVal;
      } else {
        nVal = eVal;
      }
      RetCand.push_back(nVal);
    }
    int len = fCand.size();
    std::vector<int> eList = SortingPerm<int, int>(RetCand);
    std::vector<int> nCand(len);
    for (int i = 0; i < len; i++) {
      int j = eList[i];
      nCand[i] = RetCand[j];
    }
    return nCand;
  };
  auto IsMinimal = [&](std::vector<int> const &eCand) -> bool {
    for (int mult = 2; mult <= d - 1; mult++) {
      std::vector<int> eProd;
      for (auto &x : eCand) {
        int y = x * mult;
        int res = y % d;
        eProd.push_back(res);
      }
      std::vector<int> NewCand = Canonicalization(eProd);
      if (NewCand < eCand)
        return false;
    }
    return true;
  };
  std::vector<std::vector<int>> ListCand = {{1}};
  for (int iDim = 1; iDim < n; iDim++) {
    std::vector<std::vector<int>> NewListCand;
    for (auto &eCand : ListCand) {
      int LastVal = eCand[iDim - 1];
      for (int i = LastVal; i <= MaxVal; i++) {
        std::vector<int> NewCand = ConcatenateVect(eCand, {i});
        if (IsMinimal(NewCand)) {
          NewListCand.push_back(NewCand);
        }
      }
    }
    ListCand = NewListCand;
#ifdef DEBUG_SHORTEST_CONFIG
    std::cerr << "SHORT: iDim=" << iDim << " |ListCand|=" << ListCand.size() << "\n";
#endif
  }
  return ListCand;
}

struct PrimeListAllowed {
  int p;
  bool DoWeListFeasible;
  std::vector<MyVector<int>> ListCases;
};

template <typename T>
bool IsMatchingListOfPrimes(std::vector<PrimeListAllowed> const &ListPrime,
                            MyMatrix<T> const &M) {
  std::vector<MyVector<int>> ListClasses = ComputeTranslationClasses<T, int>(M);
  int n = M.rows();
  MyMatrix<T> eInv = Inverse(M);
  auto GetOrder = [&](MyVector<int> const &eV) -> int {
    if (IsZeroVector(eV))
      return -1;
    MyVector<T> V2(n);
    for (int i = 0; i < n; i++) {
      T sum(0);
      for (int j = 0; j < n; j++) {
        sum += eV(j) * eInv(j, i);
      }
      V2(i) = sum;
    }
    int ord = 1;
    while (true) {
      bool IsProductCorr = true;
      for (int i = 0; i < n; i++) {
        T eProd = ord * V2(i);
        if (!IsInteger(eProd))
          IsProductCorr = false;
      }
      if (IsProductCorr)
        return ord;
      ord++;
    }
    return -1;
  };
  for (auto &eV : ListClasses) {
    int ord = GetOrder(eV);
    for (auto &eCaseP : ListPrime) {
      if (ord == eCaseP.p) {
        auto IsCorrectClass = [&](MyVector<int> const &W) -> bool {
          for (auto &W2 : eCaseP.ListCases) {
            if (W2 == W) {
              if (eCaseP.DoWeListFeasible)
                return true;
              else
                return false;
            }
          }
          if (eCaseP.DoWeListFeasible)
            return false;
          else
            return true;
        };
        MyVector<int> Vtest(n);
        for (int i = 0; i < n; i++) {
          T sum(0);
          for (int j = 0; j < n; j++) {
            sum += eV(j) * eInv(j, i);
          }
          sum *= ord;
#ifdef SANITY_CHECK_SHORTEST_CONFIG
          if (!IsInteger(sum)) {
            std::cerr << "SHORT: The sum should be integral\n";
            throw TerminalException{1};
          }
#endif
          int sum_i = UniversalScalarConversion<int, T>(sum);
          Vtest(i) = sum_i;
        }
        MyVector<int> VtestCan = CyclicCanonicalization_SymN_fact(Vtest, ord);
        if (!IsCorrectClass(VtestCan)) {
          return false;
        }
      }
    }
  }
  return true;
}

template <typename Tint>
Tint SHORT_GetMaximumDeterminant(MyMatrix<Tint> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  IteratorBinomial<mpz_class> eIter(nbRow, nbCol);
  std::vector<Face> ListFace = eIter.ListAllFace();
  bool IsFirst = true;
  Tint eRet = -2;
  for (auto &eFace : ListFace) {
    MyMatrix<Tint> Mred = SelectRow(M, eFace);
    Tint eDet = DeterminantMat(Mred);
    Tint eDetA = T_abs(eDet);
    if (IsFirst) {
      eRet = eDetA;
      IsFirst = false;
    } else {
      if (eDetA > eRet) {
        eRet = eDetA;
      }
    }
  }
  return eRet;
}

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
SHORT_SpannSimplicial(MyMatrix<Tint> const &M,
                      std::vector<MyMatrix<Tint>> const &ListSHVinp,
                      std::ostream &os) {
  Tint eMaxDet = SHORT_GetMaximumDeterminant(M);
  int n = M.cols();
  int nbVect = M.rows();
  std::vector<MyMatrix<Tint>> ListMatrGen =
      SHORT_GetStabilizer<T, Tint, Tgroup>(M, os);
  //
  // Building the set of inequalities
  //
  IteratorBinomial<mpz_class> eIter(nbVect, n - 1);
  std::vector<Face> ListAllFace = eIter.ListAllFace();
  std::vector<MyVector<T>> ListIneq;
  MyMatrix<T> M_T = UniversalMatrixConversion<T, Tint>(M);
  for (auto const &eFace : ListAllFace) {
    MyMatrix<T> Mred = SelectRow(M_T, eFace);
    if (RankMat(Mred) == n - 1) {
      MyVector<T> eVect(n);
      for (int i = 0; i < n; i++) {
        MyMatrix<T> eLine(n, 1);
        for (int j = 0; j < n; j++)
          eLine(j, 0) = 0;
        eLine(i, 0) = 1;
        MyMatrix<T> Mdet = Concatenate(Mred, eLine);
        T eDet = DeterminantMat(Mdet);
        eVect(i) = eDet;
      }
      MyVector<T> eIneq(n + 1), fIneq(n + 1);
      eIneq(0) = eMaxDet;
      fIneq(0) = eMaxDet;
      for (int i = 0; i < n; i++) {
        eIneq(i + 1) = eVect(i);
        fIneq(i + 1) = -eVect(i);
      }
      ListIneq.push_back(eIneq);
      ListIneq.push_back(fIneq);
    }
  }
  MyMatrix<T> FAC = MatrixFromVectorFamily(ListIneq);
  MyMatrix<T> EXT = cdd::DualDescription(FAC, os);
  std::vector<MyVector<Tint>> ListPt =
      GetListIntegralPoint<T, Tint>(FAC, EXT, os);
  //
  // Breaking into orbits
  //
  std::vector<MyMatrix<Tint>> ListGen;
  for (auto &eGen : ListMatrGen)
    ListGen.push_back(TransposedMat(eGen));
  std::function<MyVector<Tint>(MyVector<Tint> const &, MyMatrix<Tint> const &)>
      TheAct = [](MyVector<Tint> const &x,
                  MyMatrix<Tint> const &M) -> MyVector<Tint> { return M * x; };
  std::vector<MyVector<Tint>> ListRepr =
      OrbitSplittingGeneralized(ListPt, ListGen, TheAct);
  //
  // Building the set of inequalities
  //
  auto IsPresent = [&](MyMatrix<Tint> const &P) -> bool {
    for (auto &P2 : ListSHVinp) {
      std::optional<MyMatrix<Tint>> eResEquiv =
          SHORT_TestEquivalence<T, Tint, Tgroup>(P, P2, os);
      if (eResEquiv)
        return true;
    }
    return false;
  };
  auto PassFacetIsoCheck = [&](MyMatrix<Tint> const &U) -> bool {
    int len = U.rows();
    for (int i = 0; i < len; i++) {
      Face eFace(len);
      for (int j = 0; j < len; j++)
        eFace[j] = 1;
      eFace[i] = 0;
      MyMatrix<Tint> Usel = SelectRow(U, eFace);
      if (!IsPresent(Usel))
        return false;
    }
    return true;
  };
  std::vector<MyMatrix<Tint>> ListSpann;
  auto FuncInsert = [&](MyMatrix<Tint> const &Mnew) -> void {
    if (!PassFacetIsoCheck(Mnew))
      return;
    for (auto &P2 : ListSpann) {
      std::optional<MyMatrix<Tint>> eResEquiv =
          SHORT_TestEquivalence<T, Tint, Tgroup>(Mnew, P2, os);
      if (eResEquiv) {
        return;
      }
    }
    ReplyRealizability<T, Tint> eTestRes =
        SHORT_TestRealizabilityShortestFamily<T, Tint, Tgroup>(Mnew,
                                                               os);
    if (eTestRes.reply && eTestRes.replyCone)
      ListSpann.push_back(Mnew);
  };
  for (auto &ePt : ListRepr) {
    MyMatrix<Tint> eLine(n, 1);
    Tint eNorm = 0;
    for (int i = 0; i < n; i++) {
      Tint eVal = ePt(i);
      eNorm += eVal * eVal;
      eLine(i) = eVal;
    }
    if (eNorm > 0) {
      MyMatrix<Tint> eCand = Concatenate(M, eLine);
      FuncInsert(eCand);
    }
  }
  return ListSpann;
}

// clang-format off
#endif  // SRC_SHORT_SHORT_ENUMERATION_H_
// clang-format on
