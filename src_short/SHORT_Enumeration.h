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

// clang-format off
#endif  // SRC_SHORT_SHORT_ENUMERATION_H_
// clang-format on
