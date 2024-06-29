// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_CTYPE_MPI_CTYPEMPI_TYPES_H_
#define SRC_CTYPE_MPI_CTYPEMPI_TYPES_H_

// clang-format off
#include "Boost_bitset.h"
#include "COMB_Combinatorics.h"
#include "COMB_Combinatorics_buildset.h"
#include "MAT_Matrix.h"
#include "POLY_c_cddlib.h"
#include "POLY_cddlib.h"
#include "PolytopeEquiStabInt.h"
#include <cstring>
#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
// clang-format off

// #define DEBUG
//#define TIMINGS
//#define PRINT_FLIP
//#define PRINT_TRIPLE
//#define PRINT_GET_ADJ
// #define DEBUG_CRITERION_ELIMINATION

// The expression of the vectors is obtained from Section 5 of
// Dutour Sikiric, Grishukhin, Zonotopes and Parallelotopes,
// Southeast Asian Bulletin of Mathematics (2017) 41: 197–207
// That is, it is formed as { \pm e(S) for S \subset {1, ..., n} }.
template<typename T>
MyMatrix<T> GetPrincipalDomain(int const& n) {
  int n_elt = 1;
  for (int i=0; i<n; i++)
    n_elt *= 2;
  n_elt -= 1;
  MyMatrix<T> Mret(n_elt, n);
  int pos = 0;
  MyMatrix<int> M = BuildSet(n, 2);
  int n_row = M.rows();
  for (int i_row=0; i_row<n_row; i_row++) {
    int esum = 0;
    for (int i=0; i<n; i++)
      esum += M(i_row,i);
    if (esum > 0) {
      for (int i=0; i<n; i++)
        Mret(pos, i) = M(i_row, i);
      pos += 1;
    }
  }
  if (pos != n_elt) {
    std::cerr << "pos=" << pos << " n_elt=" << n_elt << "\n";
    throw TerminalException{1};
  }
  return Mret;
}

template <typename T> MyMatrix<T> ReduceExpandedMatrix(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int n = M.cols();
  int nbPair = nbRow / 2;
  MyMatrix<T> Mret(nbPair, n);
  auto is_sign_ok = [&](int const &iRow) -> bool {
    for (int i = 0; i < n; i++) {
      T eVal = M(iRow, i);
      if (eVal != 0) {
        return eVal > 0;
      }
    }
  };
  int iPair = 0;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    if (is_sign_ok(iRow)) {
      for (int i = 0; i < n; i++)
        Mret(iPair, i) = M(iRow, i);
      iPair++;
    }
  }
  return Mret;
}

template <typename T> struct TypeCtypeExch {
  MyMatrix<T> eMat;
};

template <typename T>
void PairExch_to_vectorchar(TypeCtypeExch<T> const &eCtype, int const &nbRow,
                            int const &nbCol, char *ptr_o, [[maybe_unused]] std::ostream& os) {
#ifdef ERR_LOG
  os << "PairExch_to_vectorchar, Begin\n";
#endif
  for (int i = 0; i < nbRow; i++)
    for (int j = 0; j < nbCol; j++) {
      std::memcpy(ptr_o, (char *)(&eCtype.eMat(i, j)), sizeof(T));
      ptr_o += sizeof(T);
    }
#ifdef ERR_LOG
  os << "PairExch_to_vectorchar, End\n";
#endif
}

template <typename T>
TypeCtypeExch<T> vectorchar_to_PairExch(char *ptr_i, int const &nbRow,
                                        int const &nbCol,
                                        [[maybe_unused]] std::ostream& os) {
#ifdef ERR_LOG
  os << "vectorchar_to_PairExch, Begin\n";
#endif
  MyMatrix<T> eMat(nbRow, nbCol);
  for (int i = 0; i < nbRow; i++)
    for (int j = 0; j < nbCol; j++) {
      T eVal;
      std::memcpy((char *)(&eVal), ptr_i, sizeof(T));
      ptr_i += sizeof(T);
      eMat(i, j) = eVal;
    }
#ifdef ERR_LOG
  os << "vectorchar_to_PairExch, End\n";
#endif
  return {eMat};
}

struct TypeAdjExch {
  int8_t iProc1;
  int pos1;
  int8_t iProc2;
  int pos2;
};

std::ostream &operator<<(std::ostream &os, TypeAdjExch const &obj) {
  os << "(" << int(obj.iProc1) << "-" << obj.pos1 << "|" << int(obj.iProc2)
     << "-" << obj.pos2 << ")";
  return os;
}

void TypeAdjExch_to_ptrchar(TypeAdjExch const &eR, char *ptr_o) {
  std::memcpy(ptr_o, (char *)(&eR.iProc1), sizeof(int8_t));
  ptr_o += sizeof(int8_t);
  std::memcpy(ptr_o, (char *)(&eR.pos1), sizeof(int));
  ptr_o += sizeof(int);
  //
  std::memcpy(ptr_o, (char *)(&eR.iProc2), sizeof(int8_t));
  ptr_o += sizeof(int8_t);
  std::memcpy(ptr_o, (char *)(&eR.pos2), sizeof(int));
  ptr_o += sizeof(int);
}

TypeAdjExch ptrchar_to_TypeAdjExch(char *ptr_i) {
  int8_t iProc1;
  std::memcpy((char *)(&iProc1), ptr_i, sizeof(int8_t));
  ptr_i += sizeof(int8_t);
  //
  int pos1;
  std::memcpy((char *)(&pos1), ptr_i, sizeof(int));
  ptr_i += sizeof(int);
  //
  int8_t iProc2;
  std::memcpy((char *)(&iProc2), ptr_i, sizeof(int8_t));
  ptr_i += sizeof(int8_t);
  //
  int pos2;
  std::memcpy((char *)(&pos2), ptr_i, sizeof(int));
  ptr_i += sizeof(int);
  //
  return {iProc1, pos1, iProc2, pos2};
}

std::vector<char> TypeAdjExch_to_stdvectorchar(TypeAdjExch const &eR) {
  int len = sizeof(int8_t) + 2 * (sizeof(int) + sizeof(int8_t));
  std::vector<char> eV(len);
  char *ptr_o = eV.data();
  //
  int8_t val = 1;
  std::memcpy(ptr_o, (char *)(&val), sizeof(int8_t));
  ptr_o += sizeof(int8_t);
  //
  TypeAdjExch_to_ptrchar(eR, ptr_o);
  //
  return eV;
}

template <typename T> struct TypeCtypeAdjExch {
  MyMatrix<T> eMat;
  int8_t iProc;
  int pos;
};

template <typename T>
void PairAdjExch_to_ptrchar(TypeCtypeAdjExch<T> const &eCtype, int const &nbRow,
                            int const &nbCol, char *ptr_o) {
#ifdef ERR_LOG
  std::cerr << "PairAdjExch_to_vectorchar, Begin\n";
#endif
  for (int i = 0; i < nbRow; i++)
    for (int j = 0; j < nbCol; j++) {
      std::memcpy(ptr_o, (char *)(&eCtype.eMat(i, j)), sizeof(T));
      ptr_o += sizeof(T);
    }
  std::memcpy(ptr_o, (char *)(&eCtype.iProc), sizeof(int8_t));
  ptr_o += sizeof(int8_t);
  std::memcpy(ptr_o, (char *)(&eCtype.pos), sizeof(int));
  ptr_o += sizeof(int);
#ifdef ERR_LOG
  std::cerr << "PairAdjExch_to_vectorchar, End\n";
#endif
}

template <typename T>
TypeCtypeAdjExch<T> ptrchar_to_PairAdjExch(char *ptr_i, int const &nbRow,
                                           int const &nbCol) {
#ifdef ERR_LOG
  std::cerr << "vectorchar_to_PairAdjExch, Begin\n";
#endif
  MyMatrix<T> eMat(nbRow, nbCol);
  for (int i = 0; i < nbRow; i++)
    for (int j = 0; j < nbCol; j++) {
      T eVal;
      std::memcpy((char *)(&eVal), ptr_i, sizeof(T));
      ptr_i += sizeof(T);
      eMat(i, j) = eVal;
    }
  int8_t iProc;
  std::memcpy((char *)(&iProc), ptr_i, sizeof(int8_t));
  ptr_i += sizeof(int8_t);
  //
  int pos;
  std::memcpy((char *)(&pos), ptr_i, sizeof(int));
  ptr_i += sizeof(int);
#ifdef ERR_LOG
  std::cerr << "vectorchar_to_PairAdjExch, End\n";
#endif
  return {eMat, iProc, pos};
}

template <typename T>
std::vector<char>
PairAdjExch_to_stdvectorchar(TypeCtypeAdjExch<T> const &eCtype,
                             int const &nbRow, int const &nbCol) {
  int len =
      sizeof(int8_t) + nbRow * nbCol * sizeof(T) + sizeof(int8_t) + sizeof(int);
  std::vector<char> eV(len);
  char *ptr_o = eV.data();
  //
  int8_t val = 0;
  std::memcpy(ptr_o, (char *)(&val), sizeof(int8_t));
  ptr_o += sizeof(int8_t);
  //
  PairAdjExch_to_ptrchar(eCtype, nbRow, nbCol, ptr_o);
  //
  return eV;
}

template <typename T> std::vector<MyMatrix<T>> CTYP_GetBasis(int n) {
  std::vector<MyMatrix<T>> ListSymmMat;
  for (int i = 0; i < n; i++) {
    MyMatrix<T> eMat = ZeroMatrix<T>(n, n);
    eMat(i, i) = 1;
    ListSymmMat.push_back(eMat);
  }
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++) {
      MyMatrix<T> eMat = ZeroMatrix<T>(n, n);
      eMat(i, j) = 1;
      eMat(j, i) = 1;
      ListSymmMat.push_back(eMat);
    }
  return ListSymmMat;
}

template<typename Tidx>
struct triple {
  Tidx i;
  Tidx j;
  Tidx k;
};

namespace std {
template<typename Tidx> struct hash<triple<Tidx>> {
  std::size_t operator()(const triple<Tidx> &et) const {
    size_t seed = 34;
    auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
      seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };
    combine_hash(seed, std::hash<Tidx>()(et.i));
    combine_hash(seed, std::hash<Tidx>()(et.j));
    combine_hash(seed, std::hash<Tidx>()(et.k));
    return seed;
  }
};
// clang-format off
}  // namespace std
// clang-format on

template <typename Tidx>
bool operator==(triple<Tidx> const &obj1, triple<Tidx> const &obj2) {
  return obj1.i == obj2.i && obj1.j == obj2.j && obj1.k == obj2.j;
}

template <typename T> bool CheckCoveringParityClasses(MyMatrix<T> const &M) {
  size_t n_rows = M.rows();
  size_t n_cols = M.cols();
  std::vector<int> ListStatus(n_rows, 0);
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    int pos = -1;
    int e_pow = 1;
    T eTwo = 2;
    for (size_t i = 0; i < n_cols; i++) {
      T res_T = ResInt(M(i_row, i), eTwo);
      int res = UniversalScalarConversion<int, T>(res_T);
      pos += res * e_pow;
      e_pow *= 2;
    }
    if (pos == -1)
      return false;
    ListStatus[pos] += 1;
  }
  for (size_t i_row = 0; i_row < n_rows; i_row++)
    if (ListStatus[i_row] != 1)
      return false;
  return true;
}

template <typename T, typename Tidx>
MyMatrix<T> CTYP_TheFlipping(MyMatrix<T> const &TheCtype,
                             std::vector<triple<Tidx>> const &TheInfo) {
  size_t n_rows = TheCtype.rows();
  size_t n_cols = TheCtype.cols();
  size_t n_rows_ret = n_rows / 2;
#ifdef PRINT_FLIP
  std::cerr << "CTYP_TheFlipping n_rows=" << n_rows << " n_cols=" << n_cols
            << "\n";
#endif
  Face ListIchange(n_rows);
  for (auto &e_triple : TheInfo) {
#ifdef PRINT_FLIP
    std::cerr << "e_triple=" << static_cast<int>(e_triple.i) << " , "
              << static_cast<int>(e_triple.j) << " , "
              << static_cast<int>(e_triple.k) << "\n";
#endif
    ListIchange[e_triple.i] = 1;
  }
#ifdef PRINT_FLIP
  std::cerr << "ListIchange =";
  for (size_t i_row = 0; i_row < n_rows; i_row++)
    std::cerr << " " << ListIchange[i_row];
  std::cerr << "\n";
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    if (ListIchange[i_row] == 1) {
      std::cerr << "Removed line =";
      for (size_t i = 0; i < n_cols; i++)
        std::cerr << " " << TheCtype(i_row, i);
      std::cerr << "\n";
    }
  }
#endif
  MyMatrix<T> RetMat(n_rows_ret, n_cols);
  std::vector<T> V(n_cols);
  size_t idx = 0;
  Face IsAssigned(n_rows_ret);
  auto insert_if_signok = [&]() -> void {
    for (size_t i = 0; i < n_cols; i++) {
      T eVal = V[i];
      if (eVal != 0) {
        if (eVal > 0) {
          int pos = -1;
          int e_pow = 1;
          T eTwo = 2;
          for (size_t i = 0; i < n_cols; i++) {
            T res_T = ResInt(V[i], eTwo);
            int res = UniversalScalarConversion<int, T>(res_T);
            pos += res * e_pow;
            e_pow *= 2;
          }
          if (IsAssigned[pos] == 0) {
            for (size_t i_col = 0; i_col < n_cols; i_col++)
              RetMat(idx, i_col) = V[i_col];
            idx++;
            IsAssigned[pos] = 1;
          }
        }
        return;
      }
    }
  };
  for (auto &e_triple : TheInfo) {
    Tidx j = e_triple.j;
    Tidx k = e_triple.k;
    //
    for (size_t i_col = 0; i_col < n_cols; i_col++)
      V[i_col] = -TheCtype(j, i_col) + TheCtype(k, i_col);
    insert_if_signok();
    //
    for (size_t i_col = 0; i_col < n_cols; i_col++)
      V[i_col] = TheCtype(j, i_col) - TheCtype(k, i_col);
    insert_if_signok();
  }
  for (size_t i_row = 0; i_row < n_rows; i_row++) {
    if (ListIchange[i_row] == 0 && i_row % 2 == 0) {
      for (size_t i_col = 0; i_col < n_cols; i_col++)
        RetMat(idx, i_col) = TheCtype(i_row, i_col);
      idx++;
    }
  }
#ifdef DEBUG
  if (n_rows_ret != idx) {
    std::cerr << "Incoherence in idx, n_rows_ret. idx=" << idx
              << " n_rows_ret=" << n_rows_ret << "\n";
    throw TerminalException{1};
  }
  if (!CheckCoveringParityClasses(RetMat)) {
    std::cerr << "RetMat parity classes is not correct\n";
    std::cerr << "RetMat=\n";
    WriteMatrix(std::cerr, RetMat);
    throw TerminalException{1};
  }
#endif
  return RetMat;
}

template <typename T, typename Tidx>
std::pair<std::vector<triple<Tidx>>, std::vector<Tidx>>
CTYP_GetListTriple(MyMatrix<T> const &TheCtype, [[maybe_unused]] std::ostream& os) {
  int n_edge = TheCtype.rows();
  int n_edgered = n_edge / 2;
  int n_cols = TheCtype.cols();
#ifdef PRINT_TRIPLE
  os << "n_edge=" << n_edge << " n_cols=" << n_cols << "\n";
  os << "TEST TheCtype=\n";
  WriteMatrix(os, TheCtype);
#endif
  std::vector<triple<Tidx>> ListTriples;
  std::vector<Tidx> MappingVect(n_edgered * n_edgered, -1);
  T eTwo = 2;
  auto get_position = [&](MyVector<T> const &eV, Tidx start_idx) -> Tidx {
    auto get_nature = [&](Tidx pos) -> bool {
#ifdef PRINT_TRIPLE
      os << "get_nature pos=" << pos << "\n";
#endif
      for (Tidx i_col = 0; i_col < n_cols; i_col++)
        if (TheCtype(pos, i_col) != eV(i_col))
          return false;
      return true;
    };
    auto get_value = [&]() -> Tidx {
      int pos = -1;
      int e_pow = 1;
      for (int i = 0; i < n_cols; i++) {
        T res_T = ResInt(eV(i), eTwo);
        int res = UniversalScalarConversion<int, T>(res_T);
        pos += res * e_pow;
        e_pow *= 2;
      }
#ifdef PRINT_TRIPLE
      os << "eV =";
      for (int i = 0; i < n_cols; i++)
        os << " " << eV(i);
      os << "\n";
      os << "Found pos=" << pos << "\n";
#endif
      if (pos == -1)
        return -1;
      if (get_nature(2 * pos))
        return 2 * pos;
      if (get_nature(2 * pos + 1))
        return 2 * pos + 1;
      return -1;
    };
    int pos = get_value();
#ifdef PRINT_TRIPLE
    os << "After get_value pos=" << pos << "\n";
#endif
    if (pos > start_idx)
      return pos;
    return -1;
  };
  MyVector<T> eDiff(n_cols);
  for (Tidx i = 0; i < n_edge; i++)
    for (Tidx j = i + 1; j < n_edge; j++) {
#ifdef PRINT_TRIPLE
      os << "i=" << static_cast<int>(i) << " j=" << static_cast<int>(j)
         << "\n";
#endif
      for (Tidx i_col = 0; i_col < n_cols; i_col++)
        eDiff(i_col) = -TheCtype(i, i_col) - TheCtype(j, i_col);
#ifdef PRINT_TRIPLE
      os << "We have eDiff=" << StringVectorGAP(eDiff) << "\n";
#endif
      Tidx k = get_position(eDiff, j);
      Tidx crit = -1;
#ifdef PRINT_TRIPLE
      os << "k=" << static_cast<int>(k) << "\n";
#endif
      if (k != crit) {
        ListTriples.push_back({i, j, k});
        ListTriples.push_back({j, k, i});
        ListTriples.push_back({k, i, j});
        //
        Tidx ired = i / 2;
        Tidx jred = j / 2;
        Tidx kred = k / 2;
#ifdef PRINT_TRIPLE
        os << "n_edgered=" << n_edgered << " i/j/kred=" << ired << " "
           << jred << " " << kred << "\n";
#endif
        MappingVect[ired * n_edgered + jred] = kred;
        MappingVect[jred * n_edgered + ired] = kred;
        //
        MappingVect[ired * n_edgered + kred] = jred;
        MappingVect[kred * n_edgered + ired] = jred;
        //
        MappingVect[jred * n_edgered + kred] = ired;
        MappingVect[kred * n_edgered + jred] = ired;
      }
    }
#ifdef PRINT_TRIPLE
  os << "Exiting CTYP_GetListTriple\n";
#endif
  return {std::move(ListTriples), std::move(MappingVect)};
}

template <typename T> MyMatrix<T> ExpressMatrixForCType(MyMatrix<T> const &M, [[maybe_unused]] std::ostream& os) {
  int n = M.cols();
  int nbRow = M.rows();
#ifdef PRINT_EXPRESS
  os << "n=" << n << " nbRow=" << nbRow << "\n";
#endif
  MyMatrix<T> Mret(2 * nbRow, n);
#ifdef DEBUG
  std::vector<int> ListStatus(nbRow, 0);
#endif
  T eTwo = 2;
  for (int iRow = 0; iRow < nbRow; iRow++) {
#ifdef PRINT_EXPRESS
    os << "iRow=" << iRow << "/" << nbRow << "\n";
#endif
    int pos = -1;
    int e_pow = 1;
    for (int i = 0; i < n; i++) {
      T res_T = ResInt(M(iRow, i), eTwo);
      int res = UniversalScalarConversion<int, T>(res_T);
#ifdef PRINT_EXPRESS
      os << "  i=" << i << " M(iRow,i)=" << M(iRow, i)
         << " res_T=" << res_T << " res=" << res << " e_pow=" << e_pow
         << "\n";
#endif
      pos += res * e_pow;
      e_pow *= 2;
    }
#ifdef PRINT_EXPRESS
    os << "  pos=" << pos << "\n";
#endif
#ifdef DEBUG
    ListStatus[pos] += 1;
#endif
    for (int i = 0; i < n; i++) {
      Mret(2 * pos, i) = M(iRow, i);
      Mret(2 * pos + 1, i) = -M(iRow, i);
    }
  }
#ifdef DEBUG
  for (int i = 0; i < nbRow; i++)
    if (ListStatus[i] != 1) {
      std::cerr << "Consistency error at i=" << i << "\n";
      throw TerminalException{1};
    }
#endif
  return Mret;
}

template <typename T, typename Tidx> struct DataCtypeFacet {
  MyMatrix<T> TheCtype;
  MyMatrix<T> ListInequalities;
  std::vector<std::vector<triple<Tidx>>> ListInformations;
  std::vector<int> ListIrred;
  int nb_triple;
  int nb_ineq;
  int nb_ineq_after_crit;
};

template <typename T, typename Tidx>
DataCtypeFacet<T, Tidx>
CTYP_GetConeInformation(TypeCtypeExch<T> const &TheCtypeArr, std::ostream& os) {
#ifdef TIMINGS
  MicrosecondTime time;
#endif
#ifdef PRINT_GET_ADJ
  os << "CTYP_GetConeInformation, step 1\n";
#endif
  MyMatrix<T> TheCtype = ExpressMatrixForCType(TheCtypeArr.eMat, os);
  int n_edge = TheCtype.rows();

#ifdef TIMINGS
  os << "|ExpressMatrixForCType|=" << time << "\n";
#endif
#ifdef PRINT_GET_ADJ
  os << "CTYP_GetConeInformation, step 2\n";
#endif
  std::pair<std::vector<triple<Tidx>>, std::vector<Tidx>> PairTriple =
    CTYP_GetListTriple<T, Tidx>(TheCtype, os);

#ifdef TIMINGS
  os << "|CTYP_GetListTriple|=" << time << "\n";
#endif
#ifdef PRINT_GET_ADJ
  os << "CTYP_GetConeInformation, step 3\n";
#endif
  Tidx n = TheCtype.cols();
  uint8_t n_i = static_cast<uint8_t>(n);
  Tidx tot_dim = n * (n + 1) / 2;
  // The inequality is written as
  // Q[v_i] <= Q[2 v_k + v_i]  with  v_i + v_j + v_k = 0
  // Q[v_i] <= Q[v_k - v_j]  with  v_i + v_j + v_k = 0
  // We have Q[v + w] + Q[v - w] = 2 Q[v] + 2 Q[w]
  // So, Q[v_i] + Q[v_k - v_j] = 2 Q[v_k] + 2 Q[v_j]
  // So, the inequality is equivalent to
  // Q[v_i] <= 2 Q[v_k] + 2 Q[v_j] - Q[v_i]
  // or Q[v_i] <= Q[v_k] + Q[v_j]
  // So, yes, we have the right inequality.
  auto ComputeInequality = [&](MyVector<T> const &V1,
                               MyVector<T> const &V2) -> MyVector<T> {
    MyVector<T> TheVector(tot_dim);
    int idx = 0;
    for (uint8_t i = 0; i < n_i; i++) {
      TheVector(idx) = V1(i) * V1(i) - V2(i) * V2(i);
      idx++;
    }
    for (uint8_t i = 0; i < n_i; i++)
      for (uint8_t j = i + 1; j < n_i; j++) {
        // Factor 2 removed for simplification and faster code.
        TheVector(idx) = V1(i) * V1(j) - V2(i) * V2(j);
        idx++;
      }
    return TheVector;
  };
  std::unordered_map<MyVector<T>, std::vector<triple<Tidx>>> Tot_map;
  MyVector<T> V1(n), V2(n);
  auto FuncInsertInequality = [&](Tidx i, Tidx j, Tidx k) -> void {
    for (Tidx i_col = 0; i_col < n; i_col++) {
      V1(i_col) = 2 * TheCtype(k, i_col) + TheCtype(i, i_col);
      V2(i_col) = TheCtype(i, i_col);
    }
    MyVector<T> TheVector = ComputeInequality(V1, V2);
    triple<Tidx> TheInfo = {i, j, k};
    std::vector<triple<Tidx>> &list_trip = Tot_map[TheVector];
    list_trip.push_back(TheInfo);
  };
  for (auto &e_triple : PairTriple.first)
    FuncInsertInequality(e_triple.i, e_triple.j, e_triple.k);
  int nb_ineq = Tot_map.size();
#ifdef PRINT_GET_ADJ
  os << "Input |Tot_map|=" << Tot_map.size() << "\n";
#endif

#ifdef DEBUG_CRITERION_ELIMINATION
  os << "|Tot_map|=" << Tot_map.size() << "\n";
  std::vector<int> ListNbMatch(n_edge, 0);
  for (auto &et : PairTriple.first) {
    ListNbMatch[et.i]++;
    ListNbMatch[et.j]++;
    ListNbMatch[et.k]++;
  }
  os << "ListNbMatch =";
  for (auto &eV : ListNbMatch)
    os << " " << eV;
  os << "\n";
  int j_ineq = 0;
  for (auto &kv : Tot_map) {
    os << "j_ineq=" << j_ineq << " ineq =";
    int e_dim = kv.first.size();
    for (int i = 0; i < e_dim; i++)
      os << " " << kv.first(i);
    os << " LSet =";
    for (auto &et : kv.second)
      os << " {" << static_cast<int>(et.i) << ","
         << static_cast<int>(et.j) << "," << static_cast<int>(et.k)
         << "}";
    os << "\n";
    j_ineq++;
  }
  os << "n_edge=" << n_edge << "\n";
  os << "TheCtype=\n";
  for (int i_edge = 0; i_edge < n_edge; i_edge++) {
    int n = TheCtype.cols();
    for (int i = 0; i < n; i++)
      os << " " << TheCtype(i_edge, i);
    os << "\n";
  }
  for (int i_edge = 0; i_edge < n_edge; i_edge++) {
    for (int j_edge = 0; j_edge < n_edge; j_edge++)
      os << " "
         << static_cast<int>(PairTriple.second[i_edge * n_edge + j_edge]);
    os << "\n";
  }
  int nb_triple_div3 = PairTriple.first.size() / 3;
  std::unordered_map<triple<Tidx>, int> MapTriple;
  for (int i_triple = 0; i_triple < nb_triple_div3; i_triple++) {
    MapTriple[PairTriple.first[3 * i_triple]] = i_triple;
  }
  os << "nb_triple_div3=" << nb_triple_div3 << "\n";
  for (int i_triple = 0; i_triple < nb_triple_div3; i_triple++) {
    triple et = PairTriple.first[3 * i_triple];
    os << "et=" << static_cast<int>(et.i) << " "
       << static_cast<int>(et.j) << " " << static_cast<int>(et.k)
       << "\n";
  }
#endif

#ifdef TIMINGS
  os << "|Insert inequalities|=" << time << "\n";
#endif
  int n_edgered = n_edge / 2;
#ifdef PRINT_GET_ADJ
  int nb_match = 0;
  int nb_pass = 0;
#endif
#ifdef PRINT_GET_ADJ
  for (int i_edge = 0; i_edge < n_edgered; i_edge++) {
    for (int j_edge = 0; j_edge < n_edgered; j_edge++)
      os << " "
         << static_cast<int>(PairTriple.second[i_edge * n_edgered + j_edge]);
    os << "\n";
  }
#endif
  // #define PRINT_GET_ADJ_O
  //  We apply here the 3 dimensional criterion for feasibility of C-type
  //  switches
  auto TestApplicabilityCriterion_with_e = [&](triple<Tidx> const &e_triple,
                                               Tidx const &e) -> bool {
#ifdef PRINT_GET_ADJ
    nb_pass++;
#endif
    Tidx i = e_triple.i / 2;
    Tidx j = e_triple.j / 2;
    Tidx k = e_triple.k / 2;
    Tidx crit = -1;
#ifdef PRINT_GET_ADJ_O
    os << "i=" << static_cast<int>(i) << " j=" << static_cast<int>(j)
       << " k=" << static_cast<int>(k) << " e=" << static_cast<int>(e)
       << "\n";
#endif
    //
    // testing e
    if (e == i || e == j || e == k)
      return false;
    //
    // getting f and testing it
    Tidx f = PairTriple.second[i * n_edgered + e];
#ifdef PRINT_GET_ADJ_O
    os << "f=" << static_cast<int>(f) << "\n";
#endif
    if (f == crit || f == j || f == k)
      return false;
    //
    // getting g and testing it
    Tidx g = PairTriple.second[j * n_edgered + e];
#ifdef PRINT_GET_ADJ_O
    os << "g=" << static_cast<int>(g) << "\n";
#endif
    if (g == crit || g == f || g == i || g == k)
      return false;
    //
    // getting h and testing it
    Tidx h = PairTriple.second[i * n_edgered + g];
#ifdef PRINT_GET_ADJ_O
    os << "h=" << static_cast<int>(h) << "\n";
#endif
    if (h == crit || h == f || h == e || h == j || h == k)
      return false;
    //
    // testing presence of {j,f,h}
    Tidx h2 = PairTriple.second[j * n_edgered + f];
#ifdef PRINT_GET_ADJ_O
    os << "h2=" << static_cast<int>(h2) << "\n";
#endif
    if (h2 != h)
      return false;
      //
      // We have the 7-uple.
#ifdef PRINT_GET_ADJ
    nb_match++;
#endif
    return true;
  };
  auto TestApplicabilityCriterion = [&](triple<Tidx> const &e_triple) -> bool {
    for (Tidx e = 0; e < n_edgered; e++)
      if (TestApplicabilityCriterion_with_e(e_triple, e))
        return true;
    return false;
  };
  std::vector<Tidx> ListResultCriterion(n_edge * n_edge, 0);
  for (auto &e_triple : PairTriple.first) {
    if (TestApplicabilityCriterion(e_triple)) {
      Tidx i = e_triple.i;
      Tidx j = e_triple.j;
      Tidx k = e_triple.k;
#ifdef PRINT_GET_ADJ
      os << "FOUND i=" << static_cast<int>(i)
         << " j=" << static_cast<int>(j) << " k=" << static_cast<int>(k)
         << "\n";
      os << "ENT1 = " << static_cast<int>(j) << " "
         << static_cast<int>(k) << "\n";
      os << "ENT2 = " << static_cast<int>(i) << " "
         << static_cast<int>(j) << "\n";
#endif
      ListResultCriterion[j * n_edge + k] = 1;
      ListResultCriterion[i * n_edge + j] = 1;
    }
  }
  auto TestFeasibilityListTriple =
      [&](std::vector<triple<Tidx>> const &list_triple) -> bool {
    for (auto &e_triple : list_triple) {
#ifdef PRINT_GET_ADJ
      os << "e_triple i=" << static_cast<int>(e_triple.i) << " "
         << static_cast<int>(e_triple.j) << "\n";
#endif
      if (ListResultCriterion[e_triple.i * n_edge + e_triple.j] == 1)
        return false;
    }
    return true;
  };
  // erasing the inequalities that are sure to be redundant.
#ifdef PRINT_GET_ADJ
  int nb_redund = 0;
#endif
  std::unordered_map<MyVector<T>, std::vector<triple<Tidx>>> Tot_mapB;
  for (auto &kv : Tot_map) {
    if (!TestFeasibilityListTriple(kv.second)) {
#ifdef PRINT_GET_ADJ
      nb_redund++;
#endif
    } else {
      Tot_mapB[kv.first] = kv.second;
    }
  }
#ifdef PRINT_GET_ADJ
  os << "nb_match=" << nb_match << " nb_pass=" << nb_pass << "\n";
  os << "nb_redund = " << nb_redund << "\n";
  os << "After criterion |Tot_mapB|=" << Tot_mapB.size() << "\n";
#endif

#ifdef TIMINGS
  os << "|Criterion Ineq Drop|=" << time << "\n";
#endif
#ifdef PRINT_GET_ADJ
  os << "CTYP_GetConeInformation, step 4\n";
#endif
  int nb_ineq_after_crit = Tot_mapB.size();
  MyMatrix<T> ListInequalities(nb_ineq_after_crit, tot_dim);
  std::vector<std::vector<triple<Tidx>>> ListInformations;
  size_t i_ineq = 0;
  for (auto &kv : Tot_mapB) {
    for (Tidx i_col = 0; i_col < tot_dim; i_col++)
      ListInequalities(i_ineq, i_col) = kv.first(i_col);
    i_ineq++;
    ListInformations.push_back(std::move(kv.second));
  }
#ifdef PRINT_GET_ADJ
  os << "CTYP_GetConeInformation, step 5\n";
  os << "ListInequalities=\n";
  WriteMatrix(os, ListInequalities);
#endif

#ifdef TIMINGS
  os << "|ListInformations|=" << time << "\n";
#endif
  std::vector<int> ListIrred =
      cbased_cdd::RedundancyReductionClarkson(ListInequalities);
#ifdef TIMINGS
  os << "|RedundancyReductionClarkson|=" << time << "\n";
#endif

#ifdef PRINT_GET_ADJ
  os << "|ListIrred|=" << ListIrred.size() << "\n";
  os << "ListIrred =";
  for (auto &idx : ListIrred)
    os << " " << idx;
  os << "\n";
  MyMatrix<T> ListInequalitiesIrred = SelectRow(ListInequalities, ListIrred);
  os << "ListInequalitiesIrred=\n";
  WriteMatrix(os, ListInequalitiesIrred);
#endif
  int nb_triple = PairTriple.first.size();
  return {std::move(TheCtype),
          std::move(ListInequalities),
          std::move(ListInformations),
          std::move(ListIrred),
          nb_triple,
          nb_ineq,
          nb_ineq_after_crit};
}

template <typename T, typename Tidx>
MyMatrix<T> CTYP_GetInsideGramMat(DataCtypeFacet<T, Tidx> const &data) {
  int n = data.TheCtype.cols();
  MyVector<T> V = GetSpaceInteriorPoint_Basic(data.ListInequalities);
  MyMatrix<T> eGram(n, n);
  int pos = 0;
  for (int i = 0; i < n; i++) {
    eGram(i, i) = 2 * V(pos);
    pos++;
  }
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++) {
      eGram(i, j) = V(pos);
      pos++;
    }
  return eGram;
}

template <typename T>
TypeCtypeExch<T> CTYP_GetCTypeFromGramMat(MyMatrix<T> const &eG,
                                          std::ostream &os) {
  int n = eG.rows();
  int n_elt = 1;
  for (int i = 0; i < n; i++)
    n_elt *= 2;
  n_elt -= 1;
  MyMatrix<T> Mret(n_elt, n);
  int pos = 0;
  MyMatrix<int> M = BuildSet(n, 2);
  int n_row = M.rows();
  MyVector<T> Vmid(n);
  for (int i_row = 0; i_row < n_row; i_row++) {
    int esum = 0;
    for (int i = 0; i < n; i++)
      esum += M(i_row, i);
    if (esum > 0) {
      for (int i = 0; i < n; i++)
        Vmid(i, 0) = M(i_row, i) / 2;
      resultCVP<T, T> result = CVPVallentinProgram_exact(eG, Vmid, os);
      int n_closest = result.ListVect.rows();
      if (n_closest != 2) {
        std::cerr << "n_closest=" << n_closest << "\n";
        std::cerr << "The number of closest vector is not equal to 2\n";
        throw TerminalException{1};
      }
      for (int i = 0; i < n; i++)
        Mret(pos, i) = result.ListVect(1, i) - result.ListVect(0, i);
      pos += 1;
    }
  }
  if (pos != n_elt) {
    std::cerr << "pos=" << pos << " n_elt=" << n_elt << "\n";
    throw TerminalException{1};
  }
  return Mret;
}

template <typename T, typename Tidx>
std::vector<TypeCtypeExch<T>>
CTYP_Kernel_GetAdjacentCanonicCtypes(TypeCtypeExch<T> const &TheCtypeArr,
                                     bool const &canonicalize,
                                     std::ostream &os) {
#ifdef TIMINGS
  MicrosecondTime time;
#endif
  DataCtypeFacet<T, Tidx> data = CTYP_GetConeInformation<T, Tidx>(TheCtypeArr, os);
#ifdef TIMINGS
  os << "|data|=" << time << "\n";
#endif
  std::vector<TypeCtypeExch<T>> ListCtype;
  for (auto &e_int : data.ListIrred) {
    MyMatrix<T> FlipMat =
        CTYP_TheFlipping(data.TheCtype, data.ListInformations[e_int]);
    if (canonicalize) {
      MyMatrix<T> CanMat =
          LinPolytopeAntipodalIntegral_CanonicForm(FlipMat, os);
      ListCtype.push_back({std::move(CanMat)});
    } else {
      ListCtype.push_back({std::move(FlipMat)});
    }
  }
#ifdef PRINT_GET_ADJ
  os << "CTYP_GetConeInformation, step 7\n";
#endif

#ifdef TIMINGS
  os << "|Flip + Canonic|=" << time << "\n";
#endif
  return ListCtype;
}

template <typename T>
std::vector<TypeCtypeExch<T>>
CTYP_GetAdjacentCanonicCtypes(TypeCtypeExch<T> const &TheCtypeArr,
                              std::ostream &os) {
  bool canonicalize = true;
  using Tidx = int8_t;
  return CTYP_Kernel_GetAdjacentCanonicCtypes<T, Tidx>(TheCtypeArr,
                                                       canonicalize, os);
}

struct StructuralInfo {
  int nb_triple;
  int nb_ineq;
  int nb_ineq_after_crit;
  int nb_free;
  int nb_autom;
};

namespace boost::serialization {
template <class Archive>
inline void serialize(Archive &ar, StructuralInfo &eRec,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("nb_triple", eRec.nb_triple);
  ar &make_nvp("nb_ineq", eRec.nb_ineq);
  ar &make_nvp("nb_ineq_after_crit", eRec.nb_ineq_after_crit);
  ar &make_nvp("nb_free", eRec.nb_free);
  ar &make_nvp("nb_autom", eRec.nb_autom);
}
} // namespace boost::serialization

template <typename T>
int CTYP_GetNumberFreeVectors(TypeCtypeExch<T> const &TheCtypeArr) {
  int n = TheCtypeArr.eMat.cols();
  BlockIteration blk(n, 2);
  int n_vect = TheCtypeArr.eMat.rows();
  blk.IncrementSilent();
  int nb_free = 0;
  T eTwo = 2;
  while (true) {
    std::vector<int> eVect = blk.GetVect();
    std::vector<MyVector<T>> ListVect;
    for (int i_vect = 0; i_vect < n_vect; i_vect++) {
      MyVector<T> eV = GetMatrixRow(TheCtypeArr.eMat, i_vect);
      T eSum = 0;
      for (int i = 0; i < n; i++)
        eSum += eV(i) * eVect[i];
      T res = ResInt(eSum, eTwo);
      if (res == 0) {
        ListVect.push_back(eV);
      }
    }
    MyMatrix<T> eMat = MatrixFromVectorFamily(ListVect);
    int rnk = RankMat(eMat);
    if (rnk == n - 1)
      nb_free++;
    int test = blk.IncrementShow();
    if (test == -1)
      break;
  }
  return nb_free;
}

template <typename T, typename Tgroup>
int CTYP_GetNbAutom(TypeCtypeExch<T> const &TheCtypeArr, [[maybe_unused]] std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tint = typename Tgroup::Tint;
#ifdef TIMINGS
  MicrosecondTime time;
#endif

  int nbRow = TheCtypeArr.eMat.rows();
  int n_edge = 2 * nbRow;
  std::vector<std::vector<unsigned int>> ListGen =
      LinPolytopeAntipodalIntegral_Automorphism(TheCtypeArr.eMat, os);
#ifdef TIMINGS
  os << "|LinPolytopeAntipodal_Automorphism|=" << time << "\n";
#endif

  std::vector<Tidx> v(n_edge);
  std::vector<Telt> ListGenPerm;
  for (auto &eGen : ListGen) {
    for (int i_edge = 0; i_edge < n_edge; i_edge++)
      v[i_edge] = eGen[i_edge];
    ListGenPerm.push_back(Telt(v));
  }
  Tgroup GRP(ListGenPerm, n_edge);
  Tint e_size = GRP.size();
  int nb_autom = UniversalScalarConversion<int, Tint>(e_size);
  return nb_autom;
}

template <typename T, typename Tgroup>
StructuralInfo CTYP_GetStructuralInfo(TypeCtypeExch<T> const &TheCtypeArr,
                                      std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
#ifdef TIMINGS
  MicrosecondTime time;
#endif

  DataCtypeFacet<T, Tidx> data = CTYP_GetConeInformation<T, Tidx>(TheCtypeArr, os);
#ifdef TIMINGS
  os << "|GetNumberFreeVectors|=" << time << "\n";
#endif
  int nb_triple = data.nb_triple;
  int nb_ineq = data.nb_ineq;
  int nb_ineq_after_crit = data.nb_ineq_after_crit;

  int nb_free = CTYP_GetNumberFreeVectors(TheCtypeArr);
#ifdef TIMINGS
  os << "|GetNumberFreeVectors|=" << time << "\n";
#endif

  int nb_autom = CTYP_GetNbAutom<T, Tgroup>(TheCtypeArr, os);

#ifdef TIMINGS
  os << "|NumberAutomorphism|=" << time << "\n";
#endif
  return {nb_triple, nb_ineq, nb_ineq_after_crit, nb_free, nb_autom};
}

struct TypeIndex {
  size_t iProc;
  int idxMatrix;
  int iAdj;
};

template <typename T> struct PairExch {
  TypeCtypeExch<T> eCtype;
  TypeIndex eIndex;
};

template <typename T>
std::ostream &operator<<(std::ostream &os, TypeCtypeExch<T> const &obj) {
  int nbRow = obj.eMat.rows();
  int nbCol = obj.eMat.cols();
  os << " " << nbRow << " " << nbCol;
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int iCol = 0; iCol < nbCol; iCol++)
      os << " " << obj.eMat(iRow, iCol);
  return os;
}

template <typename T>
bool operator==(TypeCtypeExch<T> const &obj1, TypeCtypeExch<T> const &obj2) {
  int nbRow1 = obj1.eMat.rows();
  int nbRow2 = obj2.eMat.rows();
  if (nbRow1 != nbRow2)
    return false;
  int nbCol1 = obj1.eMat.cols();
  int nbCol2 = obj2.eMat.cols();
  if (nbCol1 != nbCol2)
    return false;
  for (int iRow = 0; iRow < nbRow1; iRow++)
    for (int iCol = 0; iCol < nbCol1; iCol++)
      if (obj1.eMat(iRow, iCol) != obj2.eMat(iRow, iCol))
        return false;
  return true;
}

template <typename T>
bool operator==(TypeCtypeAdjExch<T> const &obj1,
                TypeCtypeAdjExch<T> const &obj2) {
  int nbRow1 = obj1.eMat.rows();
  int nbRow2 = obj2.eMat.rows();
  if (nbRow1 != nbRow2)
    return false;
  int nbCol1 = obj1.eMat.cols();
  int nbCol2 = obj2.eMat.cols();
  if (nbCol1 != nbCol2)
    return false;
  for (int iRow = 0; iRow < nbRow1; iRow++)
    for (int iCol = 0; iCol < nbCol1; iCol++)
      if (obj1.eMat(iRow, iCol) != obj2.eMat(iRow, iCol))
        return false;
  return true;
}

std::ostream &operator<<(std::ostream &os, TypeIndex const &obj) {
  os << obj.iProc << " " << obj.idxMatrix << " " << obj.iAdj;
  return os;
}

namespace std {
template <typename T> struct less<TypeCtypeExch<T>> {
  bool operator()(TypeCtypeExch<T> const &eTPE1,
                  TypeCtypeExch<T> const &eTPE2) const {
    int nbRow = eTPE1.eMat.rows();
    int nbCol = eTPE1.eMat.cols();
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 0; iCol < nbCol; iCol++) {
        if (eTPE1.eMat(iRow, iCol) < eTPE2.eMat(iRow, iCol))
          return true;
        if (eTPE1.eMat(iRow, iCol) > eTPE2.eMat(iRow, iCol))
          return false;
      }
    return false;
  }
};
// clang-format off
}  // namespace std
// clang-format on

namespace boost::serialization {

// TypeCtypeExch
template <class Archive, typename T>
inline void serialize(Archive &ar, TypeCtypeExch<T> &eRecMat,
                      [[maybe_unused]] const unsigned int version) {
  int rows = eRecMat.eMat.rows();
  int cols = eRecMat.eMat.cols();
  ar &make_nvp("rows", rows);
  ar &make_nvp("cols", cols);
  eRecMat.eMat.resize(rows, cols);
  for (int r = 0; r < rows; ++r)
    for (int c = 0; c < cols; ++c)
      ar &make_nvp("val", eRecMat.eMat(r, c));
}

// TypeCtypeExch
template <class Archive>
inline void serialize(Archive &ar, TypeIndex &eTypIdx,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("iProc", eTypIdx.iProc);
  ar &make_nvp("idxMatrix", eTypIdx.idxMatrix);
  ar &make_nvp("iAdj", eTypIdx.iAdj);
}

// PairExch
template <class Archive, typename T>
inline void serialize(Archive &ar, PairExch<T> &ePair,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("perfect", ePair.eCtype);
  ar &make_nvp("index", ePair.eIndex);
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

namespace std {
template <typename Tint> struct hash<TypeCtypeExch<Tint>> {
  std::size_t operator()(const TypeCtypeExch<Tint> &k) const {
    std::size_t h1 = 0;
    int nbRow = k.eMat.rows();
    int nbCol = k.eMat.cols();
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 0; iCol < nbCol; iCol++) {
        Tint eVal = k.eMat(iRow, iCol);
        std::size_t h2 = std::hash<Tint>()(eVal);
        h1 = h2 ^ (h1 << 1);
      }
    return h1;
  }
};
template <typename Tint> struct hash<TypeCtypeAdjExch<Tint>> {
  std::size_t operator()(const TypeCtypeAdjExch<Tint> &k) const {
    std::size_t h1 = 0;
    int nbRow = k.eMat.rows();
    int nbCol = k.eMat.cols();
    for (int iRow = 0; iRow < nbRow; iRow++)
      for (int iCol = 0; iCol < nbCol; iCol++) {
        Tint eVal = k.eMat(iRow, iCol);
        std::size_t h2 = std::hash<Tint>()(eVal);
        h1 = h2 ^ (h1 << 1);
      }
    return h1;
  }
};
// clang-format off
}  // namespace std
// clang-format on

template <typename T>
TypeCtypeExch<T> ParseStringToCtypeExch(std::string const &str) {
  std::vector<std::string> LStr = STRING_Split(str, " ");
  int nRow;
  std::istringstream(LStr[0]) >> nRow;
  int n;
  std::istringstream(LStr[1]) >> n;
  MyMatrix<T> eMat(nRow, n);
  int idx = 2;
  for (int iRow = 0; iRow < nRow; iRow++) {
    for (int iCol = 0; iCol < n; iCol++) {
      T eVal;
      std::istringstream(LStr[idx]) >> eVal;
      eMat(iRow, iCol) = eVal;
      idx++;
    }
  }
  return {eMat};
}

TypeIndex ParseStringToTypeIndex(std::string const &str) {
  std::vector<std::string> LStr = STRING_Split(str, " ");
  size_t iProc;
  std::istringstream(LStr[0]) >> iProc;
  int idxMatrixF;
  std::istringstream(LStr[1]) >> idxMatrixF;
  int iAdj;
  std::istringstream(LStr[1]) >> iAdj;
  return {iProc, idxMatrixF, iAdj};
}

// clang-format off
#endif  // SRC_CTYPE_MPI_CTYPEMPI_TYPES_H_
// clang-format on
