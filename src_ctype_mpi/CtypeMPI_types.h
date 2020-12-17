#ifndef INCLUDE_CTYPE_FUNCTIONALITY
#define INCLUDE_CTYPE_FUNCTIONALITY

#include "MAT_Matrix.h"
#include "Boost_bitset.h"
#include "POLY_cddlib.h"
#include "POLY_c_cddlib.h"
#include "Temp_PolytopeEquiStab.h"
#include <cstring>


#define DEBUG
#define TIMINGS
//#define PRINT_FLIP
//#define PRINT_TRIPLE
//#define PRINT_GET_ADJ


namespace std {
  template <typename T>
  struct hash<MyVector<T>>
  {
    std::size_t operator()(const MyVector<T>& V) const
    {
      std::size_t h1 = 12756;
      int len=V.size();
      for (int i=0; i<len; i++) {
        T eVal = V(i);
        std::size_t h2 = std::hash<T>()(eVal);
        h1 = h2 ^ ( h1 << 1);
      }
      return h1;
    }
  };
}


template<typename T>
MyMatrix<T> ReduceExpandedMatrix(MyMatrix<T> const& M)
{
  int nbRow=M.rows();
  int n=M.cols();
  int nbPair=nbRow / 2;
  MyMatrix<T> Mret(nbPair, n);
  auto is_sign_ok=[&](int const& iRow) -> bool {
    for (int i=0; i<n; i++) {
      T eVal = M(iRow,i);
      if (eVal != 0) {
        return eVal > 0;
      }
    }
  };
  int iPair=0;
  for (int iRow=0; iRow<nbRow; iRow++) {
    if (is_sign_ok(iRow)) {
      for (int i=0; i<n; i++)
        Mret(iPair,i) = M(iRow,i);
      iPair++;
    }
  }
  return Mret;
}






template<typename T>
struct TypeCtypeExch {
  MyMatrix<T> eMat;
};


template<typename T>
std::vector<MyMatrix<T>> CTYP_GetBasis(int n)
{
  std::vector<MyMatrix<T>> ListSymmMat;
  for (int i=0; i<n; i++) {
    MyMatrix<T> eMat = ZeroMatrix<T>(n,n);
    eMat(i,i) = 1;
    ListSymmMat.push_back(eMat);
  }
  for (int i=0; i<n; i++)
    for (int j=i+1; j<n; j++) {
      MyMatrix<T> eMat = ZeroMatrix<T>(n,n);
      eMat(i,j) = 1;
      eMat(j,i) = 1;
      ListSymmMat.push_back(eMat);
    }
  return ListSymmMat;
}

using Tidx = int8_t;

struct triple {
  Tidx i;
  Tidx j;
  Tidx k;
};


template<typename T>
bool CheckCoveringParityClasses(MyMatrix<T> const& M)
{
  size_t n_rows = M.rows();
  size_t n_cols = M.cols();
  std::vector<int> ListStatus(n_rows, 0);
  for (size_t i_row=0; i_row<n_rows; i_row++) {
    int pos = -1;
    int e_pow = 1;
    T eTwo = 2;
    for (size_t i=0; i<n_cols; i++) {
      T res_T = ResInt(M(i_row, i), eTwo);
      int res = UniversalTypeConversion<int,T>(res_T);
      pos += res * e_pow;
      e_pow *= 2;
    }
    if (pos == -1)
      return false;
    ListStatus[pos] += 1;
  }
  for (size_t i_row=0; i_row<n_rows; i_row++)
    if (ListStatus[i_row] != 1)
      return false;
  return true;
}




template<typename T>
MyMatrix<T> CTYP_TheFlipping(MyMatrix<T> const& TheCtype, std::vector<triple> const& TheInfo)
{
  size_t n_rows = TheCtype.rows();
  size_t n_cols = TheCtype.cols();
  size_t n_rows_ret = n_rows / 2;
#ifdef PRINT_FLIP
  std::cerr << "CTYP_TheFlipping n_rows=" << n_rows << " n_cols=" << n_cols << "\n";
#endif
  Face ListIchange(n_rows);
  for (auto & e_triple : TheInfo) {
#ifdef PRINT_FLIP
    std::cerr << "e_triple=" << (int)e_triple.i << " , " << (int)e_triple.j << " , " << (int)e_triple.k << "\n";
#endif
    ListIchange[e_triple.i] = 1;
  }
#ifdef PRINT_FLIP
  std::cerr << "ListIchange =";
  for (size_t i_row=0; i_row<n_rows; i_row++)
    std::cerr << " " << ListIchange[i_row];
  std::cerr << "\n";
  for (size_t i_row=0; i_row<n_rows; i_row++) {
    if (ListIchange[i_row] == 1) {
      std::cerr << "Removed line =";
      for (size_t i=0; i<n_cols; i++)
        std::cerr << " " << TheCtype(i_row,i);
      std::cerr << "\n";
    }
  }
#endif
  MyMatrix<T> RetMat(n_rows_ret, n_cols);
  std::vector<T> V(n_cols);
  size_t idx=0;
  Face IsAssigned(n_rows_ret);
  auto insert_if_signok=[&]() -> void {
    for (size_t i=0; i<n_cols; i++) {
      T eVal=V[i];
      if (eVal != 0) {
        if (eVal > 0) {
          int pos = -1;
          int e_pow = 1;
          T eTwo = 2;
          for (size_t i=0; i<n_cols; i++) {
            T res_T = ResInt(V[i], eTwo);
            int res = UniversalTypeConversion<int,T>(res_T);
            pos += res * e_pow;
            e_pow *= 2;
          }
          if (IsAssigned[pos] == 0) {
            for (size_t i_col=0; i_col<n_cols; i_col++)
              RetMat(idx, i_col) = V[i_col];
            idx++;
            IsAssigned[pos]=1;
          }
        }
        return;
      }
    }
  };
  for (auto & e_triple : TheInfo) {
    Tidx j = e_triple.j;
    Tidx k = e_triple.k;
    //
    for (size_t i_col=0; i_col<n_cols; i_col++)
      V[i_col] = -TheCtype(j, i_col) + TheCtype(k, i_col);
    insert_if_signok();
    //
    for (size_t i_col=0; i_col<n_cols; i_col++)
      V[i_col] =  TheCtype(j, i_col) - TheCtype(k, i_col);
    insert_if_signok();
  }
  for (size_t i_row=0; i_row<n_rows; i_row++) {
    if (ListIchange[i_row] == 0 && i_row % 2 == 0) {
      for (size_t i_col=0; i_col<n_cols; i_col++)
        RetMat(idx, i_col) =  TheCtype(i_row, i_col);
      idx++;
    }
  }
#ifdef DEBUG
  if (n_rows_ret != idx) {
    std::cerr << "Incoherence in idx, n_rows_ret. idx=" << idx << " n_rows_ret=" << n_rows_ret << "\n";
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


template<typename T>
std::pair<std::vector<triple>, std::vector<int8_t>> CTYP_GetListTriple(MyMatrix<T> const& TheCtype)
{
  int n_edge = TheCtype.rows();
  int n_cols = TheCtype.cols();
#ifdef PRINT_TRIPLE
  std::cerr << "n_edge=" << n_edge << " n_cols=" << n_cols << "\n";
#endif
  std::vector<triple> ListTriples;
  std::vector<int8_t> MappingVect(n_edge * n_edge, -1);
  auto get_position=[&](MyVector<T> const& eV, Tidx start_idx) -> Tidx {
    auto get_nature=[&](int8_t pos) -> bool {
      for (Tidx i_col=0; i_col<n_cols; i_col++)
        if (TheCtype(pos, i_col) != eV(i_col))
          return false;
      return true;
    };
    auto get_value=[&]() -> Tidx {
      int pos = -1;
      int e_pow = 1;
      T eTwo = 2;
      for (int i=0; i<n_cols; i++) {
        T res_T = ResInt(eV(i), eTwo);
        int res = UniversalTypeConversion<int,T>(res_T);
        pos += res * e_pow;
        e_pow *= 2;
      }
#ifdef PRINT_TRIPLE
      std::cerr << "eV =";
      for (int i=0; i<n_cols; i++)
        std::cerr << " " << eV(i);
      std::cerr << "\n";
      std::cerr << "Found pos=" << pos << "\n";
#endif
      if (pos == -1)
        return -1;
      if (get_nature(2*pos))
        return 2*pos;
      if (get_nature(2*pos+1))
        return 2*pos+1;
      return -1;
    };
    int pos = get_value();
    if (pos > start_idx)
      return pos;
    return -1;
  };
  for (int8_t i=0; i<n_edge; i++)
    for (int8_t j=i+1; j<n_edge; j++) {
#ifdef PRINT_TRIPLE
      std::cerr << "i=" << (int)i << " j=" << (int)j << "\n";
#endif
      MyVector<T> eDiff(n_cols);
      for (int8_t i_col=0; i_col<n_cols; i_col++)
        eDiff(i_col) = - TheCtype(i, i_col) - TheCtype(j,i_col);
#ifdef PRINT_TRIPLE
      std::cerr << "We have eDiff\n";
#endif
      int8_t k = get_position(eDiff, j);
#ifdef PRINT_TRIPLE
      std::cerr << "k=" << (int)k << "\n";
#endif
      if (k != -1) {
        ListTriples.push_back({i,j,k});
        ListTriples.push_back({j,k,i});
        ListTriples.push_back({k,i,j});
        //
        MappingVect[i * n_edge + j] = k;
        MappingVect[j * n_edge + i] = k;
        //
        MappingVect[i * n_edge + k] = j;
        MappingVect[k * n_edge + i] = j;
        //
        MappingVect[j * n_edge + k] = i;
        MappingVect[k * n_edge + j] = i;
      }
    }
  return {std::move(ListTriples), std::move(MappingVect)};
}


template<typename T>
MyMatrix<T> ExpressMatrixForCType(MyMatrix<T> const& M)
{
  int n = M.cols();
  int nbRow = M.rows();
#ifdef PRINT_EXPRESS
  std::cerr << "n=" << n << " nbRow=" << nbRow << "\n";
#endif
  MyMatrix<T> Mret(2*nbRow, n);
#ifdef DEBUG
  std::vector<int> ListStatus(nbRow,0);
#endif
  for (int iRow=0; iRow<nbRow; iRow++) {
#ifdef PRINT_EXPRESS
    std::cerr << "iRow=" << iRow << "/" << nbRow << "\n";
#endif
    int pos = -1;
    int e_pow = 1;
    T eTwo = 2;
    for (int i=0; i<n; i++) {
      T res_T = ResInt(M(iRow,i), eTwo);
      int res = UniversalTypeConversion<int,T>(res_T);
#ifdef PRINT_EXPRESS
      std::cerr << "  i=" << i << " M(iRow,i)=" << M(iRow,i) << " res_T=" << res_T << " res=" << res << " e_pow=" << e_pow << "\n";
#endif
      pos += res * e_pow;
      e_pow *= 2;
    }
#ifdef PRINT_EXPRESS
    std::cerr << "  pos=" << pos << "\n";
#endif
#ifdef DEBUG
    ListStatus[pos] += 1;
#endif
    for (int i=0; i<n; i++) {
      Mret(2*pos  , i) =  M(iRow,i);
      Mret(2*pos+1, i) = -M(iRow,i);
    }
  }
#ifdef DEBUG
  for (int i=0; i<nbRow; i++)
    if (ListStatus[i] != 1) {
      std::cerr << "Consistency error at i=" << i << "\n";
      throw TerminalException{1};
    }
#endif
  return Mret;
}




template<typename T>
std::vector<TypeCtypeExch<T>> CTYP_GetAdjacentCanonicCtypes(TypeCtypeExch<T> const& TheCtypeArr)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
#ifdef PRINT_GET_ADJ
  std::cerr << "CTYP_GetAdjacentCanonicCtypes, step 1\n";
#endif
  MyMatrix<T> TheCtype = ExpressMatrixForCType(TheCtypeArr.eMat);


#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
#endif
#ifdef PRINT_GET_ADJ
  std::cerr << "CTYP_GetAdjacentCanonicCtypes, step 2\n";
#endif
  std::pair<std::vector<triple>, std::vector<int8_t>> PairTriple = CTYP_GetListTriple(TheCtype);


#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
#endif
#ifdef PRINT_GET_ADJ
  std::cerr << "CTYP_GetAdjacentCanonicCtypes, step 3\n";
#endif
  int8_t n = TheCtype.cols();
  int8_t tot_dim = n*(n+1) / 2;
  auto ComputeInequality=[&](MyVector<T> const& V1, MyVector<T> const& V2) -> MyVector<T> {
    MyVector<T> TheVector(tot_dim);
    int idx=0;
    for (int8_t i=0; i<n; i++) {
      TheVector(idx) = V1(i) * V1(i) - V2(i) * V2(i);
      idx++;
    }
    for (int8_t i=0; i<n; i++)
      for (int8_t j=i+1; j<n; j++) {
        // Factor 2 removed for simplification and faster code.
        TheVector(idx) = V1(i) * V1(j) - V2(i) * V2(j);
        idx++;
      }
    return TheVector;
  };
  std::unordered_map<MyVector<T>, std::vector<triple>> Tot_map;
  auto FuncInsertInequality=[&](int8_t i, int8_t j, int8_t k) -> void {
    MyVector<T> V1(n), V2(n);
    for (int8_t i_col=0; i_col<n; i_col++) {
      V1(i_col) = 2 * TheCtype(k, i_col) + TheCtype(i, i_col);
      V2(i_col) = TheCtype(i, i_col);
    }
    MyVector<T> TheVector = ComputeInequality(V1, V2);
    triple TheInfo = {i,j,k};
    std::vector<triple>& list_trip = Tot_map[TheVector];
    list_trip.push_back(TheInfo);
  };
  for (auto & e_triple : PairTriple.first)
    FuncInsertInequality(e_triple.i, e_triple.j, e_triple.k);


#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
#endif
  int n_edge = TheCtype.rows();
  // We apply here the 3 dimensional criterion for feasibility of C-type switches
  auto TestApplicabilityCriterion_with_e=[&](triple const& e_triple, int8_t const& e) -> bool {
    int8_t i = e_triple.i;
    int8_t j = e_triple.j;
    int8_t k = e_triple.k;
    //
    // testing e
    if (e == i || e == j || e == k)
      return false;
    //
    // getting f and testing it
    int8_t f = PairTriple.second[i * n_edge + e];
    if (f == -1 || f == j || f == k)
      return false;
    //
    // getting g and testing it
    int8_t g = PairTriple.second[j * n_edge + e];
    if (g == -1 || g == f || g == i || g == k)
      return false;
    //
    // getting h and testing it
    int8_t h = PairTriple.second[i * n_edge + g];
    if (h == -1 || h == f || h == e || h == j || h == k)
      return false;
    //
    // testing presence of {j,f,h}
    int8_t h2 = PairTriple.second[j * n_edge + f];
    if (h2 != h)
      return false;
    //
    // We have the 7-uple.
    return true;
  };
  auto TestApplicabilityCriterion=[&](triple const& e_triple) -> bool {
    for (int8_t e=0; e<n_edge; e++)
      if (TestApplicabilityCriterion_with_e(e_triple, e))
        return true;
    return false;
  };
  std::vector<int8_t> ListResultCriterion(n_edge * n_edge, 0);
  for (auto & e_triple : PairTriple.first) {
    if (TestApplicabilityCriterion(e_triple)) {
      int8_t i = e_triple.i;
      int8_t j = e_triple.j;
      int8_t k = e_triple.k;
      ListResultCriterion[i * n_edge + k] = 1;
      ListResultCriterion[k * n_edge + j] = 1;
    }
  }
  auto TestFeasibilityListTriple=[&](std::vector<triple> const& list_triple) -> bool {
    for (auto & e_triple : list_triple)
      if (ListResultCriterion[e_triple.i * n_edge + e_triple.j] == 1)
        return false;
    return true;
  };
  // erasing the inequalities that are sure to be redundant.
  for (auto & kv : Tot_map) {
    if (!TestFeasibilityListTriple(kv.second))
      Tot_map.erase(kv.first);
  }


#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time5 = std::chrono::system_clock::now();
#endif
#ifdef PRINT_GET_ADJ
  std::cerr << "CTYP_GetAdjacentCanonicCtypes, step 4\n";
#endif
  size_t n_ineq = Tot_map.size();
  MyMatrix<T> ListInequalities(n_ineq, tot_dim);
  //  std::cerr << " n_ineq=" << n_ineq << " tot_dim=" << (int)tot_dim << "\n";
  std::vector<std::vector<triple>> ListInformations;
  size_t i_ineq=0;
  for (auto & kv : Tot_map) {
    for (int8_t i_col=0; i_col<tot_dim; i_col++)
      ListInequalities(i_ineq, i_col) = kv.first(i_col);
    i_ineq++;
    ListInformations.push_back(std::move(kv.second));
  }
#ifdef PRINT_GET_ADJ
  std::cerr << "CTYP_GetAdjacentCanonicCtypes, step 5\n";
#endif


#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time6 = std::chrono::system_clock::now();
#endif
  std::vector<int> ListIrred = cbased_cdd::RedundancyReductionClarkson(ListInequalities);
  //  std::vector<int> ListIrred = cdd::RedundancyReductionClarkson(ListInequalities);
#ifdef PRINT_GET_ADJ
  std::cerr << "|ListIrred|=" << ListIrred.size() << "\n";
#endif


#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time7 = std::chrono::system_clock::now();
#endif
  std::vector<TypeCtypeExch<T>> ListCtype;
  for (auto & e_int : ListIrred) {
    MyMatrix<T> FlipMat = CTYP_TheFlipping(TheCtype, ListInformations[e_int]);
    //    std::cerr << "FlipMat=\n";
    //    WriteMatrix(std::cerr, FlipMat);
    MyMatrix<T> CanMat = LinPolytopeAntipodalIntegral_CanonicForm(FlipMat);
    //    std::cerr << "CanMat=\n";
    //    WriteMatrix(std::cerr, CanMat);
    ListCtype.push_back({std::move(CanMat)});
  }
#ifdef PRINT_GET_ADJ
  std::cerr << "CTYP_GetAdjacentCanonicCtypes, step 7\n";
#endif


#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time8 = std::chrono::system_clock::now();
  std::cerr << "|ExpressMatrixForCType|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
  std::cerr << "|CTYP_GetListTriple|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
  std::cerr << "|Insert inequalities|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
  std::cerr << "|Criterion Ineq Drop|=" << std::chrono::duration_cast<std::chrono::microseconds>(time5 - time4).count() << "\n";
  std::cerr << "|ListInformations|=" << std::chrono::duration_cast<std::chrono::microseconds>(time6 - time5).count() << "\n";
  std::cerr << "|RedundancyReductionClarkson|=" << std::chrono::duration_cast<std::chrono::microseconds>(time7 - time6).count() << "\n";
  std::cerr << "|Flip + Canonic|=" << std::chrono::duration_cast<std::chrono::microseconds>(time8 - time7).count() << "\n";
#endif
  return ListCtype;
}


struct TypeIndex {
  size_t iProc;
  int idxMatrix;
  int iAdj;
};


template<typename T>
struct PairExch {
  TypeCtypeExch<T> eCtype;
  TypeIndex eIndex;
};

template<typename T>
void PairExch_to_vectorchar(PairExch<T> const& eP, std::vector<char> & eV)
{
  int nbRow=eP.eCtype.eMat.rows();
  int nbCol=eP.eCtype.eMat.cols();
  char* ptr_o = eV.data();
  //
  std::memcpy(ptr_o, (char*)(&nbRow), sizeof(int)); ptr_o += sizeof(int);
  std::memcpy(ptr_o, (char*)(&nbCol), sizeof(int)); ptr_o += sizeof(int);
  //
  for (int i=0; i<nbRow; i++)
    for (int j=0; j<nbCol; j++) {
      std::memcpy(ptr_o, (char*)(&eP.eCtype.eMat(i,j)), sizeof(T)); ptr_o += sizeof(T);
    }
  //
  std::memcpy(ptr_o, (char*)(&eP.eIndex.iProc), sizeof(size_t)); ptr_o += sizeof(size_t);
  std::memcpy(ptr_o, (char*)(&eP.eIndex.idxMatrix), sizeof(int)); ptr_o += sizeof(int);
  std::memcpy(ptr_o, (char*)(&eP.eIndex.iAdj), sizeof(int)); ptr_o += sizeof(int);
}

template<typename T>
PairExch<T> vectorchar_to_PairExch(std::vector<char> const& eV)
{
  const char* ptr_i=eV.data();
  int nbRow, nbCol;
  std::memcpy((char*)(&nbRow), ptr_i, sizeof(int)); ptr_i += sizeof(int);
  std::memcpy((char*)(&nbCol), ptr_i, sizeof(int)); ptr_i += sizeof(int);
  //
  MyMatrix<T> eMat(nbRow, nbCol);
  for (int i=0; i<nbRow; i++)
    for (int j=0; j<nbCol; j++) {
      T eVal;
      std::memcpy((char*)(&eVal), ptr_i, sizeof(T)); ptr_i += sizeof(T);
      eMat(i, j) = eVal;
    }
  //
  size_t iProc;
  int idxMatrix;
  int iAdj;
  std::memcpy((char*)(&iProc), ptr_i, sizeof(size_t)); ptr_i += sizeof(size_t);
  std::memcpy((char*)(&idxMatrix), ptr_i, sizeof(int)); ptr_i += sizeof(int);
  std::memcpy((char*)(&iAdj), ptr_i, sizeof(int)); ptr_i += sizeof(int);
  return {{eMat}, {iProc, idxMatrix, iAdj}};
}




template<typename T>
std::ostream& operator<<(std::ostream& os, TypeCtypeExch<T> const& obj)
{
  int nbRow=obj.eMat.rows();
  int nbCol=obj.eMat.cols();
  os << " " << nbRow << " " << nbCol;
  for (int iRow=0; iRow<nbRow; iRow++)
    for (int iCol=0; iCol<nbCol; iCol++)
      os << " " << obj.eMat(iRow, iCol);
  return os;
}


template<typename T>
bool operator==(TypeCtypeExch<T> const& obj1, TypeCtypeExch<T> const& obj2)
{
  int nbRow1=obj1.eMat.rows();
  int nbRow2=obj2.eMat.rows();
  if (nbRow1 != nbRow2)
    return false;
  int nbCol1=obj1.eMat.cols();
  int nbCol2=obj2.eMat.cols();
  if (nbCol1 != nbCol2)
    return false;
  for (int iRow=0; iRow<nbRow1; iRow++)
    for (int iCol=0; iCol<nbCol1; iCol++)
      if (obj1.eMat(iRow, iCol) != obj2.eMat(iRow, iCol))
        return false;
  return true;
}



std::ostream& operator<<(std::ostream& os, TypeIndex const& obj)
{
  os << obj.iProc << " " << obj.idxMatrix << " " << obj.iAdj;
  return os;
}






namespace std {
  template<typename T>
  struct less<TypeCtypeExch<T>> {
    bool operator()(TypeCtypeExch<T> const& eTPE1, TypeCtypeExch<T> const& eTPE2) const
    {
      int nbRow=eTPE1.eMat.rows();
      int nbCol=eTPE1.eMat.cols();
      for (int iRow=0; iRow<nbRow; iRow++)
        for (int iCol=0; iCol<nbCol; iCol++) {
          if (eTPE1.eMat(iRow,iCol) < eTPE2.eMat(iRow,iCol))
            return true;
          if (eTPE1.eMat(iRow,iCol) > eTPE2.eMat(iRow,iCol))
            return false;
        }
      return false;
    }
  };
}




namespace boost { namespace serialization {
    // TypeCtypeExch
    template<class Archive, typename T>
      inline void serialize(Archive & ar,
                            TypeCtypeExch<T> & eRecMat,
                            const unsigned int version)
      {
        int rows = eRecMat.eMat.rows();
        int cols = eRecMat.eMat.cols();
        ar & make_nvp("rows", rows);
        ar & make_nvp("cols", cols);
        eRecMat.eMat.resize(rows, cols);
        for (int r = 0; r < rows; ++r)
          for (int c = 0; c < cols; ++c)
            ar & make_nvp("val", eRecMat.eMat(r,c));
      }

    // TypeCtypeExch
    template<class Archive>
      inline void serialize(Archive & ar,
                            TypeIndex & eTypIdx,
                            const unsigned int version)
      {
        ar & make_nvp("iProc", eTypIdx.iProc);
        ar & make_nvp("idxMatrix", eTypIdx.idxMatrix);
        ar & make_nvp("iAdj", eTypIdx.iAdj);
      }

    // PairExch
    template<class Archive, typename T>
      inline void serialize(Archive & ar,
                            PairExch<T> & ePair,
                            const unsigned int version)
      {
        ar & make_nvp("perfect", ePair.eCtype);
        ar & make_nvp("index"  , ePair.eIndex);
      }

  }}

namespace std {
  template <typename Tint>
  struct hash<TypeCtypeExch<Tint>>
  {
    std::size_t operator()(const TypeCtypeExch<Tint>& k) const
    {
      std::size_t h1 = 0;
      int nbRow=k.eMat.rows();
      int nbCol=k.eMat.cols();
      for (int iRow=0; iRow<nbRow; iRow++)
        for (int iCol=0; iCol<nbCol; iCol++) {
          Tint eVal = k.eMat(iRow,iCol);
          std::size_t h2 = std::hash<Tint>()(eVal);
          h1 = h2 ^ ( h1 << 1);
        }
      return h1;
    }
  };
}




template<typename T>
TypeCtypeExch<T> ParseStringToCtypeExch(std::string const& str)
{
  std::vector<std::string> LStr = STRING_Split(str, " ");
  int nRow;
  std::istringstream(LStr[0]) >> nRow;
  int n;
  std::istringstream(LStr[1]) >> n;
  MyMatrix<T> eMat(nRow,n);
  int idx=2;
  for (int iRow=0; iRow<nRow; iRow++) {
    for (int iCol=0; iCol<n; iCol++) {
      T eVal;
      std::istringstream(LStr[idx]) >> eVal;
      eMat(iRow,iCol) = eVal;
      idx++;
    }
  }
  return {eMat};
}

TypeIndex ParseStringToTypeIndex(std::string const& str)
{
  std::vector<std::string> LStr = STRING_Split(str, " ");
  size_t iProc;
  std::istringstream(LStr[0]) >> iProc;
  int idxMatrixF;
  std::istringstream(LStr[1]) >> idxMatrixF;
  int iAdj;
  std::istringstream(LStr[1]) >> iAdj;
  return {iProc, idxMatrixF, iAdj};
}





#endif
