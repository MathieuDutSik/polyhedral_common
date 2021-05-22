#ifndef TEMP_POLYTOPE_EQUI_STAB
#define TEMP_POLYTOPE_EQUI_STAB

#include "GRAPH_bliss.h"
#include "GRAPH_traces.h"
#include "MAT_Matrix.h"
#include "Basic_string.h"
#include "Basic_file.h"
#include "GRAPH_GraphicalFunctions.h"
#include "COMB_Combinatorics_elem.h"
#include "MAT_MatrixInt.h"
#include "Boost_bitset.h"
#include "PERM_Fct.h"




#undef USE_BLISS
#define USE_TRACES

#define USE_PAIRS


//#define DEBUG
//#define TIMINGS


//
// The templatized functions
//

template<bool is_symmetric>
inline typename std::enable_if<is_symmetric,size_t>::type weightmatrix_get_nb(size_t nbRow)
{
  return (nbRow * (nbRow+1)) / 2;
}

template<bool is_symmetric>
inline typename std::enable_if<(not is_symmetric),size_t>::type weightmatrix_get_nb(size_t nbRow)
{
  return nbRow * nbRow;
}

template<bool is_symmetric>
inline typename std::enable_if<is_symmetric,size_t>::type weightmatrix_last_idx(size_t nbRow, size_t iRow)
{
  return iRow + 1;
}

template<bool is_symmetric>
inline typename std::enable_if<(not is_symmetric),size_t>::type weightmatrix_last_idx(size_t nbRow, size_t iRow)
{
  return nbRow;
}

template<bool is_symmetric>
inline typename std::enable_if<is_symmetric,size_t>::type weightmatrix_idx(size_t nbRow, size_t iRow, size_t iCol)
{
  if (iCol <= iRow) {
    return (iRow * (iRow + 1)) / 2 + iCol;
  } else {
    return (iCol * (iCol + 1)) / 2 + iRow;
  }
}

template<bool is_symmetric>
inline typename std::enable_if<(not is_symmetric),size_t>::type weightmatrix_idx(size_t nbRow, size_t iRow, size_t jRow)
{
  return iRow + nbRow * jRow;
}

//
// The template traits
//
template<class T, typename Enable = void>
struct is_vector {
  static bool const value = false;
};

template<class T>
struct is_vector<std::vector<T> > {
  static bool const value = true;
};

//
// The generation function
//
template<typename T>
inline typename std::enable_if<is_vector<T>::value,T>::type GetSymmGenerateValue(int const& rVal)
{
  using Tval = typename T::value_type;
  Tval eVal = rVal;
  T eVect;
  eVect.push_back(eVal);
  return eVect;
}

template<typename T>
inline typename std::enable_if<(not is_vector<T>::value),T>::type GetSymmGenerateValue(int const& rVal)
{
  T eVal = rVal;
  return eVal;
}





template<bool is_symmetric, typename T, typename Tidx_value_impl = int16_t>
struct WeightMatrix {
public:
  using Tidx_value = Tidx_value_impl;
  WeightMatrix()
  {
    nbRow=-1;
  }
  WeightMatrix(size_t const& inpNbRow) : nbRow(inpNbRow)
  {
    size_t nb = weightmatrix_get_nb<is_symmetric>(nbRow);
    TheMat.resize(nb);
  }
  WeightMatrix(size_t const& INP_nbRow, std::vector<Tidx_value> const& INP_TheMat, std::vector<T> const& INP_ListWeight) : nbRow(INP_nbRow), TheMat(INP_TheMat), ListWeight(INP_ListWeight)
  {
  }
  template<typename F>
  WeightMatrix(size_t const& _nbRow, F f) : nbRow(_nbRow)
  {
#ifdef TIMINGS
    std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
    TheMat.resize(nbRow * nbRow);
    std::unordered_map<T, Tidx_value> ValueMap;
    int idxWeight=0;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      size_t last_idx = weightmatrix_last_idx<is_symmetric>(nbRow, iRow);
      for (size_t iCol=0; iCol<last_idx; iCol++) {
        T val = f(iRow,iCol);
        Tidx_value & idx = ValueMap[val];
        if (idx == 0) {
          idxWeight++;
          idx = idxWeight;
          ListWeight.push_back(val);
        }
        Tidx_value pos_val = idx - 1;
        size_t pos = weightmatrix_idx<is_symmetric>(nbRow, iRow, iCol);
        TheMat[pos] = pos_val;
      }
    }
#ifdef TIMINGS
    std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
    std::cerr << "|WeightMatrix(nbRow,f)|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  }
  template<typename F1, typename F2>
  WeightMatrix(size_t const& _nbRow, F1 f1, F2 f2) : nbRow(_nbRow)
  {
#ifdef TIMINGS
    std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
    TheMat.resize(nbRow * nbRow);
    std::unordered_map<T, Tidx_value> ValueMap;
    int idxWeight=0;
    for (size_t iRow=0; iRow<nbRow; iRow++) {
      f1(iRow);
      size_t last_idx = weightmatrix_last_idx<is_symmetric>(nbRow, iRow);
      for (size_t iCol=0; iCol<last_idx; iCol++) {
        T val = f2(iCol);
        Tidx_value & idx = ValueMap[val];
        if (idx == 0) {
          idxWeight++;
          idx = idxWeight;
          ListWeight.push_back(val);
        }
        Tidx_value pos_val = idx - 1;
        size_t pos = weightmatrix_idx<is_symmetric>(nbRow, iRow, iCol);
        TheMat[pos] = pos_val;
      }
    }
#ifdef TIMINGS
    std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
    std::cerr << "|WeightMatrix(nbRow,f1,f2)|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  }
  WeightMatrix(WeightMatrix<is_symmetric,T> const& eMat)
  {
    nbRow = eMat.nbRow;
    ListWeight = eMat.ListWeight;
    TheMat = eMat.TheMat;
  }
  WeightMatrix<is_symmetric,T> operator=(WeightMatrix<is_symmetric,T> const& eMat)
  {
    nbRow = eMat.nbRow;
    ListWeight = eMat.ListWeight;
    TheMat = eMat.TheMat;
    return *this;
  }
  ~WeightMatrix()
  {
  }
  // Below is lighter stuff
  size_t rows(void) const
  {
    return nbRow;
  }
  size_t GetWeightSize(void) const
  {
    return ListWeight.size();
  }
  Tidx_value GetValue(size_t const& iRow, size_t const& iCol) const
  {
    size_t idx = weightmatrix_idx<is_symmetric>(nbRow, iRow, iCol);
    return TheMat[idx];
  }
  void intDirectAssign(size_t const& iRow, size_t const& iCol, Tidx_value const& pos)
  {
    size_t idx = weightmatrix_idx<is_symmetric>(nbRow, iRow, iCol);
    TheMat[idx]=pos;
  }
  void SetWeight(std::vector<T> const & inpWeight)
  {
    ListWeight = inpWeight;
  }
  std::vector<T> const& GetWeight() const
  {
    return ListWeight;
  }
  void ReorderingOfWeights(std::vector<Tidx_value> const& gListRev)
  {
    size_t nbEnt=ListWeight.size();
#ifdef DEBUG
    size_t siz=gListRev.size();
    if (nbEnt != siz) {
      std::cerr << "We should have nbEnt = siz\n";
      std::cerr << "nbEnt=" << nbEnt << "\n";
      std::cerr << "siz=" << siz << "\n";
      throw TerminalException{1};
    }
#endif
    size_t nb = weightmatrix_get_nb<is_symmetric>(nbRow);
    for (size_t idx=0; idx<nb; idx++) {
      Tidx_value eValue=TheMat[idx];
      Tidx_value nValue=gListRev[eValue];
      TheMat[idx]=nValue;
    }
    std::vector<T> NewListWeight(nbEnt);
    for (size_t iEnt=0; iEnt<nbEnt; iEnt++) {
      Tidx_value nEnt=gListRev[iEnt];
      NewListWeight[nEnt]=ListWeight[iEnt];
    }
    ListWeight = NewListWeight;
  }
  // Some sophisticated functionalities
  void ReorderingSetWeight()
  {
    std::map<T, int> ValueMap;
    size_t nbEnt=ListWeight.size();
    for (size_t i_w=0; i_w<ListWeight.size(); i_w++)
      ValueMap[ListWeight[i_w]] = i_w;
    std::vector<Tidx_value> g(nbEnt);
    size_t idx=0;
    for (auto& kv : ValueMap) {
      Tidx_value pos = kv.second;
      g[pos] = idx;
      idx++;
    }
    ReorderingOfWeights(g);
#ifdef DEBUG
    for (size_t iEnt=1; iEnt<nbEnt; iEnt++) {
      if (ListWeight[iEnt-1] >= ListWeight[iEnt]) {
        std::cerr << "ERROR: The ListWeightB is not increasing at iEnt=" << iEnt << "\n";
        throw TerminalException{1};
      }
    }
#endif
  }
  Tidx_value ReorderingSetWeight_specificPosition(Tidx_value specificPosition)
  {
    std::map<T, int> ValueMap;
    size_t nbEnt=ListWeight.size();
    for (size_t i_w=0; i_w<ListWeight.size(); i_w++)
      ValueMap[ListWeight[i_w]] = i_w;
    std::vector<Tidx_value> g(nbEnt);
    size_t idx=0;
    for (auto& kv : ValueMap) {
      Tidx_value pos = kv.second;
      g[pos] = idx;
      idx++;
    }
    ReorderingOfWeights(g);
#ifdef DEBUG
    for (size_t iEnt=1; iEnt<nbEnt; iEnt++) {
      if (ListWeight[iEnt-1] >= ListWeight[iEnt]) {
        std::cerr << "ERROR: The ListWeightB is not increasing at iEnt=" << iEnt << "\n";
        throw TerminalException{1};
      }
    }
#endif
    if (specificPosition == -1)
      return -1;
    return g[specificPosition];
  }
  WeightMatrix<true, T, Tidx_value> GetSymmetricWeightMatrix() const
  {
    size_t siz=ListWeight.size();
    size_t nb = nbRow * (2*nbRow + 1);
    std::vector<Tidx_value> RET_TheMat(nb);
    auto set_entry=[&](size_t i, size_t j, Tidx_value pos) -> void {
      size_t idx = weightmatrix_idx<true>(2*nbRow, i, j);
      RET_TheMat[idx] = pos;
    };
    for (size_t iRow=0; iRow<nbRow; iRow++)
      for (size_t jRow=0; jRow<nbRow; jRow++) {
        Tidx_value pos=GetValue(iRow, jRow);
        set_entry(iRow, jRow + nbRow, pos);
      }
    for (size_t iRow=0; iRow<nbRow; iRow++)
      for (size_t jRow=0; jRow<=iRow; jRow++) {
        set_entry(iRow, jRow, siz);
        set_entry(iRow + nbRow, jRow + nbRow, siz+1);
      }
    // Now the list of weights
    std::unordered_set<T> setWeight;
    std::vector<T> RET_ListWeight = ListWeight;
    for (auto& eWei : ListWeight)
      setWeight.insert(eWei);
    int iVal=1;
    for (int j=0; j<2; j++) {
      while(true) {
        T genVal = GetSymmGenerateValue<T>(iVal);
        if (setWeight.count(genVal) == 0) {
          setWeight.insert(genVal);
          RET_ListWeight.push_back(genVal);
          break;
        }
        iVal++;
      }
    }
    return WeightMatrix<true,T>(2*nbRow, RET_TheMat, RET_ListWeight);
  }
private:
  size_t nbRow;
  std::vector<Tidx_value> TheMat;
  std::vector<T> ListWeight;
};

//
// Hashes and invariants
//

namespace std {
  template <bool is_symmetric,typename T, typename Tidx_value>
  struct hash<WeightMatrix<is_symmetric,T,Tidx_value>>
  {
    std::size_t operator()(WeightMatrix<is_symmetric,T,Tidx_value> const& WMat) const
    {
      auto combine_hash=[](size_t & seed, size_t new_hash) -> void {
        seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      };
      std::vector<T> const& ListWeight = WMat.GetWeight();
      size_t nbWei = ListWeight.size();
      size_t nbRow = WMat.rows();
      std::vector<int> ListAttDiag(nbWei, 0);
      std::vector<int> ListAttOff(nbWei, 0);
      for (size_t iRow=0; iRow<nbRow; iRow++) {
        Tidx_value pos = WMat.GetValue(iRow, iRow);
        ListAttDiag[pos]++;
      }
      for (size_t iRow=0; iRow<nbRow; iRow++)
        for (size_t jRow=0; jRow<weightmatrix_last_idx<is_symmetric>(nbRow, iRow); jRow++) {
          if (iRow != jRow) {
            Tidx_value pos = WMat.GetValue(iRow, jRow);
            ListAttOff[pos]++;
          }
        }
      size_t seed = 0;
      for (size_t iWei=0; iWei<nbWei; iWei++) {
        if (ListAttDiag[iWei] > 0) {
          size_t e_hash1 = std::hash<T>()(ListWeight[iWei]);
          size_t e_hash2 = std::hash<int>()(ListAttDiag[iWei]);
          combine_hash(seed, e_hash1);
          combine_hash(seed, e_hash2);
        }
        if (ListAttOff[iWei] > 0) {
          size_t e_hash1 = std::hash<T>()(ListWeight[iWei]);
          size_t e_hash2 = std::hash<int>()(ListAttOff[iWei]);
          combine_hash(seed, e_hash1);
          combine_hash(seed, e_hash2);
        }
      }
      return seed;
    }
  };
}

template<typename T, typename Tidx_value>
size_t GetLocalInvariantWeightMatrix(WeightMatrix<true,T,Tidx_value> const& WMat, Face const& eSet)
{
  size_t n=eSet.size();
  size_t nbVert=eSet.count();
  std::vector<size_t> eList;
  int set_val;
  // We consider the set of smallest size which gain us speed.
  if (2 * nbVert < n) {
    eList.resize(nbVert);
    size_t aRow=eSet.find_first();
    for (size_t i=0; i<nbVert; i++) {
      eList[i]=aRow;
      aRow=eSet.find_next(aRow);
    }
    set_val = 1;
  } else {
    eList.resize(n - nbVert);
    size_t idx = 0;
    for (size_t i=0; i<n; i++) {
      if (eSet[i] == 0) {
        eList[idx] = i;
        idx++;
      }
    }
    set_val = 0;
  }
  size_t nbWeight=WMat.GetWeightSize();
  std::vector<int> eInv(3 * nbWeight + 1, 0);
  for (auto & aVert : eList) {
    for (size_t i=0; i<n; i++) {
      Tidx_value iWeight=WMat.GetValue(aVert, i);
      if (eSet[i] != set_val) {
        eInv[iWeight]++;
      } else {
        if (i != aVert) {
          eInv[iWeight + nbWeight]++;
        } else {
          eInv[iWeight + 2 * nbWeight]++;
        }
      }
    }
  }
  eInv[3 * nbWeight] = nbVert;
  return std::hash<std::vector<int>>()(eInv);
}

template<bool is_symmetric, typename T, typename Tidx_value>
inline typename std::enable_if<is_totally_ordered<T>::value,size_t>::type GetInvariantWeightMatrix(WeightMatrix<is_symmetric, T, Tidx_value> const& WMat)
{
  static_assert(is_totally_ordered<T>::value, "Requires T to be totally ordered");
  size_t nbVert=WMat.rows();
  std::vector<T> const& ListWeight=WMat.GetWeight();
  size_t nbWeight=ListWeight.size();
  std::vector<int> ListAttDiag(nbWeight,0);
  std::vector<int> ListAttOff(nbWeight,0);
  for (size_t iVert=0; iVert<nbVert; iVert++)
    for (size_t jVert=0; jVert<weightmatrix_last_idx<is_symmetric>(nbVert,iVert); jVert++) {
      Tidx_value iWeight=WMat.GetValue(iVert, jVert);
      if (iVert != jVert) {
        ListAttOff[iWeight]++;
      } else {
        ListAttDiag[iWeight]++;
      }
    }
  std::vector<int> eList = SortingPerm<T,int>(ListWeight);
  std::vector<int> ListAtt(2*nbWeight);
  std::vector<T> ListWeight_B(nbWeight);
  for (size_t iWeight=0; iWeight<nbWeight; iWeight++) {
    size_t jWeight=eList[iWeight];
    ListAtt[iWeight] = ListAttDiag[jWeight];
    ListAtt[iWeight + nbWeight] = ListAttOff[jWeight];
    ListWeight_B[iWeight] = ListWeight[jWeight];
  }
  size_t hash1 = std::hash<std::vector<int>>()(ListAtt);
  size_t hash2 = std::hash<std::vector<T>>()(ListWeight_B);
  size_t hash = hash1 + (hash2 << 6) + (hash2 >> 2);
  return hash;
}



//
// Reading and Writing
//

template<bool is_symmetric, typename T, typename Tidx_value>
void PrintWeightedMatrix(std::ostream &os, WeightMatrix<is_symmetric,T,Tidx_value> const&WMat)
{
  size_t siz=WMat.GetWeightSize();
  size_t nbRow=WMat.rows();
  os << "nbRow=" << nbRow << "  Weights=[";
  std::vector<int> ListValues(siz,0);
  for (size_t iRow=0; iRow<nbRow; iRow++)
    for (size_t iCol=0; iCol<nbRow; iCol++) {
      Tidx_value eVal=WMat.GetValue(iRow, iCol);
      ListValues[eVal]++;
    }
  std::vector<T> const& ListWeight=WMat.GetWeight();
  for (size_t i=0; i<siz; i++) {
    if (i>0)
      os << ", ";
    os << "(" << ListWeight[i] << "," << ListValues[i] << ")";
  }
  os << "]\n";
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    for (size_t iCol=0; iCol<nbRow; iCol++) {
      Tidx_value eVal=WMat.GetValue(iRow, iCol);
      os << " " << eVal;
    }
    os << "\n";
  }
}



template<bool is_symmetric, typename T, typename Tidx_value>
void PrintWeightedMatrixGAP(std::ostream &os, WeightMatrix<is_symmetric,T,Tidx_value> const&WMat)
{
  std::vector<T> const& ListWeight=WMat.GetWeight();
  size_t nbRow=WMat.rows();
  os << "[";
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    if (iRow > 0)
      os << ",\n";
    os << "[";
    for (size_t iCol=0; iCol<nbRow; iCol++) {
      Tidx_value eVal=WMat.GetValue(iRow, iCol);
      T eWei=ListWeight[eVal];
      if (iCol > 0)
	os << ", ";
      os << eWei;
    }
    os << "]";
  }
  os << "]";
}



template<bool is_symmetric, typename T, typename Tidx_value>
void PrintWeightedMatrixNoWeight(std::ostream &os, WeightMatrix<is_symmetric,T,Tidx_value> &WMat)
{
  size_t siz=WMat.GetWeightSize();
  os << "nbWeight=" << siz << "\n";
  size_t nbRow=WMat.rows();
  os << "nbRow=" << WMat.rows() << "\n";
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    for (size_t iCol=0; iCol<nbRow; iCol++) {
      Tidx_value eVal=WMat.GetValue(iRow, iCol);
      os << " " << eVal;
    }
    os << "\n";
  }
}

template<typename T>
WeightMatrix<false, T> ReadWeightedMatrix(std::istream &is)
{
  size_t nbRow;
  is >> nbRow;
  WeightMatrix<false, T> WMat(nbRow);
  size_t nbEnt=0;
  int eVal;
  for (size_t iRow=0; iRow<nbRow; iRow++)
    for (size_t jRow=0; jRow<nbRow; jRow++) {
      is >> eVal;
      WMat.intDirectAssign(iRow, jRow, eVal);
      if (eVal> nbEnt)
	nbEnt=eVal;
    }
  nbEnt++;
  std::vector<T> ListWeight;
  T eVal_T;
  for (size_t iEnt=0; iEnt<nbEnt; iEnt++) {
    is >> eVal_T;
    ListWeight.push_back(eVal_T);
  }
  WMat.SetWeight(ListWeight);
  return WMat;
}


//
// Other unsorted functions
//



template<bool is_symmetric, typename T, typename Tidx_value>
bool RenormalizeWeightMatrix(WeightMatrix<is_symmetric,T,Tidx_value> const& WMatRef, WeightMatrix<is_symmetric,T,Tidx_value> & WMat2)
{
  size_t nbRow=WMatRef.rows();
  size_t nbRow2=WMat2.rows();
  if (nbRow != nbRow2)
    return false;
  size_t nbEnt=WMatRef.GetWeightSize();
  size_t nbEnt2=WMat2.GetWeightSize();
  if (nbEnt != nbEnt2)
    return false;
  std::vector<T> const& ListWeightRef=WMatRef.GetWeight();
  std::vector<T> const& ListWeight=WMat2.GetWeight();
  std::vector<Tidx_value> gListRev(nbEnt);
  std::unordered_map<T,Tidx_value> map;
  for (Tidx_value i=0; i<Tidx_value(nbEnt); i++)
    map[ListWeight[i]] = i+1;
  for (size_t i=0; i<nbEnt; i++) {
    Tidx_value jFound = map[ListWeightRef[i]];
    if (jFound == 0)
      return false;
    gListRev[jFound - 1] = i;
  }
  WMat2.ReorderingOfWeights(gListRev);
#ifdef DEBUG
  std::vector<T> const& ListWeight1=WMatRef.GetWeight();
  std::vector<T> const& ListWeight2=WMat2.GetWeight();
  for (size_t iEnt=0; iEnt<nbEnt; iEnt++) {
    if (ListWeight1[iEnt] == ListWeight2[iEnt]) {
      std::cerr << "ERROR: The reordering failed\n";
      throw TerminalException{1};
    }
  }
#endif
  return true;
}









int GetNeededPower(int nb)
{
  int h=0;
  int eExpo=1;
  while(true) {
    if (nb < eExpo)
      return h;
    h++;
    eExpo *= 2;
  }
}


inline void GetBinaryExpression(int eVal, size_t h, std::vector<int> & eVect)
{
  int eWork, eExpo, eExpoB, res;
  eWork=eVal;
  eExpo=1;
  for (size_t i=0; i<h; i++) {
    eExpoB=eExpo*2;
    res=eWork % eExpoB;
    eVect[i]=res/eExpo;
    eExpo=eExpoB;
    eWork=eWork - res;
  }
}


#ifdef USE_PAIRS
/* Unfortunately, the use of pairs while allowing for a graph with
   a smaller number of vertices gives a running time that is sometimes
   larger. This happened with BLISS but not with TRACES.
*/

/*
  We are doing the coloring of the graph in following way.
  Every vertex corresponds to N vertices.
  Those N vertices are made so that they are all adjacent between those.
  (v_1, ..., v_N)
  We can encode much information between those vertices:
  --- v_i adjacent to w_i or not.
  --- v_i adjacent to w_j or not for i != j for some (i,j) \n S
  It is needed that if an automorphism happens then they map
  N-uple V to W completely.
  The way to obtain that is to ensure that the set of pairs S
  define a graph so that for each vertex i in {1, ...., N} there exist
  a j such that j\notin S.
  ---
  For N even we need to remove pairs {0,1}, {2,3}, ..., {N-2,N-1}
    Write N = 2K then nb_pair = 2K + 2K(2K-1) / 2 - K = K + K(2K-1) = 2 K*K
  For N odd we need to remove pairs {0,1}, {2,3}, ..., {N-3, N-2}, {N-2,N-1}
    Write N = 2K + 1 then nb_pair = 2K + 1 + (2K+1) 2K / 2 - (K+1) = K + (2K+1) K = 2 (K+1) K

 */
int GetNeededN(int nb_color)
{
  int N=1;
  while (true) {
    int res = N % 2;
    int nb_pair;
    int K = N / 2;
    if (res == 0) {
      nb_pair = 2 * K * K;
    } else {
      if (N == 1)
        nb_pair = 1;
      else
        nb_pair = 2 * (K + 1) * K;
    }
    int e_pow = 1 << nb_pair; // Computes 2^nb_pair
    if (e_pow >= nb_color)
      return N;
    N++;
  }
}

std::vector<int> GetListPair(int N, int nb_color)
{
  if (N == 1)
    return {0,0};
  int K = N / 2;
  int e_pow = GetNeededPower(nb_color);
  std::vector<int> V;
  int idx=0;
  for (int i=0; i<N; i++) {
    V.push_back(i);
    V.push_back(i);
    idx++;
    if (idx == e_pow)
      return V;
  }
  for (int i=0; i<N; i++)
    for (int j=i+2; j<N; j++) {
      V.push_back(i);
      V.push_back(j);
      idx++;
      if (idx == e_pow)
        return V;
    }
  for (int i=0; i<K-1; i++) {
    V.push_back(2*i+1);
    V.push_back(2*i+2);
    idx++;
    if (idx == e_pow)
      return V;
  }
  return V;
}





template<typename T, typename Tidx_value>
size_t get_total_number_vertices(WeightMatrix<true, T, Tidx_value> const& WMat)
{
  size_t nbWei=WMat.GetWeightSize();
  size_t nbMult=nbWei+2;
  size_t hS=GetNeededN(nbMult);
  size_t nbRow=WMat.rows();
  size_t nbVert=nbRow + 2;
  size_t nbVertTot=nbVert * hS;
#ifdef DEBUG
  std::cerr << "nbWei=" << nbWei << " nbMult=" << nbMult << " hS=" << hS << " nbRow=" << nbRow << " nbVertTot=" << nbVertTot << "\n";
#endif
  return nbVertTot;
}


template<typename T, typename Fcolor, typename Fadj, typename Tidx_value>
void GetGraphFromWeightedMatrix_color_adj(WeightMatrix<true, T, Tidx_value> const& WMat, Fcolor f_color, Fadj f_adj)
{
  size_t nbWei=WMat.GetWeightSize();
  size_t nbMult=nbWei+2;
  size_t hS=GetNeededN(nbMult);
  std::vector<int> V = GetListPair(hS, nbMult);
  size_t e_pow = V.size() / 2;
#ifdef DEBUG
  std::cerr << "nbWei=" << nbWei << " nbMult=" << nbMult << " hS=" << hS << " e_pow=" << e_pow << "\n";
  for (size_t i_pow=0; i_pow<e_pow; i_pow++) {
    std::cerr << "i_pow=" << i_pow << "  (" << V[2*i_pow] << " | " << V[2*i_pow+1] << ")\n";
  }
#endif
  size_t nbRow=WMat.rows();
  size_t nbVert=nbRow + 2;
  for (size_t iVert=0; iVert<nbVert; iVert++)
    for (size_t iH=0; iH<hS; iH++) {
      size_t aVert=iVert + nbVert*iH;
      f_color(aVert,iH);
    }
  for (size_t iVert=0; iVert<nbVert; iVert++)
    for (size_t iH=0; iH<hS-1; iH++)
      for (size_t jH=iH+1; jH<hS; jH++) {
	size_t aVert=iVert + nbVert*iH;
	size_t bVert=iVert + nbVert*jH;
	f_adj(aVert, bVert);
	f_adj(bVert, aVert);
      }
  std::vector<int> eVect(e_pow);
  for (size_t iVert=0; iVert<nbVert-1; iVert++)
    for (size_t jVert=iVert+1; jVert<nbVert; jVert++) {
      Tidx_value eVal;
      if (jVert == nbRow+1) {
	if (iVert == nbRow)
	  eVal = nbWei;
	else
	  eVal = nbWei+1;
      } else {
	if (jVert == nbRow)
	  eVal = WMat.GetValue(iVert, iVert);
	else
	  eVal = WMat.GetValue(iVert, jVert);
      }
      GetBinaryExpression(eVal, e_pow, eVect);
      for (size_t i_pow=0; i_pow<e_pow; i_pow++)
	if (eVect[i_pow] == 1) {
          int iH1 = V[2*i_pow];
          int iH2 = V[2*i_pow+1];
	  size_t aVert=iVert + nbVert*iH1;
	  size_t bVert=jVert + nbVert*iH2;
	  f_adj(aVert, bVert);
	  f_adj(bVert, aVert);
          if (iH1 != iH2) {
            aVert=iVert + nbVert*iH2;
            bVert=jVert + nbVert*iH1;
            f_adj(aVert, bVert);
            f_adj(bVert, aVert);
          }
	}
    }
}

#else


template<typename T, typename Tidx_valud>
size_t get_total_number_vertices(WeightMatrix<true, T, Tidx_value> const& WMat)
{
  size_t nbWei=WMat.GetWeightSize();
  size_t nbMult=nbWei+2;
#ifdef DEBUG
  std::cerr << "nbWei=" << nbWei << " nbMult=" << nbMult << "\n";
#endif
  size_t hS=GetNeededPower(nbMult);
#ifdef DEBUG
  std::cerr << "hS=" << hS << "\n";
#endif
  size_t nbRow=WMat.rows();
  size_t nbVert=nbRow + 2;
#ifdef DEBUG
  std::cerr << "nbVert=" << nbVert << "\n";
#endif
  return hS*nbVert;
}


template<typename T, typename Fcolor, typename Fadj, typename Tidx_value>
void GetGraphFromWeightedMatrix_color_adj(WeightMatrix<true, T, Tidx_value> const& WMat, Fcolor f_color, Fadj f_adj)
{
  size_t nbWei=WMat.GetWeightSize();
  size_t nbMult=nbWei+2;
#ifdef DEBUG
  std::cerr << "nbWei=" << nbWei << " nbMult=" << nbMult << "\n";
#endif
  size_t hS=GetNeededPower(nbMult);
  size_t nbRow=WMat.rows();
  size_t nbVert=nbRow + 2;
  for (size_t iVert=0; iVert<nbVert; iVert++)
    for (size_t iH=0; iH<hS; iH++) {
      size_t aVert=iVert + nbVert*iH;
      f_color(aVert,iH);
    }
  for (size_t iVert=0; iVert<nbVert; iVert++)
    for (size_t iH=0; iH<hS-1; iH++)
      for (size_t jH=iH+1; jH<hS; jH++) {
	size_t aVert=iVert + nbVert*iH;
	size_t bVert=iVert + nbVert*jH;
	f_adj(aVert, bVert);
	f_adj(bVert, aVert);
      }
  std::vector<int> eVect(hS);
  for (size_t iVert=0; iVert<nbVert-1; iVert++)
    for (size_t jVert=iVert+1; jVert<nbVert; jVert++) {
      Tidx_value eVal;
      if (jVert == nbRow+1) {
	if (iVert == nbRow)
	  eVal = nbWei;
	else
	  eVal = nbWei+1;
      } else {
	if (jVert == nbRow)
	  eVal = WMat.GetValue(iVert, iVert);
	else
	  eVal = WMat.GetValue(iVert, jVert);
      }
      GetBinaryExpression(eVal, hS, eVect);
      for (size_t iH=0; iH<hS; iH++)
	if (eVect[iH] == 1) {
	  size_t aVert=iVert + nbVert*iH;
	  size_t bVert=jVert + nbVert*iH;
	  f_adj(aVert, bVert);
	  f_adj(bVert, aVert);
	}
    }
}

#endif


template<typename T, typename Tidx_value>
bliss::Graph GetBlissGraphFromWeightedMatrix(WeightMatrix<true, T, Tidx_value> const& WMat)
{
  size_t nbVert = get_total_number_vertices(WMat);
  bliss::Graph g(nbVert);
  auto f_color=[&](size_t iVert, size_t eColor) -> void {
    g.change_color(iVert, eColor);
  };
  auto f_adj=[&](size_t iVert, size_t jVert) -> void {
    g.add_edge(iVert, jVert);
  };
  GetGraphFromWeightedMatrix_color_adj(WMat, f_color, f_adj);
  return g;
}



template<typename T, typename Tgr, typename Tidx_value>
inline typename std::enable_if<(not is_functional_graph_class<Tgr>::value),Tgr>::type GetGraphFromWeightedMatrix(WeightMatrix<true, T, Tidx_value> const& WMat)
{
  size_t nof_vertices=get_total_number_vertices(WMat);
#ifdef DEBUG
  std::cerr << "nof_vertices=" << nof_vertices << "\n";
#endif
  Tgr eGR(nof_vertices);
#ifdef DEBUG
  std::cerr << "eGR built\n";
#endif
  eGR.SetHasColor(true);
  auto f_color=[&](size_t iVert, size_t eColor) -> void {
    eGR.SetColor(iVert, eColor);
  };
  auto f_adj=[&](size_t iVert, size_t jVert) -> void {
    eGR.AddAdjacent(iVert, jVert);
  };
  GetGraphFromWeightedMatrix_color_adj(WMat, f_color, f_adj);
  return eGR;
}



template<typename Tgroup>
Tgroup GetStabilizerBlissGraph(bliss::Graph g)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  bliss::Stats stats;
  std::vector<std::vector<unsigned int>> ListGen;
  std::vector<std::vector<unsigned int>>* h = &ListGen;
  g.find_automorphisms(stats, &report_aut_vectvectint, (void *)h);
  size_t nbVert = g.get_nof_vertices();
  std::vector<Telt> generatorList;
  for (auto & eGen : ListGen) {
    std::vector<Tidx> gList(nbVert);
    for (size_t iVert=0; iVert<nbVert; iVert++)
      gList[iVert]=eGen[iVert];
    generatorList.push_back(Telt(gList));
  }
  return Tgroup(generatorList, nbVert);
}


template<typename Tidx>
std::pair<std::vector<Tidx>, std::vector<Tidx>> GetCanonicalizationVector_KernelBis(int const& nbRow, std::vector<unsigned int> const& cl)
{
  size_t nof_vertices = cl.size();
  std::vector<unsigned int> clR(nof_vertices,-1);
  for (size_t i=0; i<nof_vertices; i++)
    clR[cl[i]]=i;
  //
  size_t nbVert=nbRow+2;
  size_t hS = nof_vertices / nbVert;
#ifdef DEBUG
  if (hS * nbVert != nof_vertices) {
    std::cerr << "Error in the number of vertices\n";
    std::cerr << "hS=" << hS << " nbVert=" << nbVert << " nof_vertices=" << nof_vertices << "\n";
    throw TerminalException{1};
  }
#endif
  std::vector<int> MapVectRev(nbVert,-1);
  std::vector<int> ListStatus(nof_vertices,1);
  int posCanonic=0;
  for (size_t iCan=0; iCan<nof_vertices; iCan++) {
    if (ListStatus[iCan] == 1) {
      int iNative=clR[iCan];
      int iVertNative=iNative % nbVert;
      MapVectRev[posCanonic] = iVertNative;
      for (size_t iH=0; iH<hS; iH++) {
	int uVertNative = iVertNative + nbVert * iH;
	int jCan=cl[uVertNative];
#ifdef DEBUG
	if (ListStatus[jCan] == 0) {
	  std::cerr << "Quite absurd, should not be 0 iH=" << iH << "\n";
	  throw TerminalException{1};
	}
#endif
	ListStatus[jCan] = 0;
      }
      posCanonic++;
    }
  }
  std::vector<Tidx> MapVect2(nbRow, -1), MapVectRev2(nbRow,-1);
  int posCanonicB=0;
  for (size_t iCan=0; iCan<nbVert; iCan++) {
    int iNative=MapVectRev[iCan];
    if (iNative < nbRow) {
      MapVectRev2[posCanonicB] = iNative;
      MapVect2[iNative] = posCanonicB;
      posCanonicB++;
    }
  }
  return {std::move(MapVect2), std::move(MapVectRev2)};
}




// This function takes a matrix and returns the vector
// that canonicalize it.
// This depends on the construction of the graph from GetGraphFromWeightedMatrix
//
template<typename Tgr, typename Tidx>
std::pair<std::vector<Tidx>, std::vector<Tidx>> GetCanonicalizationVector_Kernel(Tgr const& eGR, int const& nbRow)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  unsigned int nof_vertices=eGR.GetNbVert();
  //  PrintStabilizerGroupSizes(std::cerr, eGR);

#ifdef USE_BLISS
  bliss::Graph g=GetBlissGraphFromGraph(eGR);
# ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
# endif
  //  std::cerr << "|GRP|=" << GetStabilizerBlissGraph(g).size << "\n";

  bliss::Stats stats;
  const unsigned int* cl;
  cl=g.canonical_form(stats, &report_aut_void, stderr);
# ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
# endif

#endif
  //
#ifdef USE_TRACES
# ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
# endif
  std::vector<unsigned int> cl = TRACES_GetCanonicalOrdering(eGR);
# ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
# endif
#endif
  //
  std::vector<unsigned int> clR(nof_vertices,-1);
  for (size_t i=0; i<nof_vertices; i++)
    clR[cl[i]]=i;
  //
  size_t nbVert=nbRow+2;
  size_t hS = nof_vertices / nbVert;
#ifdef DEBUG
  if (hS * nbVert != nof_vertices) {
    std::cerr << "Error in the number of vertices\n";
    std::cerr << "hS=" << hS << " nbVert=" << nbVert << " nof_vertices=" << nof_vertices << "\n";
    throw TerminalException{1};
  }
#endif
  std::vector<int> MapVectRev(nbVert,-1);
  std::vector<int> ListStatus(nof_vertices,1);
  int posCanonic=0;
  for (size_t iCan=0; iCan<nof_vertices; iCan++) {
    if (ListStatus[iCan] == 1) {
      int iNative=clR[iCan];
      int iVertNative=iNative % nbVert;
      MapVectRev[posCanonic] = iVertNative;
      for (size_t iH=0; iH<hS; iH++) {
	int uVertNative = iVertNative + nbVert * iH;
	int jCan=cl[uVertNative];
#ifdef DEBUG
	if (ListStatus[jCan] == 0) {
	  std::cerr << "Quite absurd, should not be 0 iH=" << iH << "\n";
	  throw TerminalException{1};
	}
#endif
	ListStatus[jCan] = 0;
      }
      posCanonic++;
    }
  }
  std::vector<Tidx> MapVect2(nbRow, -1), MapVectRev2(nbRow,-1);
  int posCanonicB=0;
  for (size_t iCan=0; iCan<nbVert; iCan++) {
    int iNative=MapVectRev[iCan];
    if (iNative < nbRow) {
      MapVectRev2[posCanonicB] = iNative;
      MapVect2[iNative] = posCanonicB;
      posCanonicB++;
    }
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|GetBlissGraphFromGraph|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
  std::cerr << "|canonical_form|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
  std::cerr << "|Array shuffling|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif
  return {std::move(MapVect2), std::move(MapVectRev2)};
}




// This function takes a matrix and returns the vector
// that canonicalize it.
// This depends on the construction of the graph from GetGraphFromWeightedMatrix
//
template<typename T, typename Tgr, typename Tidx, typename Tidx_value>
std::pair<std::vector<Tidx>, std::vector<Tidx>> GetCanonicalizationVector(WeightMatrix<true, T, Tidx_value> const& WMat)
{
  size_t nbRow=WMat.rows();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif

  Tgr eGR=GetGraphFromWeightedMatrix<T,Tgr>(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetGraphFromWeightedMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  return GetCanonicalizationVector_Kernel<Tgr,Tidx>(eGR, nbRow);
}


template<typename Tgr, typename Tidx>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> GetGroupCanonicalizationVector_Kernel(Tgr const& eGR, int const& nbRow)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  unsigned int nof_vertices=eGR.GetNbVert();

#ifdef USE_BLISS
  std::pair<std::vector<unsigned int>, std::vector<std::vector<unsigned int>>> ePair = BLISS_GetCanonicalOrdering_ListGenerators(eGR);
#endif
  //
#ifdef USE_TRACES
  std::pair<std::vector<unsigned int>, std::vector<std::vector<unsigned int>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators(eGR);
#endif
  //
  std::vector<std::vector<Tidx>> LGen;
  for (auto& eListI : ePair.second) {
    std::vector<Tidx> eListO(nbRow);
    for (Tidx i=0; i<Tidx(nbRow); i++)
      eListO[i] = eListI[i];
    LGen.push_back(eListO);
  }
  //
  std::vector<unsigned int> clR(nof_vertices,-1);
  for (size_t i=0; i<nof_vertices; i++)
    clR[ePair.first[i]]=i;
  //
  size_t nbVert=nbRow+2;
  size_t hS = nof_vertices / nbVert;
#ifdef DEBUG
  if (hS * nbVert != nof_vertices) {
    std::cerr << "Error in the number of vertices\n";
    std::cerr << "hS=" << hS << " nbVert=" << nbVert << " nof_vertices=" << nof_vertices << "\n";
    throw TerminalException{1};
  }
#endif
  std::vector<int> MapVectRev(nbVert,-1);
  std::vector<int> ListStatus(nof_vertices,1);
  int posCanonic=0;
  for (size_t iCan=0; iCan<nof_vertices; iCan++) {
    if (ListStatus[iCan] == 1) {
      int iNative=clR[iCan];
      int iVertNative=iNative % nbVert;
      MapVectRev[posCanonic] = iVertNative;
      for (size_t iH=0; iH<hS; iH++) {
	int uVertNative = iVertNative + nbVert * iH;
	int jCan=ePair.first[uVertNative];
#ifdef DEBUG
	if (ListStatus[jCan] == 0) {
	  std::cerr << "Quite absurd, should not be 0 iH=" << iH << "\n";
	  throw TerminalException{1};
	}
#endif
	ListStatus[jCan] = 0;
      }
      posCanonic++;
    }
  }
  std::vector<Tidx> MapVectRev2(nbRow,-1);
  int posCanonicB=0;
  for (size_t iCan=0; iCan<nbVert; iCan++) {
    int iNative=MapVectRev[iCan];
    if (iNative < nbRow) {
      MapVectRev2[posCanonicB] = iNative;
      posCanonicB++;
    }
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|GetBlissGraphFromGraph|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
  std::cerr << "|canonical_form|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
  std::cerr << "|Array shuffling|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif
  return {std::move(MapVectRev2), std::move(LGen)};
}




// This function takes a matrix and returns the vector
// that canonicalize it.
// This depends on the construction of the graph from GetGraphFromWeightedMatrix
//
template<typename T, typename Tgr, typename Tidx, typename Tidx_value>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> GetGroupCanonicalizationVector(WeightMatrix<true, T, Tidx_value> const& WMat)
{
  size_t nbRow=WMat.rows();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif

  Tgr eGR=GetGraphFromWeightedMatrix<T,Tgr>(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetGraphFromWeightedMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
#ifdef DEBUG
  std::vector<std::vector<Tidx>> LGen = GetGroupCanonicalizationVector_Kernel<Tgr,Tidx>(eGR, nbRow).second;
  for (auto & eGen : LGen) {
    for (size_t i=0; i<nbRow; i++) {
      for (size_t j=0; j<nbRow; j++) {
        int iImg = eGen[i];
        int jImg = eGen[j];
        Tidx_value pos1 = WMat.GetValue(i, j);
        Tidx_value pos2 = WMat.GetValue(iImg, jImg);
        if (pos1 != pos2) {
          std::cerr << "Inconsistency at i=" << i << " j=" <<j << "\n";
          std::cerr << "iImg=" << iImg << " jImg=" << jImg << "\n";
          throw TerminalException{1};
        }
      }
    }
  }
#endif
  return GetGroupCanonicalizationVector_Kernel<Tgr,Tidx>(eGR, nbRow);
}












template<typename T, typename Tgroup, typename Tidx_value>
Tgroup GetStabilizerWeightMatrix(WeightMatrix<true, T, Tidx_value> const& WMat)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  bliss::Stats stats;
  std::vector<std::vector<unsigned int>> ListGen;
  size_t nbRow=WMat.rows();
  GraphBitset eGR=GetGraphFromWeightedMatrix<T,GraphBitset>(WMat);
  bliss::Graph g=GetBlissGraphFromGraph(eGR);
  g.find_automorphisms(stats, &report_aut_vectvectint, (void *)(&ListGen));
  std::vector<Telt> generatorList;
  for (auto & eGen : ListGen) {
    std::vector<Tidx> gList(nbRow);
    for (size_t iVert=0; iVert<nbRow; iVert++) {
      size_t jVert=eGen[iVert];
#ifdef DEBUG
      if (jVert >= nbRow) {
	std::cerr << "jVert is too large\n";
	std::cerr << "jVert=" << jVert << "\n";
	std::cerr << "nbRow=" << nbRow << "\n";
	throw TerminalException{1};
      }
#endif
      gList[iVert]=jVert;
    }
#ifdef DEBUG
    for (size_t iRow=0; iRow<nbRow; iRow++)
      for (size_t jRow=0; jRow<nbRow; jRow++) {
	int iRowI=gList[iRow];
	int jRowI=gList[jRow];
	Tidx_value eVal1=WMat.GetValue(iRow, jRow);
	Tidx_value eVal2=WMat.GetValue(iRowI, jRowI);
	if (eVal1 != eVal2) {
	  std::cerr << "eVal1=" << eVal1 << " eVal2=" << eVal2 << "\n";
	  std::cerr << "Clear error in automorphism computation\n";
	  std::cerr << "AUT iRow=" << iRow << " jRow=" << jRow << "\n";
	  throw TerminalException{1};
	}
      }
#endif
    generatorList.push_back(Telt(gList));
  }
  return Tgroup(generatorList, nbRow);
}


template<typename T, typename Tgroup, typename Tidx_value>
Tgroup GetStabilizerAsymmetricMatrix(WeightMatrix<false, T, Tidx_value> const& WMatI)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  WeightMatrix<true, T> WMatO=WMatI.GetSymmetricWeightMatrix();
  size_t nbSHV=WMatI.rows();
  Tgroup GRP=GetStabilizerWeightMatrix<T,Tgroup>(WMatO);
  std::vector<Telt> ListGenInput = GRP.GeneratorsOfGroup();
  std::vector<Tidx> v(nbSHV);
  std::vector<Telt> ListGen;
  for (auto & eGen : ListGenInput) {
    for (size_t iSHV=0; iSHV<nbSHV; iSHV++)
      v[iSHV]=OnPoints(iSHV, eGen);
    ListGen.push_back(Telt(v));
  }
  return Tgroup(ListGen, nbSHV);
}



std::pair<std::vector<int>, std::vector<int>> GetCanonicalizationFromSymmetrized(std::pair<std::vector<int>, std::vector<int>> const& PairVectSymm)
{
  size_t nbEnt=PairVectSymm.first.size() / 2;
  std::vector<int> MapVect(nbEnt, -1), MapVectRev(nbEnt, -1);
  std::vector<int> ListStatus(2*nbEnt,1);
  size_t jEntCan=0;
  for (size_t iEntCan=0; iEntCan<2*nbEnt; iEntCan++) {
    if (ListStatus[iEntCan] == 1) {
      int iEntNative = PairVectSymm.second[iEntCan];
      int jEntNative = iEntNative % nbEnt;
      MapVectRev[jEntCan] = jEntNative;
      MapVect[jEntNative] = jEntCan;
      for (int iH=0; iH<2; iH++) {
	int iEntNativeB = jEntNative + nbEnt * iH;
	int iEntCanB=PairVectSymm.first[iEntNativeB];
#ifdef DEBUG
	if (ListStatus[iEntCanB] == 0) {
	  std::cerr << "Quite absurd, should not be 0 iH=" << iH << "\n";
	  throw TerminalException{1};
	}
#endif
	ListStatus[iEntCanB] = 0;
      }
      jEntCan++;
    }
  }
  return {std::move(MapVect), std::move(MapVectRev)};
}





template<typename T, typename Tidx_value>
EquivTest<std::vector<unsigned int>> TestEquivalenceWeightMatrix_norenorm(WeightMatrix<true, T, Tidx_value> const& WMat1, WeightMatrix<true, T, Tidx_value> const& WMat2)
{
  //  using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
  Tgr eGR1=GetGraphFromWeightedMatrix<T,Tgr>(WMat1);
  Tgr eGR2=GetGraphFromWeightedMatrix<T,Tgr>(WMat2);
  unsigned int nof_vertices1=eGR1.GetNbVert();
  unsigned int nof_vertices2=eGR2.GetNbVert();
  if (nof_vertices1 != nof_vertices2)
    return {false, {}};
  unsigned int nof_vertices=nof_vertices1;
  size_t nbRow=WMat1.rows();
#ifdef USE_BLISS
  std::vector<unsigned int> cl1 = BLISS_GetCanonicalOrdering(eGR1);
  std::vector<unsigned int> cl2 = BLISS_GetCanonicalOrdering(eGR2);
#endif
#ifdef USE_TRACES
  std::vector<unsigned int> cl1 = TRACES_GetCanonicalOrdering(eGR1);
  std::vector<unsigned int> cl2 = TRACES_GetCanonicalOrdering(eGR2);
#endif
  std::vector<unsigned int> clR2(nof_vertices);
  for (unsigned int i=0; i<nof_vertices; i++)
    clR2[cl2[i]]=i;
  std::vector<unsigned int> TheEquivExp(nof_vertices, -1);
  for (unsigned int iVert=0; iVert<nof_vertices; iVert++) {
    unsigned int jVert = clR2[cl1[iVert]];
#ifdef DEBUG
    unsigned int iBlock = iVert / (nbRow + 2);
    unsigned int jBlock = jVert / (nbRow + 2);
    if (iBlock != jBlock) {
      std::cerr << "Not repecting block structure\n";
      throw TerminalException{1};
    }
#endif
    TheEquivExp[iVert] = jVert;
  }
  for (unsigned int iVert=0; iVert<nof_vertices; iVert++) {
    unsigned int jVert=TheEquivExp[iVert];
    if (eGR1.GetColor(iVert) != eGR2.GetColor(jVert))
      return {false, {}};
  }
  for (unsigned int iVert1=0; iVert1<nof_vertices; iVert1++) {
    unsigned int iVert2=TheEquivExp[iVert1];
    for (unsigned int jVert1=0; jVert1<nof_vertices; jVert1++) {
      unsigned int jVert2=TheEquivExp[jVert1];
      if (eGR1.IsAdjacent(iVert1,jVert1) != eGR2.IsAdjacent(iVert2,jVert2) )
	return {false, {}};
    }
  }
  std::vector<unsigned int> TheEquiv(nbRow);
  for (unsigned int i=0; i<nbRow; i++)
    TheEquiv[i]=TheEquivExp[i];
  return {true, std::move(TheEquiv)};
}

template<typename T, typename Telt, typename Tidx_value>
EquivTest<Telt> TestEquivalenceWeightMatrix_norenorm_perm(WeightMatrix<true, T, Tidx_value> const& WMat1, WeightMatrix<true, T, Tidx_value> const& WMat2)
{
  EquivTest<std::vector<unsigned int>> ePair = TestEquivalenceWeightMatrix_norenorm(WMat1, WMat2);
  size_t len = ePair.TheEquiv.size();
  std::vector<int> eList(len);
  for (size_t i=0; i<len; i++)
    eList[i] = ePair.TheEquiv[i];
  Telt ePerm(eList);
  return {ePair.TheReply, std::move(ePerm)};
}

template<typename T, typename Telt, typename Tidx_value>
EquivTest<Telt> TestEquivalenceWeightMatrix(WeightMatrix<true, T, Tidx_value> const& WMat1, WeightMatrix<true, T, Tidx_value> &WMat2)
{
  bool test=RenormalizeWeightMatrix(WMat1, WMat2);
  if (!test)
    return {false, {}};
  return TestEquivalenceWeightMatrix_norenorm_perm<T,Telt>(WMat1, WMat2);
}





template<typename T, typename Telt, typename Tidx_value>
EquivTest<Telt> GetEquivalenceAsymmetricMatrix(WeightMatrix<false, T, Tidx_value> const& WMat1, WeightMatrix<false, T, Tidx_value> const& WMat2)
{
  using Tidx = typename Telt::Tidx;
  WeightMatrix<true, T, Tidx_value> WMatO1=WMat1.GetSymmetricWeightMatrix();
  WeightMatrix<true, T, Tidx_value> WMatO2=WMat2.GetSymmetricWeightMatrix();
  EquivTest<Telt> eResEquiv=TestEquivalenceWeightMatrix<T,Telt>(WMatO1, WMatO2);
  if (!eResEquiv.TheReply)
    return eResEquiv;
  size_t nbSHV=WMat1.rows();
  std::vector<Tidx> v(nbSHV);
  for (size_t i=0; i<nbSHV; i++)
    v[i]=eResEquiv.TheEquiv.at(i);
  return {true, std::move(Telt(v))};
}






#endif
