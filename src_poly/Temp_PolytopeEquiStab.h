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



using Tidx_value = int16_t;


#undef USE_BLISS
#define USE_TRACES

#define USE_PAIRS


//#define DEBUG
//#define TIMINGS


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





template<bool is_symmetric, typename T>
struct WeightMatrix {
public:
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
    nbRow = eMat.rows();
    ListWeight = eMat.GetWeight();
    TheMat = eMat.TheMat;
  }
  WeightMatrix<is_symmetric,T> operator=(WeightMatrix<is_symmetric,T> const& eMat)
  {
    nbRow = eMat.rows();
    ListWeight = eMat.GetWeight();
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
  std::vector<T> GetWeight() const
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
private:
  size_t nbRow;
  std::vector<Tidx_value> TheMat;
  std::vector<T> ListWeight;
};


namespace std {
  template <bool is_symmetric,typename T>
  struct hash<WeightMatrix<is_symmetric,T>>
  {
    std::size_t operator()(WeightMatrix<is_symmetric,T> const& WMat) const
    {
      auto combine_hash=[](size_t & seed, size_t new_hash) -> void {
        seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      };
      std::vector<T> ListWeight = WMat.GetWeight();
      size_t nbWei = ListWeight.size();
      size_t nbRow = WMat.rows();
      std::vector<int> ListAttDiag(nbWei, 0);
      std::vector<int> ListAttOff(nbWei, 0);
      for (size_t iRow=0; iRow<nbRow; iRow++) {
        Tidx_value pos = WMat.GetValue(iRow, iRow);
        ListAttDiag[pos]++;
      }
      for (size_t iRow=0; iRow<nbRow; iRow++)
        for (size_t jRow=0; jRow<nbRow; jRow++) {
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




template<bool is_symmetric, typename T>
void PrintWeightedMatrix(std::ostream &os, WeightMatrix<is_symmetric,T> const&WMat)
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
  std::vector<T> ListWeight=WMat.GetWeight();
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



template<bool is_symmetric, typename T>
void PrintWeightedMatrixGAP(std::ostream &os, WeightMatrix<is_symmetric,T> const&WMat)
{
  std::vector<T> ListWeight=WMat.GetWeight();
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



template<bool is_symmetric, typename T>
void PrintWeightedMatrixNoWeight(std::ostream &os, WeightMatrix<is_symmetric,T> &WMat)
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



template<bool is_symmetric, typename T>
bool RenormalizeWeightMatrix(WeightMatrix<is_symmetric,T> const& WMatRef, WeightMatrix<is_symmetric,T> & WMat2)
{
  size_t nbRow=WMatRef.rows();
  size_t nbRow2=WMat2.rows();
  if (nbRow != nbRow2)
    return false;
  size_t nbEnt=WMatRef.GetWeightSize();
  size_t nbEnt2=WMat2.GetWeightSize();
  if (nbEnt != nbEnt2)
    return false;
  std::vector<T> ListWeightRef=WMatRef.GetWeight();
  std::vector<T> ListWeight=WMat2.GetWeight();
  std::vector<Tidx_value> gListRev(nbEnt);
  for (size_t i=0; i<nbEnt; i++) {
    Tidx_value jFound=-1;
    for (size_t j=0; j<nbEnt; j++)
      if (ListWeightRef[i] == ListWeight[j])
	jFound=j;
    if (jFound == -1)
      return false;
    gListRev[jFound]=i;
  }
  WMat2.ReorderingOfWeights(gListRev);
#ifdef DEBUG
  std::vector<T> ListWeight1=WMatRef.GetWeight();
  std::vector<T> ListWeight2=WMat2.GetWeight();
  for (size_t iEnt=0; iEnt<nbEnt; iEnt++) {
    if (ListWeight1[iEnt] == ListWeight2[iEnt]) {
      std::cerr << "ERROR: The reordering failed\n";
      throw TerminalException{1};
    }
  }
#endif
  return true;
}




template<typename T>
WeightMatrix<false,T> T_TranslateToMatrix(MyMatrix<T> const& eMat)
{
  size_t nbRow=eMat.rows();
  auto f=[&](size_t iRow, size_t iCol) -> T {
    return eMat(iRow,iCol);
  };
  return WeightMatrix<false,T>(nbRow, f);
}



// The matrices in ListMat do not have to be symmetric.
template<typename T, typename Tint>
WeightMatrix<false,std::vector<T>> T_TranslateToMatrix_ListMat_SHV(std::vector<MyMatrix<T>> const& ListMat, MyMatrix<Tint> const& SHV)
{
  size_t nbRow=SHV.rows();
  size_t n = SHV.cols();
  size_t nbMat=ListMat.size();
  std::vector<MyVector<T>> ListV(nbMat);
  auto f1=[&](size_t iRow) -> void {
    for (size_t iMat=0; iMat<nbMat; iMat++) {
      MyVector<T> V(n);
      for (size_t i=0; i<n; i++) {
        T eVal=0;
        for (size_t j=0; j<n; j++)
          eVal += ListMat[iMat](j,i) * SHV(iRow, j);
        V(i) = eVal;
      }
      ListV[iMat] = V;
    }
  };
  std::vector<T> ListScal(nbMat);
  auto f2=[&](size_t iCol) -> std::vector<T> {
    for (size_t iMat=0; iMat<nbMat; iMat++) {
      T eScal=0;
      for (size_t i=0; i<n; i++)
        eScal += ListV[iMat](i)*SHV(iCol,i);
      ListScal[iMat] = eScal;
    }
    return ListScal;
  };
  return WeightMatrix<false,std::vector<T>>(nbRow, f1, f2);
}



template<typename T, typename Tint>
WeightMatrix<true,T> T_TranslateToMatrix_QM_SHV(MyMatrix<T> const& qMat, MyMatrix<Tint> const& SHV)
{
  size_t nbRow=SHV.rows();
  size_t n=qMat.rows();
  size_t INP_nbRow=nbRow;
  size_t nbPair=nbRow / 2;
  size_t nb = nbPair * (2*nbPair + 1);
  std::vector<Tidx_value> INP_TheMat(nb);
  std::vector<T> INP_ListWeight;
  std::unordered_map<T, Tidx_value> ValueMap;
  Tidx_value idxWeight = 0;
  //
  auto set_entry=[&](size_t iRow, size_t iCol, Tidx_value val) -> void {
    size_t idx = weightmatrix_idx<true>(nbRow, iRow, iCol);
    INP_TheMat[idx] = val;
  };
  for (size_t iPair=0; iPair<nbPair; iPair++) {
    MyVector<T> V(n);
    for (size_t i=0; i<n; i++) {
      T eVal=0;
      for (size_t j=0; j<n; j++)
	eVal += qMat(j,i) * SHV(2*iPair, j);
      V(i) = eVal;
    }
    for (size_t jPair=iPair; jPair<=iPair; jPair++) {
      T eScal=0;
      for (size_t i=0; i<n; i++)
	eScal += V(i)*SHV(2*jPair,i);
      Tidx_value& value1 = ValueMap[eScal];
      if (value1 == 0) { // This is a missing value
        idxWeight++;
        value1 = idxWeight;
        INP_ListWeight.push_back(eScal);
      }
      Tidx_value& value2 = ValueMap[-eScal];
      if (value2 == 0) { // This is a missing value
        idxWeight++;
        value2 = idxWeight;
        INP_ListWeight.push_back(-eScal);
      }
      Tidx_value pos1 = value1 - 1;
      Tidx_value pos2 = value2 - 1;
      set_entry(2*iPair  , 2*jPair  , pos1);
      set_entry(2*iPair+1, 2*jPair  , pos2);
      set_entry(2*iPair  , 2*jPair+1, pos2);
      set_entry(2*iPair+1, 2*jPair+1, pos1);
    }
  }
  return WeightMatrix<true,T>(INP_nbRow, INP_TheMat, INP_ListWeight);
}







template<typename T>
MyMatrix<T> Kernel_GetQmatrix(MyMatrix<T> const& TheEXT)
{
  static_assert(is_ring_field<T>::value, "Requires T to be a field");
  size_t nbRow=TheEXT.rows();
  size_t nbCol=TheEXT.cols();
  MyMatrix<T> QMat(nbCol, nbCol);
  for (size_t iCol=0; iCol<nbCol; iCol++)
    for (size_t jCol=0; jCol<nbCol; jCol++) {
      T eSum=0;
      for (size_t iRow=0; iRow<nbRow; iRow++)
	eSum += TheEXT(iRow, jCol) * TheEXT(iRow, iCol);
      QMat(iCol, jCol)=eSum;
    }
  return Inverse_destroy(QMat);
}

template<typename T>
inline typename std::enable_if<is_ring_field<T>::value,MyMatrix<T>>::type GetQmatrix(MyMatrix<T> const& TheEXT)
{
  return Kernel_GetQmatrix(TheEXT);
}


template<typename T>
inline typename std::enable_if<(not is_ring_field<T>::value),MyMatrix<T>>::type GetQmatrix(MyMatrix<T> const& TheEXT)
{
  using Tfield=typename overlying_field<T>::field_type;
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif

  MyMatrix<Tfield> TheEXT_F = ConvertMatrixUniversal<Tfield,T>(TheEXT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|ConvertMatrixUniversal1|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif

  MyMatrix<Tfield> Q_F = Kernel_GetQmatrix(TheEXT_F);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|Kernel_GetQmatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif

  MyMatrix<Tfield> Q_F_red = RemoveFractionMatrix(Q_F);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|RemoveFractionMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif

  MyMatrix<T> RetMat = ConvertMatrixUniversal<T,Tfield>(Q_F_red);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time5 = std::chrono::system_clock::now();
  std::cerr << "|ConvertMatrixUniversal2|=" << std::chrono::duration_cast<std::chrono::microseconds>(time5 - time4).count() << "\n";
#endif
  return RetMat;
}






template<typename T>
WeightMatrix<true, T> GetSimpleWeightMatrix(MyMatrix<T> const& TheEXT, MyMatrix<T> const& Qmat)
{
  size_t nbRow=TheEXT.rows();
  size_t nbCol=TheEXT.cols();
  MyVector<T> V(nbCol);
  auto f1=[&](size_t iRow) -> void {
    for (size_t iCol=0; iCol<nbCol; iCol++) {
      T eSum=0;
      for (size_t jCol=0; jCol<nbCol; jCol++)
        eSum += Qmat(iCol,jCol) * TheEXT(iRow, jCol);
      V(iCol) = eSum;
    }
  };
  auto f2=[&](size_t jRow) -> T {
    T eSum=0;
    for (size_t iCol=0; iCol<nbCol; iCol++)
      eSum += V(iCol) * TheEXT(jRow, iCol);
    return eSum;
  };
  return WeightMatrix<true,T>(nbRow, f1, f2);
}



// ListMat is assumed to be symmetric
template<typename T>
WeightMatrix<true,std::vector<T>> GetWeightMatrix_ListMat_Subset(MyMatrix<T> const& TheEXT, std::vector<MyMatrix<T>> const& ListMat, Face const& eSubset)
{
#ifdef DEBUG
  for (auto & eMat : ListMat) {
    if (!IsSymmetricMatrix(eMat)) {
      std::cerr << "The matrix eMat should be symmetric\n";
      throw TerminalException{1};
    }
  }
#endif
  size_t nbRow=TheEXT.rows();
  size_t nbCol=TheEXT.cols();
  size_t nMat = ListMat.size();
  //
  MyMatrix<T> MatV(nMat, nbCol);
  std::vector<T> LScal(nMat + 1);
  size_t iRow_stor = 0;
  auto f1=[&](size_t iRow) -> void {
    for (size_t iMat=0; iMat<nMat; iMat++) {
      for (size_t iCol=0; iCol<nbCol; iCol++) {
        T eSum=0;
        for (size_t jCol=0; jCol<nbCol; jCol++)
          eSum += ListMat[iMat](jCol,iCol) * TheEXT(iRow, jCol);
        MatV(iMat, iCol) = eSum;
      }
    }
    iRow_stor = iRow;
  };
  auto f2=[&](size_t jRow) -> std::vector<T> {
    for (size_t iMat=0; iMat<nMat; iMat++) {
      T eSum=0;
      for (size_t iCol=0; iCol<nbCol; iCol++)
        eSum += MatV(iMat, iCol) * TheEXT(jRow, iCol);
      LScal[iMat] = eSum;
    }
    Tidx_value eVal = 0;
    if (iRow_stor == jRow)
      eVal = eSubset[jRow];
    LScal[nMat] = eVal;
    return LScal;
  };
  return WeightMatrix<true,std::vector<T>>(nbRow, f1, f2);
}













template<typename T>
struct WeightMatrixAbs {
  Tidx_value positionZero;
  Face ArrSigns;
  WeightMatrix<true,T> WMat;
};


template<typename T>
WeightMatrixAbs<T> GetSimpleWeightMatrixAntipodal_AbsTrick(MyMatrix<T> const& TheEXT, MyMatrix<T> const& Qmat)
{
  static_assert(is_totally_ordered<T>::value, "Requires T to be a totally ordered field");
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  size_t nbPair=TheEXT.rows();
  size_t nbCol=TheEXT.cols();
  size_t eProd = nbPair * nbPair;
  std::vector<Tidx_value> INP_TheMat(eProd);
  Face ArrSigns(eProd);
  std::vector<T> INP_ListWeight;
  std::unordered_map<T, Tidx_value> ValueMap;
  Tidx_value idxWeight = 0;
  Tidx_value positionZero = -1;
  //
  auto set_entry=[&](size_t iRow, size_t jRow, Tidx_value pos, bool eChg) -> void {
    size_t idx = weightmatrix_idx<true>(nbPair, iRow, jRow);
    INP_TheMat[idx] = pos;
    size_t idx1 = iRow + nbPair*jRow;
    size_t idx2 = jRow + nbPair*iRow;
    ArrSigns[idx1] = eChg;
    ArrSigns[idx2] = eChg;
  };
  MyVector<T> V(nbCol);
  for (size_t iPair=0; iPair<nbPair; iPair++) {
    for (size_t iCol=0; iCol<nbCol; iCol++) {
      T eSum=0;
      for (size_t jCol=0; jCol<nbCol; jCol++)
        eSum += Qmat(iCol,jCol) * TheEXT(iPair, jCol);
      V(iCol) = eSum;
    }
    for (size_t jPair=0; jPair<=iPair; jPair++) {
      T eScal=0;
      for (size_t iCol=0; iCol<nbCol; iCol++)
        eScal += V(iCol) * TheEXT(jPair, iCol);
      bool ChgSign=false;
      if (eScal < 0) {
        eScal = -eScal;
        ChgSign = true;
      }
      Tidx_value& value = ValueMap[eScal];
      if (value == 0) { // This is a missing value
        if (positionZero == -1 && eScal == 0)
          positionZero = idxWeight;
        idxWeight++;
        value = idxWeight;
        INP_ListWeight.push_back(eScal);
      }
      Tidx_value pos = value - 1;
      set_entry(iPair  , jPair  , pos, ChgSign);
    }
  }
  /* This cannot be simplified be a classic constructor WeightMatrix(nbRow,f1,f2)
     because we also need to compute the positionZero and the ArrSigns. */
  WeightMatrix<true,T> WMat(nbPair, INP_TheMat, INP_ListWeight);
#ifdef DEBUG
  std::cerr << "Before positionZero=" << positionZero << "\n";
#endif
  positionZero = WMat.ReorderingSetWeight_specificPosition(positionZero);
#ifdef DEBUG
  std::cerr << "Afeter positionZero=" << positionZero << "\n";
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  return {positionZero, std::move(ArrSigns), std::move(WMat)};
}





template<typename T>
WeightMatrix<true, T> GetSimpleWeightMatrixAntipodal(MyMatrix<T> const& TheEXT, MyMatrix<T> const& Qmat)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  size_t nbPair=TheEXT.rows();
  size_t nbCol=TheEXT.cols();
  size_t INP_nbRow = 2*nbPair;
  std::vector<Tidx_value> INP_TheMat(INP_nbRow * INP_nbRow);
  std::vector<T> INP_ListWeight;
  std::unordered_map<T, Tidx_value> ValueMap;
  Tidx_value idxWeight = 0;
  //
  auto set_entry=[&](size_t iRow, size_t jRow, Tidx_value pos) -> void {
    size_t idx = iRow + INP_nbRow*jRow;
    INP_TheMat[idx] = pos;
  };
  MyVector<T> V(nbCol);
  for (size_t iPair=0; iPair<nbPair; iPair++) {
    for (size_t iCol=0; iCol<nbCol; iCol++) {
      T eSum=0;
      for (size_t jCol=0; jCol<nbCol; jCol++)
        eSum += Qmat(iCol,jCol) * TheEXT(iPair, jCol);
      V(iCol) = eSum;
    }
    for (size_t jPair=0; jPair<=iPair; jPair++) {
      T eSum1=0;
      for (size_t iCol=0; iCol<nbCol; iCol++)
        eSum1 += V(iCol) * TheEXT(jPair, iCol);
      T eSum2 = -eSum1;
      Tidx_value& value1 = ValueMap[eSum1];
      if (value1 == 0) { // This is a missing value
        idxWeight++;
        value1 = idxWeight;
        INP_ListWeight.push_back(eSum1);
      }
      Tidx_value& value2 = ValueMap[eSum2];
      if (value2 == 0) { // This is a missing value
        idxWeight++;
        value2 = idxWeight;
        INP_ListWeight.push_back(eSum2);
      }
      Tidx_value pos1 = value1 - 1;
      Tidx_value pos2 = value2 - 1;
      set_entry(2*iPair  , 2*jPair  , pos1);
      set_entry(2*iPair+1, 2*jPair  , pos2);
      set_entry(2*iPair  , 2*jPair+1, pos2);
      set_entry(2*iPair+1, 2*jPair+1, pos1);
      if (iPair != jPair) {
        set_entry(2*jPair  , 2*iPair  , pos1);
        set_entry(2*jPair+1, 2*iPair  , pos2);
        set_entry(2*jPair  , 2*iPair+1, pos2);
        set_entry(2*jPair+1, 2*iPair+1, pos1);
      }
    }
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetSimpleWeightMatrixAntipodal|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  return WeightMatrix<true,T>(INP_nbRow, INP_TheMat, INP_ListWeight);
}



template<typename T>
WeightMatrix<true, T> GetWeightMatrix(MyMatrix<T> const& TheEXT)
{
  MyMatrix<T> Qmat=GetQmatrix(TheEXT);
  return GetSimpleWeightMatrix(TheEXT, Qmat);
}





template<typename T>
WeightMatrix<true, T> GetWeightMatrixAntipodal(MyMatrix<T> const& TheEXT)
{
  MyMatrix<T> Qmat=GetQmatrix(TheEXT);
  return GetSimpleWeightMatrixAntipodal(TheEXT, Qmat);
}



template<bool is_symmetric, typename T>
WeightMatrix<is_symmetric, std::vector<T>> GetWeightMatrix_ListComm(MyMatrix<T> const& TheEXT, MyMatrix<T> const&GramMat, std::vector<MyMatrix<T>> const& ListComm)
{
  size_t nbRow=TheEXT.rows();
  size_t nbCol=TheEXT.cols();
  size_t nbComm=ListComm.size();
  std::vector<MyMatrix<T>> ListProd;
  ListProd.push_back(GramMat);
  for (size_t iComm=0; iComm<nbComm; iComm++) {
    MyMatrix<T> eProd=ListComm[iComm]*GramMat;
    ListProd.push_back(eProd);
  }
  std::vector<T> eVectSum(nbComm+1);
  auto f=[&](size_t iRow, size_t jRow) -> std::vector<T> {
    for (size_t iMat=0; iMat<=nbComm; iMat++)
      eVectSum[iMat] = 0;
    for (size_t iCol=0; iCol<nbCol; iCol++)
      for (size_t jCol=0; jCol<nbCol; jCol++) {
        T eProd=TheEXT(iRow, iCol) * TheEXT(jRow, jCol);
        for (size_t iMat=0; iMat<=nbComm; iMat++)
          eVectSum[iMat] += eProd * ListProd[iMat](iCol, jCol);
      }
    return eVectSum;
  };
  return WeightMatrix<false,std::vector<T>>(nbRow, f);
}



template<typename T>
WeightMatrix<false,std::vector<T>> GetWeightMatrix_ListMatrix(std::vector<MyMatrix<T>> const& ListMatrix, MyMatrix<T> const& TheEXT)
{
  size_t nbRow=TheEXT.rows();
  size_t nbCol=TheEXT.cols();
  size_t nbMat=ListMatrix.size();
  std::vector<T> eVectScal(nbMat);
  auto f=[&](size_t iRow, size_t jRow) -> std::vector<T> {
    for (size_t iMat=0; iMat<nbMat; iMat++)
      eVectScal[iMat] = 0;
    for (size_t iCol=0; iCol<nbCol; iCol++)
      for (size_t jCol=0; jCol<nbCol; jCol++) {
	T eProd=TheEXT(iRow, iCol) * TheEXT(jRow, jCol);
        for (size_t iMat=0; iMat<nbMat; iMat++)
          eVectScal[iMat] += eProd * ListMatrix[iMat](iCol, jCol);
      }
    return eVectScal;
  };
  return WeightMatrix<false,std::vector<T>>(nbRow, f);
}



template<typename T>
WeightMatrix<true, T> GetWeightMatrixGramMatShort(MyMatrix<T> const& TheGramMat, MyMatrix<int> const& ListShort)
{
  size_t nbShort=ListShort.rows();
  size_t n=TheGramMat.rows();
  MyVector<T> V(n);
  auto f1=[&](size_t iShort) -> void {
    for (size_t i=0; i<n; i++) {
      T eSum = 0;
      for (size_t j=0; j<n; j++)
        eSum += TheGramMat(i,j) * ListShort(iShort, j);
      V(i) = eSum;
    }
  };
  auto f2=[&](size_t jShort) -> T {
    T eScal = 0;
    for (size_t i=0; i<n; i++)
      eScal += V(i) * ListShort(jShort, i);
    return eScal;
  };
  return WeightMatrix<true,T>(nbShort, f1, f2);
}



template<typename T>
WeightMatrix<true, T> GetWeightMatrixGramMatShort_Fast(MyMatrix<T> const& TheGramMat, MyMatrix<int> const& ListShort)
{
  // TO REMOVE
  size_t nbShort=ListShort.rows();
  size_t n=TheGramMat.rows();
  auto GetValue=[&](size_t const&iShort, size_t const&jShort) -> T {
    T eScal=0;
    for (size_t i=0; i<n; i++)
      for (size_t j=0; j<n; j++) {
	int eVal12=ListShort(iShort, i) * ListShort(jShort, j);
	eScal += eVal12 * TheGramMat(i,j);
      }
    return eScal;
  };
  MyMatrix<T> ScalMat(nbShort,nbShort);
  std::set<T> setWeight;
  for (size_t iShort=0; iShort<nbShort; iShort++)
    for (size_t jShort=0; jShort<=iShort; jShort++) {
      T eScal=GetValue(iShort,jShort);
      ScalMat(iShort,jShort)=eScal;
      ScalMat(jShort,iShort)=eScal;
      setWeight.insert(eScal);
    }
  struct PairData {
    T x;
    size_t idx;
  };
  auto comp=[&](PairData const& a, PairData const& b) -> bool {
    if (a.x < b.x)
      return true;
    if (a.x > b.x)
      return false;
    return false;
  };
  std::set<PairData,decltype(comp)> setWeightIdx(comp);
  std::vector<T> INP_ListWeight;
  size_t idx=0;
  for (auto & eX : setWeight) {
    setWeightIdx.insert({eX,idx});
    INP_ListWeight.push_back(eX);
    idx++;
  }
  std::vector<Tidx_value> INP_TheMat(nbShort*nbShort);
  for (size_t iShort=0; iShort<nbShort; iShort++)
    for (size_t jShort=0; jShort<=iShort; jShort++) {
      T eScal=GetValue(iShort,jShort);
      PairData test{eScal,0};
      auto iter=setWeightIdx.find(test);
#ifdef DEBUG
      if (iter == setWeightIdx.end()) {
	std::cerr << "Without a doubt a bug\n";
	throw TerminalException{1};
      }
#endif
      Tidx_value idxret=iter->idx;
      size_t pos1=iShort + nbShort*jShort;
      size_t pos2=jShort + nbShort*iShort;
      INP_TheMat[pos1]=idxret;
      INP_TheMat[pos2]=idxret;
    }
  return WeightMatrix<true,T>(nbShort, INP_TheMat, INP_ListWeight);
}


template<class T, typename Enable = void>
struct is_vector {
  static bool const value = false;
};


template<class T>
struct is_vector<std::vector<T> > {
  static bool const value = true;
};



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



template<typename T>
void GetSymmGenerateValue(int const& rVal, std::vector<T> & eVect)
{
  T eVal = rVal;
  eVect.push_back(eVal);
}



template<typename T>
WeightMatrix<true, T> GetSymmetricWeightMatrix(WeightMatrix<false, T> const& WMatI)
{
  std::set<T> setWeight;
  std::vector<T> ListWeight;
  size_t nbRow=WMatI.rows();
  WeightMatrix<true,T> WMatO(2*nbRow);
  size_t siz=WMatI.GetWeightSize();
  for (size_t iRow=0; iRow<nbRow; iRow++)
    for (size_t jRow=0; jRow<nbRow; jRow++) {
      Tidx_value pos=WMatI.GetValue(iRow, jRow);
      WMatO.intDirectAssign(iRow, jRow+nbRow, pos);
      WMatO.intDirectAssign(jRow+nbRow, iRow, pos);
    }
  for (size_t iRow=0; iRow<nbRow; iRow++)
    for (size_t jRow=0; jRow<nbRow; jRow++) {
      WMatO.intDirectAssign(iRow      , jRow      , siz);
      WMatO.intDirectAssign(iRow+nbRow, jRow+nbRow, siz+1);
    }
  for (auto& eWei : WMatI.GetWeight())
    setWeight.insert(eWei);
  ListWeight=WMatI.GetWeight();
  int iVal=1;
  for (int j=0; j<2; j++) {
    while(true) {
      T genVal = GetSymmGenerateValue<T>(iVal);
      typename std::set<T>::iterator iterTEST=setWeight.find(genVal);
      if (iterTEST == setWeight.end()) {
	setWeight.insert(genVal);
	ListWeight.push_back(genVal);
	break;
      }
      iVal++;
    }
  }
  WMatO.SetWeight(ListWeight);
  return WMatO;
}




template<bool is_symmetric, typename T, typename Tout>
std::vector<Tout> GetLocalInvariantWeightMatrix(WeightMatrix<is_symmetric,T> const&WMat, Face const& eSet)
{
  size_t nbVert=eSet.count();
  std::vector<size_t> eList(nbVert);
  size_t aRow=eSet.find_first();
  for (size_t i=0; i<nbVert; i++) {
    eList[i]=aRow;
    aRow=eSet.find_next(aRow);
  }
  size_t nbWeight=WMat.GetWeightSize();
  std::vector<Tout> eInv(nbWeight,0);
  for (size_t iVert=0; iVert<nbVert; iVert++) {
    size_t aVert=eList[iVert];
    for (size_t jVert=0; jVert<nbVert; jVert++) {
      size_t bVert=eList[jVert];
      size_t iWeight=WMat.GetValue(aVert, bVert);
      eInv[iWeight]++;
    }
  }
  return eInv;
}



template<typename T, typename Tgroup>
WeightMatrix<true, T> WeightMatrixFromPairOrbits(Tgroup const& GRP, std::ostream & os)
{
#ifdef DEBUG
  bool IsDiag;
#endif
  size_t n=GRP.n_act();
  WeightMatrix<true, T> WMat(n);
  for (size_t i=0; i<n; i++)
    for (size_t j=0; j<n; j++)
      WMat.intDirectAssign(i,j,-1);
  auto GetUnset=[&]() -> std::pair<int,int> {
    for (size_t i=0; i<n; i++)
      for (size_t j=0; j<n; j++) {
	Tidx_value eVal=WMat.GetValue(i,j);
	if (eVal == -1) {
	  return {i,j};
	}
      }
    return {-1,-1};
  };
  struct VectorListPair {
    std::vector<std::pair<int,int>> ListWorkingPair;
    size_t nbWorkingPair=0;
  };
  VectorListPair VLP0, VLP1;
  auto FuncInsert=[&](VectorListPair & VLP, std::pair<int,int> const& ePair) -> void {
    if (VLP.nbWorkingPair < VLP.ListWorkingPair.size()) {
      VLP.ListWorkingPair[VLP.nbWorkingPair]=ePair;
    } else {
      VLP.ListWorkingPair.push_back(ePair);
    }
    VLP.nbWorkingPair++;
    int i=ePair.first;
    int j=ePair.second;
    WMat.intDirectAssign(i,j,-2);
  };
  auto FuncInsertIChoice=[&](int const& iChoice, std::pair<int,int> const& ePair) -> void {
    if (iChoice == 0)
      FuncInsert(VLP0, ePair);
    else
      FuncInsert(VLP1, ePair);
  };
  auto SetZero=[&](int const& iChoice) -> void {
    if (iChoice == 0)
      VLP0.nbWorkingPair = 0;
    else
      VLP1.nbWorkingPair = 0;
  };
  auto GetNbWorkingPair=[&](int const& iChoice) ->size_t {
    if (iChoice == 0)
      return VLP0.nbWorkingPair;
    else
      return VLP1.nbWorkingPair;
  };
  auto GetEntry=[&](int const& iChoice, int const& iPair) -> std::pair<int,int> {
    if (iChoice == 0)
      return VLP0.ListWorkingPair[iPair];
    else
      return VLP1.ListWorkingPair[iPair];
  };
  int iChoice=0;
  int iOrbit=0;
  std::vector<T> ListWeight;
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> ListGen = GRP.GeneratorsOfGroup();
  while(true) {
    std::pair<int,int> eStart=GetUnset();
    if (eStart.first == -1)
      break;
#ifdef DEBUG
    IsDiag=false;
    if (eStart.first == eStart.second)
      IsDiag=true;
    os << "iOrbit=" << iOrbit << " eStart=" << eStart.first << " , " << eStart.second << "\n";
    std::cerr << "iOrbit=" << iOrbit << " eStart=" << eStart.first << " , " << eStart.second << "\n";
    std::cerr << "  IsDiag=" << IsDiag << "\n";
#endif
    T insVal=iOrbit;
    ListWeight.push_back(insVal);
    FuncInsertIChoice(iChoice, eStart);
    size_t orbSize=0;
    while(true) {
      int iChoiceB = 1 - iChoice;
      int nbPair=GetNbWorkingPair(iChoice);
      orbSize += nbPair;
      if (nbPair == 0)
	break;
      for (int iPair=0; iPair<nbPair; iPair++) {
	std::pair<int,int> ePair=GetEntry(iChoice,iPair);
	int i=ePair.first;
	int j=ePair.second;
	WMat.intDirectAssign(i,j,iOrbit);
	for (auto & eGen : ListGen) {
	  int iImg = OnPoints(i, eGen);
	  int jImg = OnPoints(j, eGen);
	  auto aInsert=[&](int const& u, int const& v) -> void {
	    int eVal1 = WMat.GetValue(u,v);
	    if (eVal1 == -1)
	      FuncInsertIChoice(iChoiceB, {u, v});
	  };
	  aInsert(iImg, jImg);
	  aInsert(jImg, iImg);
	}
      }
      SetZero(iChoice);
      iChoice=iChoiceB;
    }
#ifdef DEBUG
    os << "     size=" << GRP.size() << " orbSize=" << orbSize << "\n";
#endif
    iOrbit++;
  }
  WMat.SetWeight(ListWeight);
#ifdef DEBUG
  for (size_t i=0; i<n; i++)
    os << "i=" << i << "/" << n << " val=" << WMat.GetValue(i,i) << "\n";
#endif
  return WMat;
}



struct LocalInvInfo {
  int nbDiagCoeff;
  int nbOffCoeff;
  std::vector<int> MapDiagCoeff;
  std::vector<int> MapOffCoeff;
  MyMatrix<int> ListChosenTriple;
  WeightMatrix<true,int> WMatInt;
};



template<typename T, typename Tgroup>
LocalInvInfo ComputeLocalInvariantStrategy(WeightMatrix<true, T> const&WMat, Tgroup const& GRP, std::string const& strat, std::ostream & os)
{
  //  os << "ComputeLocalInvariantStrategy, step 1\n";
  int nbNeed=0;
  bool UsePairOrbit=false;
  std::vector<std::string> LStr=STRING_Split(strat, ":");
  //  os << "ComputeLocalInvariantStrategy, step 2\n";
  if (LStr[0] != "pairinv") {
    std::cerr << "Now we have only pairinv as simple invariant\n";
    throw TerminalException{1};
  }
  //  os << "ComputeLocalInvariantStrategy, step 3\n";
  for (size_t iStr=1; iStr<LStr.size(); iStr++) {
    std::string eStr=LStr[iStr];
    std::vector<std::string> LStrB=STRING_Split(eStr, "_");
    if (LStrB[0] == "target") {
      int nbNeed;
      std::istringstream(LStrB[1]) >> nbNeed;
    }
    if (eStr == "use_pair_orbit") {
      UsePairOrbit=true;
    }
  }
  //  os << "ComputeLocalInvariantStrategy, step 4\n";
  //  os << "nbNeed=" << nbNeed << "\n";
  //  os << "UsePairOrbit=" << UsePairOrbit << "\n";
  //
  WeightMatrix<true,int> WMatInt;
  if (UsePairOrbit) {
    WMatInt = WeightMatrixFromPairOrbits<int,Tgroup>(GRP, os);
  }
  else {
    WMatInt = NakedWeightedMatrix(WMat);
  }
  //  os << "ComputeLocalInvariantStrategy, step 5\n";
  //
  size_t nbRow=WMatInt.rows();
  size_t nbWeight=WMatInt.GetWeightSize();
  //  os << "nbRow=" << nbRow << " nbWeight=" << nbWeight << "\n";
  std::vector<int> StatusDiag(nbWeight,0);
  std::vector<int> StatusOff(nbWeight,0);
  //  os << "ComputeLocalInvariantStrategy, step 5.1\n";
  for (size_t i=0; i<nbRow; i++) {
    size_t iWeight=WMatInt.GetValue(i,i);
    StatusDiag[iWeight]=1;
  }
  //  os << "ComputeLocalInvariantStrategy, step 5.2\n";
  for (size_t i=0; i<nbRow-1; i++)
    for (size_t j=i+1; j<nbRow; j++) {
      size_t iWeight=WMatInt.GetValue(i,j);
      StatusOff[iWeight]=1;
    }
  //  os << "ComputeLocalInvariantStrategy, step 6\n";
  int nbDiagCoeff=VectorSum(StatusDiag);
  int nbOffCoeff=VectorSum(StatusOff);
  std::vector<int> MapDiagCoeff(nbWeight,-1);
  std::vector<int> MapOffCoeff(nbWeight,-1);
  int idxA=0;
  for (size_t i=0; i<nbWeight; i++)
    if (StatusDiag[i] == 1) {
      MapDiagCoeff[i]=idxA;
      idxA++;
    }
  //  os << "ComputeLocalInvariantStrategy, step 7\n";
  int idxB=0;
  for (size_t i=0; i<nbWeight; i++)
    if (StatusOff[i] == 1) {
      MapOffCoeff[i]=idxB;
      idxB++;
    }
  //  os << "ComputeLocalInvariantStrategy, step 8\n";
  LocalInvInfo eInv;
  eInv.nbDiagCoeff=nbDiagCoeff;
  eInv.nbOffCoeff=nbOffCoeff;
  eInv.MapDiagCoeff=MapDiagCoeff;
  eInv.MapOffCoeff=MapOffCoeff;
  eInv.WMatInt=WMatInt;
  //
  int nbNeedTriple=nbNeed - nbDiagCoeff - 2*nbOffCoeff;
  if (nbNeedTriple <= 0)
    return eInv;
  std::set<std::vector<int>> ListEnt;
  std::map<std::vector<int>,int> MapMult;
  auto InsertEntry=[&](std::vector<int> const& rVect) -> void {
    auto iter=ListEnt.find(rVect);
    if (iter == ListEnt.end()) {
      ListEnt.insert(rVect);
      MapMult[rVect]=1;
    } else {
      MapMult[rVect]++;
    }
  };
  std::vector<std::vector<int>> ListIndices{{0,1},{0,2},{1,2}};
  std::vector<int> eVect=BinomialStdvect_First(3);
  std::vector<int> colorVect(3);
  while(true) {
    for (int i=0; i<3; i++) {
      int iIdx=ListIndices[i][0];
      int jIdx=ListIndices[i][1];
      int iVert=eVect[iIdx];
      int jVert=eVect[jIdx];
      int iWeight=WMat.GetValue(iVert,jVert);
      colorVect[i]=iWeight;
    }
    std::sort(colorVect.begin(), colorVect.end());
    InsertEntry(colorVect);
    bool test=BinomialStdvect_Increment(nbRow,3,eVect);
    if (!test)
      break;
  }
  struct pairEntInt {
    std::vector<int> eEnt;
    int mult;
  };
  std::function<bool(pairEntInt const&,pairEntInt const&)> secondComp=[&](pairEntInt const& x,pairEntInt const& y) -> bool {
    if (x.mult < y.mult)
      return true;
    if (x.mult > y.mult)
      return false;
    return x.eEnt < y.eEnt;
  };
  std::set<pairEntInt,decltype(secondComp)> ListEntSecond(secondComp);
  for (auto & uVect : ListEnt) {
    int eMult=MapMult[uVect];
    ListEntSecond.insert({uVect,eMult});
  }
  auto iter=ListEntSecond.end();
  MyMatrix<int> ListChosenTriple(nbNeedTriple,3);
  for (int idx=0; idx<nbNeedTriple; idx++) {
    iter--;
    for (int i=0; i<3; i++)
      ListChosenTriple(idx,i) = iter->eEnt[i];
  }
  eInv.ListChosenTriple=ListChosenTriple;
  return eInv;
}



template<typename Tout>
std::vector<Tout> GetLocalInvariantWeightMatrix_Enhanced(LocalInvInfo const& LocalInv, Face const& eSet)
{
  size_t nbVert=eSet.count();
  size_t nbRow=LocalInv.WMatInt.rows();
  size_t nbVertCompl=nbRow - nbVert;
  std::vector<int> eList(nbVert);
  std::vector<int> eListCompl(nbVertCompl);
  int idx1=0;
  int idx2=0;
  for (size_t i=0; i<nbRow; i++) {
    if (eSet[i] == 1) {
      eList[idx1]=i;
      idx1++;
    }
    else {
      eListCompl[idx2]=i;
      idx2++;
    }
  }
  int nbDiagCoeff=LocalInv.nbDiagCoeff;
  int nbOffCoeff=LocalInv.nbOffCoeff;
  size_t nbTriple=LocalInv.ListChosenTriple.rows();
  std::vector<Tout> eInv(nbDiagCoeff + 2*nbOffCoeff + nbTriple,0);
  int posShift=0;
  for (size_t iVert=0; iVert<nbVert; iVert++) {
    int aVert=eList[iVert];
    int iWeight=LocalInv.WMatInt.GetValue(aVert, aVert);
    int iMap=LocalInv.MapDiagCoeff[iWeight];
    eInv[posShift + iMap]++;
  }
  posShift += nbDiagCoeff;
  for (size_t iVert=0; iVert<nbVert-1; iVert++) {
    int aVert=eList[iVert];
    for (size_t jVert=iVert+1; jVert<nbVert; jVert++) {
      int bVert=eList[jVert];
      int iWeight=LocalInv.WMatInt.GetValue(aVert, bVert);
      int iMap=LocalInv.MapOffCoeff[iWeight];
      eInv[posShift + iMap]++;
    }
  }
  posShift += nbOffCoeff;
  for (size_t iVert=0; iVert<nbVert; iVert++) {
    int aVert=eList[iVert];
    for (size_t jVert=0; jVert<nbVertCompl; jVert++) {
      int bVert=eListCompl[jVert];
      int iWeight=LocalInv.WMatInt.GetValue(aVert,bVert);
      int iMap=LocalInv.MapOffCoeff[iWeight];
      eInv[posShift + iMap]++;
    }
  }
  posShift += nbOffCoeff;
  if (nbTriple == 0)
    return eInv;
  std::vector<int> colorVect(3);
  auto FindITriple=[&]() -> int {
    for (size_t iTriple=0; iTriple<nbTriple; iTriple++) {
      bool IsMatch=true;
      for (int i=0; i<3; i++)
	if (colorVect[i] != LocalInv.ListChosenTriple(iTriple,i))
	  IsMatch=false;
      if (IsMatch)
	return iTriple;
    }
    return -1;
  };
  std::vector<std::vector<int>> ListIndices{{0,1},{0,2},{1,2}};
  std::vector<int> eVect=BinomialStdvect_First(3);
  while(true) {
    for (int i=0; i<3; i++) {
      int iIdx=ListIndices[i][0];
      int jIdx=ListIndices[i][1];
      int iVert=eVect[iIdx];
      int jVert=eVect[jIdx];
      int aVert=eList[iVert];
      int bVert=eList[jVert];
      int iWeight=LocalInv.WMatInt.GetValue(aVert,bVert);
      colorVect[i]=iWeight;
    }
    std::sort(colorVect.begin(), colorVect.end());
    int iMatchTriple=FindITriple();
    if (iMatchTriple != -1)
      eInv[posShift + iMatchTriple]++;
    bool test=BinomialStdvect_Increment(nbVert,3,eVect);
    if (!test)
      break;
  }
  return eInv;
}


// The matrix should be square and the output does not depends on the ordering of the coefficients.
template<bool is_symmetric, typename T>
inline typename std::enable_if<is_totally_ordered<T>::value,T>::type GetInvariantWeightMatrix(WeightMatrix<is_symmetric, T> const& WMat)
{
  static_assert(is_totally_ordered<T>::value, "Requires T to be totally ordered");
  size_t nbInv=3;
  size_t nbVert=WMat.rows();
  std::vector<int> ListM(nbInv);
  ListM[0]=2;
  ListM[1]=5;
  ListM[2]=23;
  size_t nbWeight=WMat.GetWeightSize();
  std::vector<T> ListWeight=WMat.GetWeight();
  std::vector<int> ListAtt(nbWeight,0);
  for (size_t iVert=0; iVert<nbVert; iVert++)
    for (size_t jVert=0; jVert<nbVert; jVert++) {
      Tidx_value iWeight=WMat.GetValue(iVert, jVert);
      ListAtt[iWeight]++;
    }
  std::vector<int> eList = SortingPerm<T,int>(ListWeight);
#ifdef DEBUG
  for (size_t jWeight=1; jWeight<nbWeight; jWeight++) {
    size_t iWeight=jWeight-1;
    int i=eList[iWeight];
    int j=eList[jWeight];
    if (ListWeight[i] > ListWeight[j]) {
      std::cerr << "Logical error in the comparison\n";
      throw TerminalException{1};
    }
  }
#endif
  T eMainInv=0;
  for (size_t iInv=0; iInv<nbInv; iInv++) {
    T eInv=0;
    T ePow=ListM[iInv];
    for (size_t iWeight=0; iWeight<nbWeight; iWeight++) {
      size_t jWeight=eList[iWeight];
      T prov2=ListAtt[jWeight]*ListWeight[jWeight];
      eInv *= ePow;
      eInv += prov2;
    }
    eMainInv += eInv;
  }
  return eMainInv;
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
   a smaller number of vertices gives a running time that is larger.
   Thus this idea while good looking is actually not a good one.
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





template<typename T>
size_t get_total_number_vertices(WeightMatrix<true, T> const& WMat)
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


template<typename T, typename Fcolor, typename Fadj>
void GetGraphFromWeightedMatrix_color_adj(WeightMatrix<true, T> const& WMat, Fcolor f_color, Fadj f_adj)
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
      int eVal;
      if (jVert == nbRow+1) {
	if (iVert == nbRow)
	  eVal=nbWei;
	else
	  eVal=nbWei+1;
      }
      else {
	if (jVert == nbRow)
	  eVal=WMat.GetValue(iVert, iVert);
	else
	  eVal=WMat.GetValue(iVert, jVert);
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


template<typename T>
size_t get_total_number_vertices(WeightMatrix<true, T> const& WMat)
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


template<typename T, typename Fcolor, typename Fadj>
void GetGraphFromWeightedMatrix_color_adj(WeightMatrix<true, T> const& WMat, Fcolor f_color, Fadj f_adj)
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
      int eVal;
      if (jVert == nbRow+1) {
	if (iVert == nbRow)
	  eVal=nbWei;
	else
	  eVal=nbWei+1;
      }
      else {
	if (jVert == nbRow)
	  eVal=WMat.GetValue(iVert, iVert);
	else
	  eVal=WMat.GetValue(iVert, jVert);
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


template<typename T>
bliss::Graph GetBlissGraphFromWeightedMatrix(WeightMatrix<true, T> const& WMat)
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



template<typename T, typename Tgr>
inline typename std::enable_if<(not is_functional_graph_class<Tgr>::value),Tgr>::type GetGraphFromWeightedMatrix(WeightMatrix<true, T> const& WMat)
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




template<typename T, typename Tgr>
inline typename std::enable_if<is_functional_graph_class<Tgr>::value,Tgr>::type GetGraphFromWeightedMatrix(WeightMatrix<true, T> const& WMat)
{
  size_t nbMult=WMat.GetWeightSize()+2;
  size_t hS=GetNeededPower(nbMult);
  size_t nbRow=WMat.rows();
  size_t nbVert=nbRow + 2;
  size_t nof_vertices=hS*nbVert;
  std::function<bool(int,int)> fAdj=[=](size_t const& aVert, size_t const& bVert) -> bool {
    int eVal;
    size_t iVert=aVert % nbVert;
    size_t iH=(aVert - iVert)/nbVert;
    size_t jVert=bVert % nbVert;
    size_t jH=(aVert - iVert)/nbVert;
    std::vector<int> eVect(hS);
    if (iVert == jVert) {
      if (iH != jH) {
	return true;
      }
      else {
	return false;
      }
    }
    if (iH == jH) {
      if (iVert > jVert) {
        std::swap(iVert, jVert);
      }
      if (jVert == nbRow+1) {
	if (iVert == nbRow)
	  eVal=nbMult;
	else
	  eVal=nbMult+1;
      }
      else {
	if (jVert == nbRow)
	  eVal=WMat.GetValue(iVert, iVert);
	else
	  eVal=WMat.GetValue(iVert, jVert);
      }
      GetBinaryExpression(eVal, hS, eVect);
      if (eVect[iH] == 1) {
	return true;
      }
      else {
	return false;
      }
    }
    return false;
  };
  std::function<int(int)> fColor=[=](int const& aVert) -> int {
    size_t iVert=aVert % nbVert;
    size_t iH=(aVert - iVert)/nbVert;
    return iH;
  };
  Tgr eGR(nof_vertices, fAdj);
  eGR.SetFColor(fColor);
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


template<typename Tgroup>
Tgroup GetGroupListGen(std::vector<std::vector<unsigned int>> const& ListGen, size_t const& nbVert)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> generatorList;
  std::vector<Tidx> gList(nbVert);
  for (auto & eGen : ListGen) {
    for (size_t iVert=0; iVert<nbVert; iVert++)
      gList[iVert]=eGen[iVert];
    generatorList.push_back(Telt(gList));
  }
  return Tgroup(generatorList, nbVert);
}


template<typename Tgr>
bool CheckListGenerators(std::vector<std::vector<unsigned int>> const& ListGen, Tgr const& eGR)
{
  size_t nbVert = eGR.GetNbVert();
  for (auto & eGen : ListGen) {
    for (size_t iVert=0; iVert<nbVert; iVert++) {
      int eColor = eGR.GetColor(iVert);
      int iVert_img = eGen[iVert];
      int fColor = eGR.GetColor(iVert_img);
      if (eColor != fColor) {
        return false;
      }
      //
      for (auto & jVert : eGR.Adjacency(iVert)) {
        int jVert_img = eGen[jVert];
        bool test = eGR.IsAdjacent(iVert_img, jVert_img);
        if (!test)
          return false;
      }
    }
  }
  return true;
}


template<typename Tgr, typename Tgroup>
void PrintStabilizerGroupSizes(std::ostream& os, Tgr const& eGR)
{
  bliss::Graph g=GetBlissGraphFromGraph(eGR);
  size_t nbVert=eGR.GetNbVert();
  std::vector<std::vector<unsigned int>> ListGen1 = BLISS_GetListGenerators(eGR);
  std::vector<std::vector<unsigned int>> ListGen2 = TRACES_GetListGenerators(eGR);
  auto siz1 = GetGroupListGen<Tgroup>(ListGen1, nbVert).size();
  auto siz2 = GetGroupListGen<Tgroup>(ListGen2, nbVert).size();
  bool test1 = CheckListGenerators(ListGen1, eGR);
  bool test2 = CheckListGenerators(ListGen2, eGR);
  os << "|GRP bliss|=" << siz1 << " |GRP traces|=" << siz2 << " test1=" << test1 << " test2=" << test2 << "\n";
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
template<typename T, typename Tgr, typename Tidx>
std::pair<std::vector<Tidx>, std::vector<Tidx>> GetCanonicalizationVector(WeightMatrix<true, T> const& WMat)
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
template<typename T, typename Tgr, typename Tidx>
std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> GetGroupCanonicalizationVector(WeightMatrix<true, T> const& WMat)
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
        Tidx pos1 = WMat.GetValue(i, j);
        Tidx pos2 = WMat.GetValue(iImg, jImg);
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



template<typename T>
void SignRenormalizationMatrix(MyMatrix<T> & M)
{
  size_t nbRow = M.rows();
  size_t n=M.cols();
  auto get_need_chgsign=[&](int const& iRow) -> bool {
    for (size_t i=0; i<n; i++) {
      T eVal = M(iRow,i);
      if (eVal != 0) {
        return eVal < 0;
      }
    }
    return false;
  };
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    if (get_need_chgsign(iRow)) {
      for (size_t i=0; i<n; i++)
        M(iRow,i) = - M(iRow,i);
    }
  }
}

template<typename T>
MyMatrix<T> ExpandReducedMatrix(MyMatrix<T> const& M)
{
  size_t nbPair=M.rows();
  size_t n=M.cols();
  MyMatrix<T> Mret(2*nbPair, n);
  for (size_t iPair=0; iPair<nbPair; iPair++)
    for (size_t i=0; i<n; i++) {
      Mret(2*iPair  , i) =  M(iPair, i);
      Mret(2*iPair+1, i) = -M(iPair, i);
    }
  return Mret;
}


/*
  Consider the case of the A2 root system with vectors
  \pm (1,0), \pm (0,1), \pm (1,1).
  If we consider the automorphisms of this vector configuration what we get is:
  ---Rotation by 2pi / 6 : Define subgroup of size 6
  ---Symmetry by axis.
  All together the group is of size 12.
  ----
  If we consider the absolute graph formed by the 3 vectors: (1,0), (0,1) and (1,1)
  then we get that this system defined a complete graph on 3 elements. So the group
  is of size 6. So, we indeed have the equality G = {\pm Id} x G_{abs}.
  ---
  The following holds:
  ---The construction of the weight matrix and so on means that orthogonal
  transformation on the vectors are not a problem.
  ---Since the absolute graph is the complete graph, we obtain that any ordering
  of the vector is possible by the canonicalization code.
  ---Thus if we put the vectors (1,0), (0,1) and (1,1)
  then the absolute canonicalization code may return us
  {(1,0), (0,1), (1,1)} or {(1,0), (1,1), (0,1)}.
  I think the hermite normal form of those are different.
  So, the method does not work.
  ---But we may be able to do something better. We still have the signs
  to be assigned.
*/
template<typename Tint>
EquivTest<MyMatrix<Tint>> LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick(MyMatrix<Tint> const& EXT, MyMatrix<Tint> const& Qmat)
{
  size_t nbRow= EXT.rows();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif

  WeightMatrixAbs<Tint> WMatAbs = GetSimpleWeightMatrixAntipodal_AbsTrick(EXT, Qmat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  //  std::cerr << "WMatAbs.positionZero=" << WMatAbs.positionZero << "\n";
  //  std::cerr << "WMatAbs.WMat=\n";
  //  PrintWeightedMatrix(std::cerr, WMatAbs.WMat);

  GraphBitset eGR=GetGraphFromWeightedMatrix<Tint,GraphBitset>(WMatAbs.WMat);
  //  GRAPH_PrintOutputGAP_vertex_colored("GAP_graph", eGR);

  //  std::cerr << "WMatAbs.WMat : ";
  //  PrintStabilizerGroupSizes(std::cerr, eGR);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|GetGraphFromWeightedMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif

#ifdef USE_BLISS
  std::pair<std::vector<unsigned int>, std::vector<std::vector<unsigned int>>> ePair = BLISS_GetCanonicalOrdering_ListGenerators(eGR);
#endif
#ifdef USE_TRACES
  std::pair<std::vector<unsigned int>, std::vector<std::vector<unsigned int>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators(eGR);
#endif
#ifdef DEBUG
  //  PrintStabilizerGroupSizes(std::cerr, eGR);
  std::string eExpr = GetCanonicalForm_string(eGR, ePair.first);
  mpz_class eHash1 = MD5_hash_mpz(eExpr);
  std::cerr << "eHash1=" << eHash1 << "\n";
  //
  size_t hS = eGR.GetNbVert() / (nbRow + 2);
  std::cerr << "|eGR|=" << eGR.GetNbVert() << " nbRow=" << nbRow << " hS=" << hS << "\n";
  for (auto & eGen : ePair.second) {
    std::vector<unsigned int> eGenRed(nbRow+2);
    for (size_t i=0; i<nbRow+2; i++) {
      unsigned int val = eGen[i];
      if (val >= nbRow+2) {
        std::cerr << "At i=" << i << " we have val=" << val << " nbRow=" << nbRow << "\n";
        throw TerminalException{1};
      }
      eGenRed[i] = val;
    }
    for (size_t i=nbRow; i<nbRow+2; i++) {
      if (eGenRed[i] != i) {
        std::cerr << "Point is not preserved\n";
        throw TerminalException{1};
      }
    }
    for (size_t iH=0; iH<hS; iH++) {
      for (size_t i=0; i<nbRow+2; i++) {
        unsigned int val1 = eGen[i + iH * (nbRow+2)];
        unsigned int val2 = iH * (nbRow+2) + eGenRed[i];
        if (val1 != val2) {
          std::cerr << "val1=" << val1 << " val2=" << val2 << "\n";
          std::cerr << "iH" << iH << " i=" << i << " hS=" << hS << " nbRow=" << nbRow << "\n";
          throw TerminalException{1};
        }
      }
      for (size_t i=0; i<nbRow; i++) {
        for (size_t j=0; j<nbRow; j++) {
          int iImg = eGenRed[i];
          int jImg = eGenRed[j];
          Tidx_value pos1 = WMatAbs.WMat.GetValue(i, j);
          Tidx_value pos2 = WMatAbs.WMat.GetValue(iImg, jImg);
          if (pos1 != pos2) {
            std::cerr << "Inconsistency at i=" << i << " j=" <<j << "\n";
            std::cerr << "iImg=" << iImg << " jImg=" << jImg << "\n";
            throw TerminalException{1};
          }
        }
      }
    }

  }
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|GetListGenerators|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif

  // We check if the Generating vector eGen can be mapped from the absolute
  // graph to the original one.
  auto TestExistSignVector=[&](std::vector<unsigned int> const& eGen) -> bool {
    /* We map a vector v_i to another v_j with sign +-1
       V[i] = 0 for unassigned
              1 for positive sign
              2 for negative sign
              3 for positive sign and treated
              4 for negative sign and treated
     */
    std::vector<uint8_t> V(nbRow, 0);
    V[0] = 1;
    while(true) {
      bool IsFinished = true;
      for (size_t i=0; i<nbRow; i++) {
        uint8_t val = V[i];
        if (val < 3 && val != 0) {
          IsFinished=false;
          V[i] = val + 2;
          size_t iImg = eGen[i];
          for (size_t j=0; j<nbRow; j++) {
            size_t jImg = eGen[j];
            Tidx_value pos = WMatAbs.WMat.GetValue(i, j);
            if (pos != WMatAbs.positionZero) {
              bool ChgSign1 = WMatAbs.ArrSigns[i + nbRow * j];
              bool ChgSign2 = WMatAbs.ArrSigns[iImg + nbRow * jImg];
              bool ChgSign = ChgSign1 ^ ChgSign2; // true if ChgSign1 != ChgSign2
              uint8_t valJ;
              if ((ChgSign && val == 1) || (!ChgSign && val == 2))
                valJ = 2;
              else
                valJ = 1;
              if (V[j] == 0) {
                V[j] = valJ;
              } else {
                if ((valJ % 2) != (V[j] % 2)) {
                  return false;
                }
              }
            }
          }
        }
      }
      if (IsFinished)
        break;
    }
    return true;
  };
  auto IsCorrectListGen=[&]() -> bool {
    for (auto& eGen : ePair.second) {
      bool test = TestExistSignVector(eGen);
      //      std::cerr << "test=" << test << "\n";
      if (!test)
        return false;
    }
    return true;
  };
  if (!IsCorrectListGen()) {
    return {false, {}};
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time5 = std::chrono::system_clock::now();
  std::cerr << "|Check Generators|=" << std::chrono::duration_cast<std::chrono::microseconds>(time5 - time4).count() << "\n";
#endif
  //
  std::pair<std::vector<int>, std::vector<int>> PairCanonic = GetCanonicalizationVector_KernelBis<int>(nbRow, ePair.first);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time6 = std::chrono::system_clock::now();
  std::cerr << "|GetCanonicalizationVector_Kernel|=" << std::chrono::duration_cast<std::chrono::microseconds>(time6 - time4).count() << "\n";
#endif

  size_t n_cols=EXT.cols();
  MyMatrix<Tint> EXTreord(nbRow, n_cols);
  std::vector<int> ListSigns(nbRow,0);
  ListSigns[0]=1;
#ifdef DEBUG
  std::string strAssign;
  std::cerr << "positionZero=" << WMatAbs.positionZero << "\n";
#endif
  auto SetSign=[&](size_t const& i_row) -> void {
    int i_row_orig = PairCanonic.second[i_row];
    for (size_t k_row=0; k_row<nbRow; k_row++) {
      if (k_row != i_row && ListSigns[k_row] != 0) {
        int k_row_orig = PairCanonic.second[k_row];
        if (WMatAbs.WMat.GetValue(i_row_orig, k_row_orig) != WMatAbs.positionZero) {
          bool ChgSign = WMatAbs.ArrSigns[i_row_orig + nbRow * k_row_orig];
          int ValSign = 1 - 2*int(ChgSign);
          int RetSign = ValSign * ListSigns[k_row];
          ListSigns[i_row] = RetSign;
#ifdef DEBUG
          strAssign += " (" + std::to_string(i_row) + " / " + std::to_string(k_row) + ")";
#endif
          return;
        }
      }
    }
  };
  while(true) {
    int nbUndone=0;
    for (size_t i_row=0; i_row<nbRow; i_row++)
      if (ListSigns[i_row] == 0) {
        nbUndone++;
        SetSign(i_row);
      }
    if (nbUndone == 0)
      break;
  };
#ifdef DEBUG
  mpz_class eHash2 = MD5_hash_mpz(strAssign);
  std::cerr << "strAssign=" << strAssign << "\n";
  std::cerr << "eHash2=" << eHash2 << "\n";
#endif
#ifdef DEBUG
  std::string strWMat;
  for (size_t i_row=0; i_row<nbRow; i_row++) {
    int i_rowC = PairCanonic.second[i_row];
    for (size_t j_row=0; j_row<nbRow; j_row++) {
      int j_rowC = PairCanonic.second[j_row];
      Tidx_value pos = WMatAbs.WMat.GetValue(i_rowC, j_rowC);
      strWMat += " " + std::to_string(pos);
    }
  }
  for (auto & eVal : WMatAbs.WMat.GetWeight()) {
    strWMat += " " + std::to_string(eVal);
  }
  mpz_class eHash3 = MD5_hash_mpz(strWMat);
  std::cerr << "eHash3=" << eHash3 << "\n";
#endif

  for (size_t i_row=0; i_row<nbRow; i_row++) {
    int j_row = PairCanonic.second[i_row];
    int eSign = ListSigns[i_row];
    for (size_t i_col=0; i_col<n_cols; i_col++)
      EXTreord(i_row, i_col) = eSign * EXT(j_row, i_col);
  }
#ifdef DEBUG
  std::cerr << "EXTreord=\n";
  WriteMatrix(std::cerr, EXTreord);
  WriteMatrixGAP(std::cerr, EXTreord);
  std::cerr << "\n";
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time7 = std::chrono::system_clock::now();
  std::cerr << "|EXTreord|=" << std::chrono::duration_cast<std::chrono::microseconds>(time7 - time6).count() << "\n";
#endif

  MyMatrix<Tint> RedMat = ComputeColHermiteNormalForm(EXTreord).second;
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time8 = std::chrono::system_clock::now();
  std::cerr << "|ComputeColHermiteNormalForm|=" << std::chrono::duration_cast<std::chrono::microseconds>(time8 - time7).count() << "\n";
#endif
  SignRenormalizationMatrix(RedMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time9 = std::chrono::system_clock::now();
  std::cerr << "|SignRenormalizationMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time9 - time8).count() << "\n";
#endif

  return {true, std::move(RedMat)};
}





template<typename Tint>
EquivTest<std::vector<std::vector<unsigned int>>> LinPolytopeAntipodalIntegral_Automorphism_AbsTrick(MyMatrix<Tint> const& EXT, MyMatrix<Tint> const& Qmat)
{
  size_t nbRow= EXT.rows();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif

  WeightMatrixAbs<Tint> WMatAbs = GetSimpleWeightMatrixAntipodal_AbsTrick(EXT, Qmat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  //  std::cerr << "WMatAbs.positionZero=" << WMatAbs.positionZero << "\n";
  //  std::cerr << "WMatAbs.WMat=\n";
  //  PrintWeightedMatrix(std::cerr, WMatAbs.WMat);

  GraphBitset eGR=GetGraphFromWeightedMatrix<Tint,GraphBitset>(WMatAbs.WMat);
  //  GRAPH_PrintOutputGAP_vertex_colored("GAP_graph", eGR);

  //  std::cerr << "WMatAbs.WMat : ";
  //  PrintStabilizerGroupSizes(std::cerr, eGR);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|GetGraphFromWeightedMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif

#ifdef USE_BLISS
  std::vector<std::vector<unsigned int>> ListGen = BLISS_GetListGenerators(eGR);
#endif
#ifdef USE_TRACES
  std::vector<std::vector<unsigned int>> ListGen = TRACES_GetListGenerators(eGR);
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|GetListGenerators|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif

  // We check if the Generating vector eGen can be mapped from the absolute
  // graph to the original one.
  std::vector<std::vector<unsigned int>> ListGenRet;
  auto TestExistSignVector=[&](std::vector<unsigned int> const& eGen) -> bool {
    /* We map a vector v_i to another v_j with sign +-1
       V[i] = 0 for unassigned
              1 for positive sign
              2 for negative sign
              3 for positive sign and treated
              4 for negative sign and treated
     */
    std::vector<uint8_t> V(nbRow, 0);
    std::vector<unsigned int> eGenRet(2*nbRow,0);
    auto setSign=[&](int const& idx, uint8_t const& val) -> void {
      if (val == 1) {
        eGenRet[idx        ] = eGen[idx];
        eGenRet[idx + nbRow] = eGen[idx] + nbRow;
      } else {
        eGenRet[idx        ] = eGen[idx] + nbRow;
        eGenRet[idx + nbRow] = eGen[idx];
      }
      V[idx] = val;
    };
    setSign(0, 1);
    while(true) {
      bool IsFinished = true;
      for (size_t i=0; i<nbRow; i++) {
        uint8_t val = V[i];
        if (val < 3 && val != 0) {
          IsFinished=false;
          V[i] = val + 2;
          size_t iImg = eGen[i];
          for (size_t j=0; j<nbRow; j++) {
            size_t jImg = eGen[j];
            Tidx_value pos = WMatAbs.WMat.GetValue(i, j);
            if (pos != WMatAbs.positionZero) {
              bool ChgSign1 = WMatAbs.ArrSigns[i + nbRow * j];
              bool ChgSign2 = WMatAbs.ArrSigns[iImg + nbRow * jImg];
              bool ChgSign = ChgSign1 ^ ChgSign2; // true if ChgSign1 != ChgSign2
              uint8_t valJ;
              if ((ChgSign && val == 1) || (!ChgSign && val == 2))
                valJ = 2;
              else
                valJ = 1;
              if (V[j] == 0) {
                setSign(j, valJ);
              } else {
                if ((valJ % 2) != (V[j] % 2)) {
                  return false;
                }
              }
            }
          }
        }
      }
      if (IsFinished)
        break;
    }
    ListGenRet.push_back(eGenRet);
    return true;
  };
  auto IsCorrectListGen=[&]() -> bool {
    for (auto& eGen : ListGen) {
      bool test = TestExistSignVector(eGen);
      //      std::cerr << "test=" << test << "\n";
      if (!test)
        return false;
    }
    return true;
  };
  if (!IsCorrectListGen()) {
    return {false, {}};
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time5 = std::chrono::system_clock::now();
  std::cerr << "|Check Generators|=" << std::chrono::duration_cast<std::chrono::microseconds>(time5 - time4).count() << "\n";
#endif
  //
  std::vector<unsigned int> AntipodalGen(2*nbRow,0);
  for (size_t iRow=0; iRow<nbRow; iRow++) {
    AntipodalGen[iRow] = iRow + nbRow;
    AntipodalGen[nbRow + iRow] = iRow;
  }
  ListGenRet.push_back(AntipodalGen);
  //
  return {true, std::move(ListGenRet)};
}











template<typename T, typename Tgroup>
Tgroup GetStabilizerWeightMatrix(WeightMatrix<true, T> const& WMat)
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


template<typename T, typename Tgroup>
Tgroup GetStabilizerAsymmetricMatrix(WeightMatrix<false, T> const& WMatI)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  WeightMatrix<true, T> WMatO=GetSymmetricWeightMatrix(WMatI);
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





template<typename T>
EquivTest<std::vector<unsigned int>> TestEquivalenceWeightMatrix_norenorm(WeightMatrix<true, T> const& WMat1, WeightMatrix<true, T> const& WMat2)
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
#ifdef DEBUG
  for (unsigned int iVert1=0; iVert1<nbRow; iVert1++) {
    unsigned int iVert2=TheEquiv[iVert1];
    for (unsigned int jVert1=0; jVert1<nbRow; jVert1++) {
      unsigned int jVert2=TheEquiv[jVert1];
      Tidx_value eVal1=WMat1.GetValue(iVert1, jVert1);
      Tidx_value eVal2=WMat2.GetValue(iVert2, jVert2);
      if (eVal1 != eVal2) {
	std::cerr << "nbRow=" << nbRow << "\n";
	std::cerr << "nof_vertices=" << nof_vertices << "\n";
	std::cerr << "iVert1=" << iVert1 << " jVert1=" << jVert1 << "\n";
	std::cerr << "iVert2=" << iVert2 << " jVert2=" << jVert2 << "\n";
	std::cerr << "eVal1=" << eVal1 << " eVal2=" << eVal2 << "\n";
	std::cerr << "Our reduction technique is broken\n";
	std::cerr << "Please panic and debug\n";
	throw TerminalException{1};
      }
    }
  }
#endif
  return {true, std::move(TheEquiv)};
}

template<typename T, typename Telt>
EquivTest<Telt> TestEquivalenceWeightMatrix_norenorm_perm(WeightMatrix<true, T> const& WMat1, WeightMatrix<true, T> const& WMat2)
{
  EquivTest<std::vector<unsigned int>> ePair = TestEquivalenceWeightMatrix_norenorm(WMat1, WMat2);
  size_t len = ePair.TheEquiv.size();
  std::vector<int> eList(len);
  for (size_t i=0; i<len; i++)
    eList[i] = ePair.TheEquiv[i];
  Telt ePerm(eList);
  return {ePair.TheReply, std::move(ePerm)};
}

template<typename T, typename Telt>
EquivTest<Telt> TestEquivalenceWeightMatrix(WeightMatrix<true, T> const& WMat1, WeightMatrix<true, T> &WMat2)
{
  bool test=RenormalizeWeightMatrix(WMat1, WMat2);
  if (!test)
    return {false, {}};
  return TestEquivalenceWeightMatrix_norenorm_perm<T,Telt>(WMat1, WMat2);
}





template<typename T, typename Telt>
EquivTest<Telt> GetEquivalenceAsymmetricMatrix(WeightMatrix<false, T> const& WMat1, WeightMatrix<false,T> const& WMat2)
{
  using Tidx = typename Telt::Tidx;
  WeightMatrix<true, T> WMatO1=GetSymmetricWeightMatrix<T>(WMat1);
  WeightMatrix<true, T> WMatO2=GetSymmetricWeightMatrix<T>(WMat2);
  EquivTest<Telt> eResEquiv=TestEquivalenceWeightMatrix<T,Telt>(WMatO1, WMatO2);
  if (!eResEquiv.TheReply)
    return eResEquiv;
  size_t nbSHV=WMat1.rows();
  std::vector<Tidx> v(nbSHV);
  for (size_t i=0; i<nbSHV; i++)
    v[i]=eResEquiv.TheEquiv.at(i);
  return {true, std::move(Telt(v))};
}





template<typename T, typename Telt>
EquivTest<Telt> TestEquivalenceSubset(WeightMatrix<true, T> const& WMat, Face const& f1, Face const& f2)
{
  using Tidx = typename Telt::Tidx;
  size_t siz=WMat.GetWeightSize();
  size_t n=WMat.rows();
  auto g=[&](Face const& f, size_t iRow, size_t iCol) -> int {
     if (iRow < n && iCol < n)
       return WMat.GetValue(iRow,iCol);
     if (iRow == n && iCol == n)
       return siz + 2;
     if (iRow == n) {
       if (f[iCol] == 0)
         return siz;
       else
         return siz + 1;
     }
     // Last case: Necessarily we have iCol == n && iRow < n
     if (f[iRow] == 0)
       return siz;
     else
       return siz + 1;
  };
  WeightMatrix<true,int> WMat1(n+1,[&](size_t iRow, size_t iCol) -> int {
    return g(f1, iRow, iCol);
  });
  WeightMatrix<true,int> WMat2(n+1,[&](size_t iRow, size_t iCol) -> int {
    return g(f2, iRow, iCol);
  });
  EquivTest<Telt> test=TestEquivalenceWeightMatrix_norenorm_perm<int,Telt>(WMat1, WMat2);
  if (!test.TheReply)
    return {false, {}};
  std::vector<Tidx> eList(n);
  for (size_t i=0; i<n; i++) {
    int eVal=test.TheEquiv.at(i);
    eList[i] = eVal;
  }
  return {true, std::move(Telt(eList))};
}



template<typename T, typename Tgroup>
Tgroup StabilizerSubset(WeightMatrix<true,T> const& WMat, Face const& f)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  size_t siz=WMat.GetWeightSize();
  size_t n=WMat.rows();
  auto g=[&](size_t iRow, size_t iCol) -> int {
     if (iRow < n && iCol < n)
       return WMat.GetValue(iRow,iCol);
     if (iRow == n && iCol == n)
       return siz + 2;
     if (iRow == n) {
       if (f[iCol] == 0)
         return siz;
       else
         return siz + 1;
     }
     // Last case: Necessarily we have iCol == n && iRow < n
     if (f[iRow] == 0)
       return siz;
     else
       return siz + 1;
  };
  WeightMatrix<true,int> WMatW(n+1, g);
  Tgroup GRP=GetStabilizerWeightMatrix<T,Tgroup>(WMatW);
  std::vector<Telt> ListPerm;
  for (auto & ePerm : GRP.GeneratorsOfGroup()) {
    std::vector<Tidx> eList(n);
    for (size_t i=0; i<n; i++)
      eList[i]=OnPoints(i, ePerm);
    ListPerm.push_back(Telt(eList));
  }
  return Tgroup(ListPerm, n);
}



template<bool is_symmetric,typename T>
WeightMatrix<is_symmetric,int> NakedWeightedMatrix(WeightMatrix<is_symmetric,T> const& WMat)
{
  size_t n=WMat.rows();
  WeightMatrix<is_symmetric, int> WMatNaked(n);
  for (size_t i=0; i<n; i++)
    for (size_t j=0; j<n; j++) {
      Tidx_value eVal=WMat.GetValue(i,j);
      WMatNaked.intDirectAssign(i,j,eVal);
    }
  size_t nbWeight=WMat.GetWeightSize();
  std::vector<int> ListWeight(nbWeight);
  for (size_t i=0; i<nbWeight; i++)
    ListWeight[i]=i;
  WMatNaked.SetWeight(ListWeight);
  return WMatNaked;
}


template<typename T, typename Telt>
Telt CanonicalizeWeightMatrix(WeightMatrix<true, T> const& WMat)
{
  using Tidx = typename Telt::Tidx;
  GraphBitset eGR=GetGraphFromWeightedMatrix<T,GraphBitset>(WMat);
  bliss::Graph g=GetBlissGraphFromGraph(eGR);
  unsigned int nof_vertices=eGR.GetNbVert();
  bliss::Stats stats;
  const unsigned int* cl=g.canonical_form(stats, &report_aut_void, stderr);
  std::vector<int> clR(nof_vertices);
  for (unsigned int i=0; i<nof_vertices; i++)
    clR[cl[i]]=i;
  unsigned int nbVert=WMat.GetNbVert();
  std::vector<unsigned int> eVect_R(nof_vertices);
  std::vector<Tidx> eVect(nbVert);
  for (unsigned int iVert=0; iVert<nbVert; iVert++) {
    unsigned int i=cl[iVert]; // or clR? this needs to be debugged
    eVect_R[i]=iVert;
    eVect[iVert] = i;
  }
  // Need to check that code.
  return Telt(eVect);
}



template<typename T, typename Tgroup>
Tgroup LinPolytope_Automorphism(MyMatrix<T> const & EXT)
{
  MyMatrix<T> EXTred=ColumnReduction(EXT);
  WeightMatrix<true,T> WMat=GetWeightMatrix(EXTred);
  return GetStabilizerWeightMatrix<T,Tgroup>(WMat);
}


template<typename T>
size_t GetInvariant_ListMat_Subset(MyMatrix<T> const& EXT, std::vector<MyMatrix<T>> const&ListMat, Face const& eSubset)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif


  WeightMatrix<true,std::vector<T>> WMat = GetWeightMatrix_ListMat_Subset(EXT, ListMat, eSubset);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetWeightMatrix_ListMatrix_Subset|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif


  WMat.ReorderingSetWeight();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|ReorderingSetWeight|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif


  size_t e_hash = std::hash<WeightMatrix<true,std::vector<T>>>()(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|hash|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif
  return e_hash;
}

template<typename T>
std::vector<std::vector<unsigned int>> GetListGenAutomorphism_ListMat_Subset(MyMatrix<T> const& EXT, std::vector<MyMatrix<T>> const&ListMat, Face const& eSubset)
{
  size_t nbRow = EXT.rows();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif


  WeightMatrix<true,std::vector<T>> WMat = GetWeightMatrix_ListMat_Subset(EXT, ListMat, eSubset);
  // No need to reorder in autom case
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetWeightMatrix_ListMatrix_Subset|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif

  //  using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
  Tgr eGR=GetGraphFromWeightedMatrix<std::vector<T>,Tgr>(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|GetGraphFromWeightMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif


#ifdef USE_BLISS
  std::vector<std::vector<unsigned int>> ListGenTot = BLISS_GetListGenerators(eGR);
#endif
#ifdef USE_TRACES
  std::vector<std::vector<unsigned int>> ListGenTot = TRACES_GetListGenerators(eGR);
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|GetListGenerators|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif


  std::vector<std::vector<unsigned int>> ListGen;
  for (auto & eGen : ListGenTot) {
    std::vector<unsigned int> eGenRed(nbRow);
    for (size_t i=0; i<nbRow; i++) {
      unsigned int val = eGen[i];
#ifdef DEBUG
      unsigned int nbRow_ui = nbRow;
      if (val >= nbRow_ui) {
        std::cerr << "At i=" << i << " we have val=" << val << " nbRow=" << nbRow << "\n";
        throw TerminalException{1};
      }
#endif
      eGenRed[i] = val;
    }
    ListGen.push_back(eGenRed);
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time5 = std::chrono::system_clock::now();
  std::cerr << "|ListGen|=" << std::chrono::duration_cast<std::chrono::microseconds>(time5 - time4).count() << "\n";
#endif
  return ListGen;
}






template<typename T>
EquivTest<std::vector<unsigned int>> TestEquivalence_ListMat_Subset(MyMatrix<T> const& EXT1, std::vector<MyMatrix<T>> const&ListMat1, Face const& eSubset1,
                                                                    MyMatrix<T> const& EXT2, std::vector<MyMatrix<T>> const&ListMat2, Face const& eSubset2)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif

  WeightMatrix<true,std::vector<T>> WMat1 = GetWeightMatrix_ListMat_Subset(EXT1, ListMat1, eSubset1);
  WeightMatrix<true,std::vector<T>> WMat2 = GetWeightMatrix_ListMat_Subset(EXT2, ListMat2, eSubset2);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetWeightMatrix_ListMatrix_Subset|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif


  WMat1.ReorderingSetWeight();
  WMat2.ReorderingSetWeight();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|ReorderingSetWeight|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif


  EquivTest<std::vector<unsigned int>> PairTest = TestEquivalenceWeightMatrix_norenorm(WMat1, WMat2);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|TestEquivalence_ListMat_Subset|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif
  return PairTest;
}







template<typename Tint>
MyMatrix<Tint> LinPolytopeIntegral_CanonicForm(MyMatrix<Tint> const& EXT)
{
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif

  WeightMatrix<true,Tint> WMat=GetWeightMatrix(EXT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
#endif

  WMat.ReorderingSetWeight();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
#endif

  std::pair<std::vector<int>, std::vector<int>> PairCanonic = GetCanonicalizationVector<Tint,GraphBitset,int>(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
#endif

  MyMatrix<Tint> EXTreord(n_rows, n_cols);
  for (size_t i_row=0; i_row<n_rows; i_row++) {
    size_t j_row = PairCanonic.second[i_row];
    for (size_t i_col=0; i_col<n_cols; i_col++)
      EXTreord(i_row, i_col) = EXT(j_row, i_col);
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time5 = std::chrono::system_clock::now();
#endif

  MyMatrix<Tint> RedMat = ComputeColHermiteNormalForm(EXTreord).second;
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time6 = std::chrono::system_clock::now();
  std::cerr << "|GetWeightMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
  std::cerr << "|ReorderingSetWeight|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
  std::cerr << "|GetCanonicalizationVector|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
  std::cerr << "|EXTreord|=" << std::chrono::duration_cast<std::chrono::microseconds>(time5 - time4).count() << "\n";
  std::cerr << "|ComputeColHermiteNormalForm|=" << std::chrono::duration_cast<std::chrono::microseconds>(time6 - time5).count() << "\n";
#endif
  return RedMat;
}





template<typename Tint>
MyMatrix<Tint> LinPolytopeAntipodalIntegral_CanonicForm(MyMatrix<Tint> const& EXT)
{
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  MyMatrix<Tint> Qmat=GetQmatrix(EXT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetQmatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif

  EquivTest<MyMatrix<Tint>> EauivTest = LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick(EXT, Qmat);
  if (EauivTest.TheReply) {
    return EauivTest.TheEquiv;
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|LinPolytopeAntipodalIntegral_CanonicForm_AbsTrick|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif

  WeightMatrix<true,Tint> WMat=GetWeightMatrixAntipodal(EXT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|GetWeightMatrixAntipodal|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif
  //  std::cerr << "After direct construction WMat=\n";
  //  PrintWeightedMatrix(std::cerr, WMat);

  WMat.ReorderingSetWeight();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time5 = std::chrono::system_clock::now();
  std::cerr << "|ReorderingSetWeight|=" << std::chrono::duration_cast<std::chrono::microseconds>(time5 - time4).count() << "\n";
#endif

  std::pair<std::vector<int>, std::vector<int>> PairCanonic = GetCanonicalizationVector<Tint,GraphBitset,int>(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time6 = std::chrono::system_clock::now();
  std::cerr << "|GetCanonicalizationVector|=" << std::chrono::duration_cast<std::chrono::microseconds>(time6 - time5).count() << "\n";
#endif

  MyMatrix<Tint> EXTreord(n_rows, n_cols);
  size_t idx=0;
  Face IsIncluded(n_rows);
  for (size_t i_row=0; i_row<2*n_rows; i_row++) {
    int j_row = PairCanonic.second[i_row];
    int res = j_row % 2;
    int pos = j_row / 2;
    if (res == 0) {
      if (IsIncluded[pos] == 0) {
        IsIncluded[pos]=1;
        for (size_t i_col=0; i_col<n_cols; i_col++)
          EXTreord(idx, i_col) = EXT(pos, i_col);
        idx++;
      }
    } else {
      if (IsIncluded[pos] == 0) {
        IsIncluded[pos]=1;
        for (size_t i_col=0; i_col<n_cols; i_col++)
          EXTreord(idx, i_col) = -EXT(pos, i_col);
        idx++;
      }
    }
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time7 = std::chrono::system_clock::now();
  std::cerr << "|EXTreord 2|=" << std::chrono::duration_cast<std::chrono::microseconds>(time7 - time6).count() << "\n";
#endif

  MyMatrix<Tint> RedMat = ComputeColHermiteNormalForm(EXTreord).second;
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time8 = std::chrono::system_clock::now();
  std::cerr << "|ComputeColHermiteNormalForm 2|=" << std::chrono::duration_cast<std::chrono::microseconds>(time8 - time7).count() << "\n";
#endif

  SignRenormalizationMatrix(RedMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time9 = std::chrono::system_clock::now();
  std::cerr << "|SignRenormalizationMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time9 - time8).count() << "\n";
#endif
  return RedMat;
}




template<typename Tint>
std::vector<std::vector<unsigned int>> LinPolytopeAntipodalIntegral_Automorphism(MyMatrix<Tint> const& EXT)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  MyMatrix<Tint> Qmat=GetQmatrix(EXT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetQmatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif

  EquivTest<std::vector<std::vector<unsigned int>>> EquivTest = LinPolytopeAntipodalIntegral_Automorphism_AbsTrick(EXT, Qmat);
  if (EquivTest.TheReply) {
    return EquivTest.TheEquiv;
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|LinPolytopeAntipodalIntegral_Automorphism_AbsTrick|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif

  WeightMatrix<true,Tint> WMat=GetWeightMatrixAntipodal(EXT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|GetWeightMatrixAntipodal|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif

  GraphBitset eGR=GetGraphFromWeightedMatrix<Tint,GraphBitset>(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time5 = std::chrono::system_clock::now();
  std::cerr << "|GetGraphFromWeightedMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time5 - time4).count() << "\n";
#endif

#ifdef USE_BLISS
  std::vector<std::vector<unsigned int>> ListGen = BLISS_GetListGenerators(eGR);
#endif
#ifdef USE_TRACES
  std::vector<std::vector<unsigned int>> ListGen = TRACES_GetListGenerators(eGR);
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time6 = std::chrono::system_clock::now();
  std::cerr << "|*_GetListGenerators|=" << std::chrono::duration_cast<std::chrono::microseconds>(time6 - time5).count() << "\n";
#endif
  return ListGen;
}











template<typename T, typename Telt>
MyMatrix<T> RepresentVertexPermutation(MyMatrix<T> const& EXT1, MyMatrix<T> const& EXT2, Telt const& ePerm)
{
  SelectionRowCol<T> eSelect=TMat_SelectRowCol(EXT1);
  std::vector<int> ListRowSelect=eSelect.ListRowSelect;
  MyMatrix<T> M1=SelectRow(EXT1, ListRowSelect);
  MyMatrix<T> M1inv=Inverse(M1);
  size_t nbRow=ListRowSelect.size();
  std::vector<int> ListRowSelectImg(nbRow);
  for (size_t iRow=0; iRow<nbRow; iRow++)
    ListRowSelectImg[iRow]=ePerm.at(iRow);
  MyMatrix<T> M2=SelectRow(EXT2, ListRowSelectImg);
  return M1inv*M2;
}



template<typename T, typename Telt>
Telt GetPermutationOnVectors(MyMatrix<T> const& EXT1, MyMatrix<T> const& EXT2)
{
  using Tidx = typename Telt::Tidx;
  size_t nbVect=EXT1.rows();
  std::vector<MyVector<T>> EXTrow1(nbVect), EXTrow2(nbVect);
  for (size_t iVect=0; iVect<nbVect; iVect++) {
    EXTrow1[iVect]=GetMatrixRow(EXT1, iVect);
    EXTrow2[iVect]=GetMatrixRow(EXT2, iVect);
  }
  Telt ePerm1=Telt(SortingPerm<MyVector<T>,Tidx>(EXTrow1));
  Telt ePerm2=Telt(SortingPerm<MyVector<T>,Tidx>(EXTrow2));
  Telt ePermRet=(~ePerm1) * ePerm2;
#ifdef DEBUG
  for (size_t iVect=0; iVect<nbVect; iVect++) {
    size_t jVect=ePermRet.at(iVect);
    if (EXTrow2[jVect] != EXTrow1[iVect]) {
      std::cerr << "iVect=" << iVect << " jVect=" << jVect << "\n";
      std::cerr << "Id:";
      for (size_t k=0; k<nbVect; k++)
	std::cerr << " " << k;
      std::cerr << "\n";
      std::cerr << " p:";
      for (size_t k=0; k<nbVect; k++)
	std::cerr << " " << ePermRet.at(k);
      std::cerr << "\n";

      std::cerr << "perm1:";
      for (size_t k=0; k<nbVect; k++)
	std::cerr << " " << ePerm1.at(k);
      std::cerr << "\n";
      std::cerr << "perm2:";
      for (size_t k=0; k<nbVect; k++)
	std::cerr << " " << ePerm2.at(k);
      std::cerr << "\n";


      std::cerr << "EXTrow1[iVect]=";
      WriteVector(std::cerr, EXTrow1[iVect]);
      std::cerr << "EXTrow2[jVect]=";
      WriteVector(std::cerr, EXTrow2[jVect]);
      std::cerr << "EXTrow1=\n";
      WriteMatrix(std::cerr, EXT1);
      std::cerr << "EXTrow2=\n";
      WriteMatrix(std::cerr, EXT2);
      std::cerr << "Error in GetPermutationOnVectors\n";
      throw TerminalException{1};
    }
  }
#endif
  return ePermRet;
}





#endif
