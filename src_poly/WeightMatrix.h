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


namespace std {
  template <bool is_symmetric,typename T, typename Tidx_value>
  struct hash<WeightMatrix<is_symmetric,T,Tidx_value>>
  {
    std::size_t operator()(WeightMatrix<is_symmetric,T,Tidx_value> const& WMat) const
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



template<bool is_symmetric, typename T, typename Tidx_value>
void PrintWeightedMatrixGAP(std::ostream &os, WeightMatrix<is_symmetric,T,Tidx_value> const&WMat)
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










template<bool is_symmetric, typename T, typename Tout, typename Tidx_value>
std::vector<Tout> GetLocalInvariantWeightMatrix(WeightMatrix<is_symmetric,T,Tidx_value> const&WMat, Face const& eSet)
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
  using Tidx_value = typename WeightMatrix<true, T>::Tidx_value;
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



template<typename T, typename Tgroup, typename Tidx_value>
LocalInvInfo ComputeLocalInvariantStrategy(WeightMatrix<true,T,Tidx_value> const&WMat, Tgroup const& GRP, std::string const& strat, std::ostream & os)
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
  } else {
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
    } else {
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
template<bool is_symmetric, typename T, typename Tidx_value>
inline typename std::enable_if<is_totally_ordered<T>::value,T>::type GetInvariantWeightMatrix(WeightMatrix<is_symmetric, T, Tidx_value> const& WMat)
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
      int eVal;
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
      int eVal;
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




template<typename T, typename Tgr, typename Tidx_value>
inline typename std::enable_if<is_functional_graph_class<Tgr>::value,Tgr>::type GetGraphFromWeightedMatrix(WeightMatrix<true, T, Tidx_value> const& WMat)
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
      } else {
	return false;
      }
    }
    if (iH == jH) {
      if (iVert > jVert) {
        std::swap(iVert, jVert);
      }
      if (jVert == nbRow+1) {
	if (iVert == nbRow)
	  eVal = nbMult;
	else
	  eVal = nbMult+1;
      } else {
	if (jVert == nbRow)
	  eVal = WMat.GetValue(iVert, iVert);
	else
	  eVal = WMat.GetValue(iVert, jVert);
      }
      GetBinaryExpression(eVal, hS, eVect);
      if (eVect[iH] == 1) {
	return true;
      } else {
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





template<typename T, typename Telt, typename Tidx_value>
EquivTest<Telt> TestEquivalenceSubset(WeightMatrix<true, T, Tidx_value> const& WMat, Face const& f1, Face const& f2)
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



template<typename T, typename Tgroup, typename Tidx_value>
Tgroup StabilizerSubset(WeightMatrix<true, T, Tidx_value> const& WMat, Face const& f)
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



template<bool is_symmetric,typename T, typename Tidx_value>
WeightMatrix<is_symmetric,int,Tidx_value> NakedWeightedMatrix(WeightMatrix<is_symmetric,T,Tidx_value> const& WMat)
{
  size_t n=WMat.rows();
  auto f=[&](size_t iRow, size_t iCol) -> int {
    return WMat.GetValue(iRow, iCol);
  };
  return WeightMatrix<is_symmetric, int, Tidx_value>(n, f);
}


template<typename T, typename Telt, typename Tidx_value>
Telt CanonicalizeWeightMatrix(WeightMatrix<true, T, Tidx_value> const& WMat)
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



#endif
