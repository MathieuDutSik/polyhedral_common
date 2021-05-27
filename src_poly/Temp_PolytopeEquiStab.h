#ifndef INCLUDE_TEMP_POLYTOPE_EQUI_STAB_H
#define INCLUDE_TEMP_POLYTOPE_EQUI_STAB_H


#include "WeightMatrix.h"
#include "WeightMatrixSpecified.h"

//
// Equivalence of subsets and stabilizer of a WeightMatrix
//


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
     if (iRow == n) { // Thus iCol < n.
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
  WeightMatrix<true,int,Tidx_value> WMat1(n+1,[&](size_t iRow, size_t iCol) -> int {
    return g(f1, iRow, iCol);
  });
  WeightMatrix<true,int,Tidx_value> WMat2(n+1,[&](size_t iRow, size_t iCol) -> int {
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
  using Tgr = GraphListAdj;
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
  WeightMatrix<true,int,Tidx_value> WMatW(n+1, g);
  Tgroup GRP=GetStabilizerWeightMatrix<T,Tgr,Tgroup,Tidx_value>(WMatW);
  std::vector<Telt> ListPerm;
  for (auto & ePerm : GRP.GeneratorsOfGroup()) {
    std::vector<Tidx> eList(n);
    for (size_t i=0; i<n; i++)
      eList[i]=OnPoints(i, ePerm);
    ListPerm.push_back(Telt(eList));
  }
  return Tgroup(ListPerm, n);
}

//
// PairOrbits as powerful invariants of subsets extracted from the group action on pairs.
//

template<typename Tgroup, typename Tidx_value>
WeightMatrix<true, int,Tidx_value> WeightMatrixFromPairOrbits(Tgroup const& GRP)
{
  using Telt = typename Tgroup::Telt;
  Tidx_value miss_val = std::numeric_limits<Tidx_value>::max();
  size_t n=GRP.n_act();
  WeightMatrix<true, int, Tidx_value> WMat(n);
  for (size_t i=0; i<n; i++)
    for (size_t j=0; j<n; j++)
      WMat.intDirectAssign(i,j,miss_val);
  auto GetUnset=[&]() -> std::pair<int,int> {
    for (size_t i=0; i<n; i++)
      for (size_t j=0; j<n; j++) {
	Tidx_value eVal=WMat.GetValue(i,j);
	if (eVal == miss_val) {
	  return {i,j};
	}
      }
    return {-1,-1};
  };
  int iOrbit=0;
  std::vector<int> ListWeight;
  std::vector<Telt> ListGen = GRP.GeneratorsOfGroup();
  while(true) {
    std::pair<int,int> eStart=GetUnset();
    if (eStart.first == -1)
      break;
    ListWeight.push_back(iOrbit);
    std::vector<std::pair<int,int>> eList{eStart};
    size_t orbSize = 0;
    while(true) {
      int nbPair = eList.size();
      if (nbPair == 0)
	break;
      orbSize += nbPair;
      std::vector<std::pair<int,int>> fList;
      for (auto & ePair : eList) {
	int i=ePair.first;
	int j=ePair.second;
	WMat.intDirectAssign(i,j,iOrbit);
	for (auto & eGen : ListGen) {
	  int iImg = OnPoints(i, eGen);
	  int jImg = OnPoints(j, eGen);
          Tidx_value eVal1 = WMat.GetValue(iImg,jImg);
          if (eVal1 == miss_val)
            fList.push_back({iImg,jImg});
	}
      }
      eList = std::move(fList);
    }
    iOrbit++;
  }
  WMat.SetWeight(ListWeight);
  return WMat;
}


//
// The antipodal configurations and the absolute trick
//





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


template<typename T, typename Tidx_value>
struct WeightMatrixAbs {
  Tidx_value positionZero;
  Face ArrSigns;
  WeightMatrix<true, T, Tidx_value> WMat;
};


template<typename T, typename Tidx_value>
WeightMatrixAbs<T, Tidx_value> GetSimpleWeightMatrixAntipodal_AbsTrick(MyMatrix<T> const& TheEXT, MyMatrix<T> const& Qmat)
{
  static_assert(is_totally_ordered<T>::value, "Requires T to be a totally ordered field");
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  size_t nbPair=TheEXT.rows();
  size_t nbCol=TheEXT.cols();
  size_t n_ent = (nbPair * (nbPair + 1)) / 2;
  std::vector<Tidx_value> INP_TheMat(n_ent);
  Face ArrSigns(n_ent);
  std::vector<T> INP_ListWeight;
  std::unordered_map<T, Tidx_value> ValueMap;
  Tidx_value idxWeight = 0;
  Tidx_value positionZero = -1;
  //
  auto set_entry=[&](size_t iRow, size_t jRow, Tidx_value pos, bool eChg) -> void {
    size_t idx = weightmatrix_idx<true>(nbPair, iRow, jRow);
    INP_TheMat[idx] = pos;
    ArrSigns[idx] = eChg;
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
  WeightMatrix<true, T, Tidx_value> WMat(nbPair, INP_TheMat, INP_ListWeight);
#ifdef DEBUG
  std::cerr << "Before positionZero=" << positionZero << "\n";
#endif
  positionZero = WMat.ReorderingSetWeight_specificPosition(positionZero);
#ifdef DEBUG
  std::cerr << "After positionZero=" << positionZero << "\n";
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  return {positionZero, std::move(ArrSigns), std::move(WMat)};
}




template<typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value> GetSimpleWeightMatrixAntipodal(MyMatrix<T> const& TheEXT, MyMatrix<T> const& Qmat)
{
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  size_t nbPair=TheEXT.rows();
  size_t nbCol=TheEXT.cols();
  size_t INP_nbRow = 2*nbPair;
  size_t nb = nbPair * (2 * nbPair + 1);
  std::vector<Tidx_value> INP_TheMat(nb);
  std::vector<T> INP_ListWeight;
  std::unordered_map<T, Tidx_value> ValueMap;
  Tidx_value idxWeight = 0;
  //
  auto set_entry=[&](size_t iRow, size_t jRow, Tidx_value pos) -> void {
    size_t idx = weightmatrix_idx<true>(2 * nbPair, iRow, jRow);
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
    }
  }
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetSimpleWeightMatrixAntipodal|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  return WeightMatrix<true, T, Tidx_value>(INP_nbRow, INP_TheMat, INP_ListWeight);
}

template<typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value> GetSimpleWeightMatrix(MyMatrix<T> const& TheEXT, MyMatrix<T> const& Qmat)
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
  return WeightMatrix<true, T, Tidx_value>(nbRow, f1, f2);
}


template<typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value> GetWeightMatrix(MyMatrix<T> const& TheEXT)
{
  MyMatrix<T> Qmat=GetQmatrix(TheEXT);
  return GetSimpleWeightMatrix<T,Tidx_value>(TheEXT, Qmat);
}





template<typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value> GetWeightMatrixAntipodal(MyMatrix<T> const& TheEXT)
{
  MyMatrix<T> Qmat=GetQmatrix(TheEXT);
  return GetSimpleWeightMatrixAntipodal<T,Tidx_value>(TheEXT, Qmat);
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
  using Tidx_value = int16_t;
  using Tgr=GraphBitset;
  size_t nbRow= EXT.rows();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  WeightMatrixAbs<Tint,Tidx_value> WMatAbs = GetSimpleWeightMatrixAntipodal_AbsTrick<Tint,Tidx_value>(EXT, Qmat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  //  std::cerr << "WMatAbs.positionZero=" << WMatAbs.positionZero << "\n";
  //  std::cerr << "WMatAbs.WMat=\n";
  //  PrintWeightedMatrix(std::cerr, WMatAbs.WMat);

  Tgr eGR=GetGraphFromWeightedMatrix<Tint,Tgr>(WMatAbs.WMat);
  //  GRAPH_PrintOutputGAP_vertex_colored("GAP_graph", eGR);

  //  std::cerr << "WMatAbs.WMat : ";
  //  PrintStabilizerGroupSizes(std::cerr, eGR);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|GetGraphFromWeightedMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif

  using Tidx=unsigned int;
  int nbVert_G = eGR.GetNbVert();
#ifdef USE_BLISS
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> ePair = BLISS_GetCanonicalOrdering_ListGenerators<Tgr,Tidx>(eGR, nbVert_G);
#endif
#ifdef USE_TRACES
  std::pair<std::vector<Tidx>, std::vector<std::vector<Tidx>>> ePair = TRACES_GetCanonicalOrdering_ListGenerators<Tgr,Tidx>(eGR, nbVert_G);
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
              size_t idx1 = weightmatrix_idx<true>(nbRow, i, j);
              size_t idx2 = weightmatrix_idx<true>(nbRow, iImg, jImg);
              bool ChgSign1 = WMatAbs.ArrSigns[idx1];
              bool ChgSign2 = WMatAbs.ArrSigns[idx2];
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
  std::vector<Tidx> CanonicOrd = GetCanonicalizationVector_KernelBis<Tidx>(nbRow, ePair.first);
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
    int i_row_orig = CanonicOrd[i_row];
    for (size_t k_row=0; k_row<nbRow; k_row++) {
      if (k_row != i_row && ListSigns[k_row] != 0) {
        int k_row_orig = CanonicOrd[k_row];
        if (WMatAbs.WMat.GetValue(i_row_orig, k_row_orig) != WMatAbs.positionZero) {
          size_t idx = weightmatrix_idx<true>(nbRow, i_row_orig, k_row_orig);
          bool ChgSign = WMatAbs.ArrSigns[idx];
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
    int i_rowC = CanonicOrd[i_row];
    for (size_t j_row=0; j_row<nbRow; j_row++) {
      int j_rowC = CanonicOrd[j_row];
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
    int j_row = CanonicOrd[i_row];
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
MyMatrix<Tint> LinPolytopeAntipodalIntegral_CanonicForm(MyMatrix<Tint> const& EXT)
{
  using Tidx_value = int16_t;
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

  WeightMatrix<true, Tint, Tidx_value> WMat=GetWeightMatrixAntipodal<Tint, Tidx_value>(EXT);
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

  std::vector<int> CanonicOrd = GetCanonicalizationVector<Tint,GraphBitset,int>(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time6 = std::chrono::system_clock::now();
  std::cerr << "|GetCanonicalizationVector|=" << std::chrono::duration_cast<std::chrono::microseconds>(time6 - time5).count() << "\n";
#endif

  MyMatrix<Tint> EXTreord(n_rows, n_cols);
  size_t idx=0;
  Face IsIncluded(n_rows);
  for (size_t i_row=0; i_row<2*n_rows; i_row++) {
    int j_row = CanonicOrd[i_row];
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
EquivTest<std::vector<std::vector<unsigned int>>> LinPolytopeAntipodalIntegral_Automorphism_AbsTrick(MyMatrix<Tint> const& EXT, MyMatrix<Tint> const& Qmat)
{
  using Tidx_value = int16_t;
  using Tgr=GraphBitset;
  size_t nbRow= EXT.rows();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif
  WeightMatrixAbs<Tint,Tidx_value> WMatAbs = GetSimpleWeightMatrixAntipodal_AbsTrick<Tint,Tidx_value>(EXT, Qmat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetSimpleWeightMatrixAntipodal_AbsTrick|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif
  //  std::cerr << "WMatAbs.positionZero=" << WMatAbs.positionZero << "\n";
  //  std::cerr << "WMatAbs.WMat=\n";
  //  PrintWeightedMatrix(std::cerr, WMatAbs.WMat);

  Tgr eGR=GetGraphFromWeightedMatrix<Tint,Tgr>(WMatAbs.WMat);
  //  GRAPH_PrintOutputGAP_vertex_colored("GAP_graph", eGR);

  //  std::cerr << "WMatAbs.WMat : ";
  //  PrintStabilizerGroupSizes(std::cerr, eGR);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|GetGraphFromWeightedMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif

  using Tidx = unsigned int;
  int nbVert_G = eGR.GetNbVert();
#ifdef USE_BLISS
  std::vector<std::vector<Tidx>> ListGen = BLISS_GetListGenerators<Tgr,Tidx>(eGR, nbVert_G);
#endif
#ifdef USE_TRACES
  std::vector<std::vector<Tidx>> ListGen = TRACES_GetListGenerators<Tgr,Tidx>(eGR, nbVert_G);
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
              size_t idx1 = weightmatrix_idx<true>(nbRow, i, j);
              size_t idx2 = weightmatrix_idx<true>(nbRow, iImg, jImg);
              bool ChgSign1 = WMatAbs.ArrSigns[idx1];
              bool ChgSign2 = WMatAbs.ArrSigns[idx2];
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



template<typename Tint>
std::vector<std::vector<unsigned int>> LinPolytopeAntipodalIntegral_Automorphism(MyMatrix<Tint> const& EXT)
{
  using Tidx_value = int16_t;
  using Tgr = GraphBitset;
  int nbRow = EXT.rows();
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

  WeightMatrix<true, Tint, Tidx_value> WMat=GetWeightMatrixAntipodal<Tint, Tidx_value>(EXT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|GetWeightMatrixAntipodal|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif

  Tgr eGR=GetGraphFromWeightedMatrix<Tint,Tgr>(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time5 = std::chrono::system_clock::now();
  std::cerr << "|GetGraphFromWeightedMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time5 - time4).count() << "\n";
#endif

  using Tidx=unsigned int;
#ifdef USE_BLISS
  std::vector<std::vector<Tidx>> ListGen = BLISS_GetListGenerators<Tgr,Tidx>(eGR, nbRow);
#endif
#ifdef USE_TRACES
  std::vector<std::vector<Tidx>> ListGen = TRACES_GetListGenerators<Tgr,Tidx>(eGR, nbRow);
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time6 = std::chrono::system_clock::now();
  std::cerr << "|GetListGenerators|=" << std::chrono::duration_cast<std::chrono::microseconds>(time6 - time5).count() << "\n";
#endif
  return ListGen;
}



//
// Various construction of weighted matrix
//










template<typename T, typename Tidx_value>
WeightMatrix<false, T, Tidx_value> T_TranslateToMatrix(MyMatrix<T> const& eMat)
{
  size_t nbRow=eMat.rows();
  auto f=[&](size_t iRow, size_t iCol) -> T {
    return eMat(iRow,iCol);
  };
  return WeightMatrix<false, T, Tidx_value>(nbRow, f);
}



// The matrices in ListMat do not have to be symmetric.
template<typename T, typename Tint, typename Tidx_value>
WeightMatrix<false, std::vector<T>, Tidx_value> T_TranslateToMatrix_ListMat_SHV(std::vector<MyMatrix<T>> const& ListMat, MyMatrix<Tint> const& SHV)
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
  return WeightMatrix<false, std::vector<T>, Tidx_value>(nbRow, f1, f2);
}




template<bool is_symmetric, typename T, typename Tidx_value>
WeightMatrix<is_symmetric, std::vector<T>, Tidx_value> GetWeightMatrix_ListComm(MyMatrix<T> const& TheEXT, MyMatrix<T> const&GramMat, std::vector<MyMatrix<T>> const& ListComm)
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
  MyMatrix<T> M(nbComm+1, nbCol);
  auto f1=[&](size_t iRow) -> void {
    for (size_t iMat=0; iMat<=nbComm; iMat++) {
      for (size_t iCol=0; iCol<nbCol; iCol++) {
        T eSum = 0;
        for (size_t jCol=0; jCol<nbCol; jCol++)
          eSum += ListProd[iMat](iCol, jCol) * TheEXT(iRow, jCol);
        M(iMat, iCol) = eSum;
      }
    }
  };
  std::vector<T> eVectSum(nbComm+1);
  auto f2=[&](size_t jRow) -> std::vector<T> {
    for (size_t iMat=0; iMat<=nbComm; iMat++) {
      T eSum=0;
      for (size_t iCol=0; iCol<nbCol; iCol++)
        eSum += TheEXT(jRow, iCol) * M(iMat, iCol);
      eVectSum[iMat] = eSum;
    }
    return eVectSum;
  };
  return WeightMatrix<false, std::vector<T>, Tidx_value>(nbRow, f1, f2);
}



template<typename T, typename Tidx_value>
WeightMatrix<false, std::vector<T>, Tidx_value> GetWeightMatrix_ListMatrix(std::vector<MyMatrix<T>> const& ListMatrix, MyMatrix<T> const& TheEXT)
{
  size_t nbRow=TheEXT.rows();
  size_t nbCol=TheEXT.cols();
  size_t nbMat=ListMatrix.size();
  MyMatrix<T> M(nbMat, nbCol);
  auto f1=[&](size_t iRow) -> void {
    for (size_t iMat=0; iMat<nbMat; iMat++) {
      for (size_t iCol=0; iCol<nbCol; iCol++) {
        T eSum = 0;
        for (size_t jCol=0; jCol<nbCol; jCol++)
          eSum += ListMatrix[iMat](iCol, jCol) * TheEXT(iRow, jCol);
        M(iMat, iCol) = eSum;
      }
    }
  };
  std::vector<T> eVectScal(nbMat);
  auto f2=[&](size_t jRow) -> std::vector<T> {
    for (size_t iMat=0; iMat<nbMat; iMat++) {
      T eSum = 0;
      for (size_t iCol=0; iCol<nbCol; iCol++)
        eSum += TheEXT(jRow, iCol) * M(iMat, iCol);
      eVectScal[iMat] = eSum;
    }
    return eVectScal;
  };
  return WeightMatrix<false, std::vector<T>, Tidx_value>(nbRow, f1, f2);
}



template<typename T, typename Tidx_value>
WeightMatrix<true, T, Tidx_value> GetWeightMatrixGramMatShort(MyMatrix<T> const& TheGramMat, MyMatrix<int> const& ListShort)
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
  return WeightMatrix<true, T, Tidx_value>(nbShort, f1, f2);
}




template<typename T, typename Tint, typename Tidx_value>
WeightMatrix<true, T, Tidx_value> T_TranslateToMatrix_QM_SHV(MyMatrix<T> const& qMat, MyMatrix<Tint> const& SHV)
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
  return WeightMatrix<true, T, Tidx_value>(INP_nbRow, INP_TheMat, INP_ListWeight);
}






template<typename T, typename Tgroup>
Tgroup LinPolytope_Automorphism(MyMatrix<T> const & EXT)
{
  using Tgr = GraphListAdj;
  using Tidx_value = int16_t;
  MyMatrix<T> EXTred=ColumnReduction(EXT);
  WeightMatrix<true, T, Tidx_value> WMat=GetWeightMatrix<T,Tidx_value>(EXTred);
  return GetStabilizerWeightMatrix<T,Tgr,Tgroup,Tidx_value>(WMat);
}


// ListMat is assumed to be symmetric
template<typename T, typename Tidx, typename Tidx_value>
WeightMatrix<true, std::vector<T>, Tidx_value> GetWeightMatrix_ListMat_Subset(MyMatrix<T> const& TheEXT, std::vector<MyMatrix<T>> const& ListMat, Face const& eSubset)
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
  using Tfield = typename overlying_field<T>::field_type;
  auto f3=[&](std::vector<int> const& Vsubset) -> bool {
    if (Vsubset.size() < nbCol)
      return false;
    auto f=[&](MyMatrix<Tfield> & M, size_t eRank, size_t iRow) -> void {
      for (size_t iCol=0; iCol<nbCol; iCol++)
        M(eRank, iCol) = UniversalTypeConversion<Tfield,T>(TheEXT(Vsubset[iRow], iCol));
    };
    SelectionRowCol<Tfield> TheSol = TMat_SelectRowCol_Kernel(Vsubset.size(), nbCol, f);
    return TheSol.TheRank == nbCol;
  };
  auto f4=[&](std::vector<Tidx> const& Vsubset, std::vector<Tidx> const& Vin) -> EquivTest<std::vector<Tidx>> {
    auto g1=[&](size_t iRow) -> MyVector<T> {
      MyVector<T> V(nbCol);
      for (size_t iCol=0; iCol<nbCol; iCol++)
        V(iCol) = TheEXT(Vsubset[iRow], iCol);
      return V;
    };
    EquivTest<MyMatrix<Tfield>> test1 = RepresentVertexPermutationTest(Vsubset.size(), nbCol, g1, g1, Vin);
    if (!test1.TheReply)
      return {false, {}};
    EquivTest<std::vector<Tidx>> test2 = RepresentVertexPermutationTest(TheEXT, TheEXT, test1.TheEquiv);
    return test2;
  };
  return WeightMatrix<true, std::vector<T>, Tidx_value>(nbRow, f1, f2);
}



template<typename T>
size_t GetInvariant_ListMat_Subset(MyMatrix<T> const& EXT, std::vector<MyMatrix<T>> const&ListMat, Face const& eSubset)
{
  using Tidx_value = int16_t;
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif


  WeightMatrix<true, std::vector<T>, Tidx_value> WMat = GetWeightMatrix_ListMat_Subset<T,Tidx_value>(EXT, ListMat, eSubset);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetWeightMatrix_ListMatrix_Subset|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif


  WMat.ReorderingSetWeight();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|ReorderingSetWeight|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif


  size_t e_hash = std::hash<WeightMatrix<true, std::vector<T>, Tidx_value>>()(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|hash|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif
  return e_hash;
}

template<typename T>
std::vector<std::vector<unsigned int>> GetListGenAutomorphism_ListMat_Subset(MyMatrix<T> const& EXT, std::vector<MyMatrix<T>> const&ListMat, Face const& eSubset)
{
  using Tidx_value = int16_t;
  //  using Tgr = GraphBitset;
  using Tgr = GraphListAdj;
  size_t nbRow = EXT.rows();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif


  WeightMatrix<true, std::vector<T>, Tidx_value> WMat = GetWeightMatrix_ListMat_Subset<T,Tidx_value>(EXT, ListMat, eSubset);
  // No need to reorder in autom case
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
  std::cerr << "|GetWeightMatrix_ListMatrix_Subset|=" << std::chrono::duration_cast<std::chrono::microseconds>(time2 - time1).count() << "\n";
#endif

  Tgr eGR=GetGraphFromWeightedMatrix<std::vector<T>,Tgr>(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
  std::cerr << "|GetGraphFromWeightMatrix|=" << std::chrono::duration_cast<std::chrono::microseconds>(time3 - time2).count() << "\n";
#endif

  using Tidx = unsigned int;
#ifdef USE_BLISS
  std::vector<std::vector<Tidx>> ListGen = BLISS_GetListGenerators<Tgr,Tidx>(eGR, nbRow);
#endif
#ifdef USE_TRACES
  std::vector<std::vector<Tidx>> ListGen = TRACES_GetListGenerators<Tgr,Tidx>(eGR, nbRow);
#endif
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
  std::cerr << "|GetListGenerators|=" << std::chrono::duration_cast<std::chrono::microseconds>(time4 - time3).count() << "\n";
#endif
  return ListGen;
}






template<typename T>
EquivTest<std::vector<unsigned int>> TestEquivalence_ListMat_Subset(
                                       MyMatrix<T> const& EXT1, std::vector<MyMatrix<T>> const&ListMat1, Face const& eSubset1,
                                       MyMatrix<T> const& EXT2, std::vector<MyMatrix<T>> const&ListMat2, Face const& eSubset2)
{
  using Tidx_value = int16_t;
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif

  WeightMatrix<true, std::vector<T>, Tidx_value> WMat1 = GetWeightMatrix_ListMat_Subset<T,Tidx_value>(EXT1, ListMat1, eSubset1);
  WeightMatrix<true, std::vector<T>, Tidx_value> WMat2 = GetWeightMatrix_ListMat_Subset<T,Tidx_value>(EXT2, ListMat2, eSubset2);
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
  using Tidx_value = int16_t;
  using Tgr=GraphBitset;
  size_t n_rows = EXT.rows();
  size_t n_cols = EXT.cols();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time1 = std::chrono::system_clock::now();
#endif

  WeightMatrix<true, Tint, Tidx_value> WMat=GetWeightMatrix<Tint,Tidx_value>(EXT);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time2 = std::chrono::system_clock::now();
#endif

  WMat.ReorderingSetWeight();
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time3 = std::chrono::system_clock::now();
#endif

  std::vector<int> CanonicOrd = GetCanonicalizationVector<Tint,Tgr,int>(WMat);
#ifdef TIMINGS
  std::chrono::time_point<std::chrono::system_clock> time4 = std::chrono::system_clock::now();
#endif

  MyMatrix<Tint> EXTreord(n_rows, n_cols);
  for (size_t i_row=0; i_row<n_rows; i_row++) {
    size_t j_row = CanonicOrd[i_row];
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







#endif
