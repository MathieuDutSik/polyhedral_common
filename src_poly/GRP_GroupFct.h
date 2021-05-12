#ifndef INCLUDE_TEMP_GROUP_FCT_H
#define INCLUDE_TEMP_GROUP_FCT_H

#include "Temp_common.h"
#include "Basic_file.h"
#include "Basic_string.h"
#include "hash_functions.h"
#include "COMB_Stor.h"
#include "NumberTheory.h"
#include "Boost_bitset.h"
#include "Temp_PolytopeEquiStab.h"



//
// permutation functions
//

template<typename Telt>
std::vector<int> PermutationOrbit(Telt const& ePerm)
{
  int siz=ePerm.size();
  std::vector<int> StatusOrbit(siz,-1);
  int idxOrbit=0;
  auto GetUnsetPoint=[&](void) -> int {
    for (int i=0; i<siz; i++)
      if (StatusOrbit[i] == -1)
	return i;
    return -1;
  };
  while(true) {
    int iPoint=GetUnsetPoint();
    if (iPoint == -1)
      break;
    int iPointWork=iPoint;
    while(true) {
      StatusOrbit[iPointWork] = idxOrbit;
      iPointWork = OnPoints(iPointWork, ePerm);
      if (iPointWork == iPoint)
	break;
    }
    idxOrbit++;
  }
  return StatusOrbit;
}


template<typename Telt>
void WritePermutationGAP(std::ostream&os, Telt const& ePerm)
{
  if (ePerm.isIdentity() ) {
    os << "()";
    return;
  }
  std::vector<int> eVectOrbit=PermutationOrbit(ePerm);
  int nbOrbit=VectorMax(eVectOrbit) + 1;
  int siz=eVectOrbit.size();
  for (int iOrbit=0; iOrbit<nbOrbit; iOrbit++) {
    int nbMatch=0;
    int eFirst=-1;
    for (int i=0; i<siz; i++)
      if (eVectOrbit[i] == iOrbit) {
	if (nbMatch == 0) {
	  eFirst=i;
	}
	nbMatch++;
      }
    std::vector<int> TheList(nbMatch);
    TheList[0]=eFirst;
    for (int i=1; i<nbMatch; i++)
      TheList[i] = OnPoints(TheList[i-1], ePerm);
    os << "(";
    for (int i=0; i<nbMatch; i++) {
      if (i>0)
	os << ",";
      int val=TheList[i]+1;
      os << val;
    }
    os << ")";
  }
}

//
// Template general code for input output of groups
//

template<typename Tgroup>
Tgroup ReadGroup(std::istream &is)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  if (!is.good()) {
    std::cerr << "ReadGroup operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  int nbGen;
  int n;
  is >> n;
  is >> nbGen;
  std::cerr << "n=" << n << " nbGen=" << nbGen << "\n";
  std::vector<Telt> ListGen;
  for (int iGen=0; iGen<nbGen; iGen++) {
    std::vector<Tidx> v(n);
    for (int i=0; i<n; i++) {
      int eVal;
      is >> eVal;
      if (eVal > n-1) {
	std::cerr << "Error in ReadGroup function\n";
	std::cerr << "Number of elements acted on n=" << n << "\n";
	std::cerr << "But eVal=" << eVal << "\n";
	throw TerminalException{1};
      }
      v[i]=eVal;
    }
    ListGen.push_back(Telt(v));
  }
  return Tgroup(ListGen, n);
}

template<typename Tgroup>
void WriteGroup(std::ostream &os, Tgroup const& TheGRP)
{
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  int nbGen=ListGen.size();
  os << TheGRP.n_act() << " " << nbGen << "\n";
  for (auto & eGen : ListGen) {
    for (int i=0; i<TheGRP.n_act(); i++) {
      int eVal=OnPoints(i, eGen);
      os << " " << eVal;
    }
    os << "\n";
  }
}


template<typename Tgroup>
void WriteGroupMakeUp(std::ostream &os, Tgroup const& TheGRP)
{
  using Telt = typename Tgroup::Telt;
  os << "nb acting element=" << TheGRP.n_act() << "\n";
  std::vector<Telt> ListGen=TheGRP.GeneratorsOfGroup();
  int nbGen=0;
  for (auto & eGen : ListGen) {
    for (std::size_t i=0; i<TheGRP.n_act(); i++) {
      int eVal = OnPoints(i, eGen);
      os << " " << eVal;
    }
    os << "\n";
    nbGen++;
  }
  os << "nbGen=" << nbGen << "\n";
  os << "size=" << TheGRP.size() << "\n";
}

template<typename Tgroup>
void WriteGroupGAP(std::ostream &os, Tgroup const& TheGRP)
{
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> ListGen=TheGRP.GeneratorsOfGroup();
  os << "local eListList, ListGen, GRP;\n";
  os << "eListList:=[\n";
  bool IsFirst=true;
  for (auto & eGen : ListGen) {
    if (!IsFirst)
      os << ",\n";
    IsFirst=false;
    os << "[";
    for (int i=0; i<TheGRP.n_act(); i++) {
      int eVal = 1 + OnPoints(i, eGen);
      if (i>0)
	os << ",";
      os << eVal;
    }
    os << "]";
  }
  os << "];\n";
  os << "ListGen:=List(eListList, PermList);\n";
  os << "GRP:=Group(ListGen);\n";
  os << "SetSize(GRP, " << TheGRP.size() << ");\n";
  os << "return GRP;\n";
}

//
// group combinatorial algorithms
//


template<typename Tgroup>
std::vector<int> OrbitIntersection(Tgroup const& TheGRP, std::vector<int> const& gList)
{
  std::vector<int> eListReturn;
  int n=TheGRP.n_act();
  std::vector<int> rList = gList;
  auto LGen = TheGRP.GeneratorsOfGroup();
  while(true) {
    std::size_t eSumPrev=0;
    for (int i=0; i<n; i++)
      eSumPrev += rList[i];
    for (auto & eGen : LGen) {
      for (int i=0; i<n; i++) {
	int j = OnPoints(i, eGen);
	if (rList[i] == 0)
	  rList[j]=0;
      }
    }
    std::size_t eSum=0;
    for (int i=0; i<n; i++)
      eSum += rList[i];
    if (eSum == eSumPrev)
      break;
  }
  return rList;
}


template<typename Tgroup>
std::vector<int> OrbitUnion(Tgroup const& TheGRP, std::vector<int> const& gList)
{
  int n=TheGRP.n_act();
  std::vector<int> gListB(n);
  for (int i=0; i<n; i++)
    gListB[i] = 1 - gList[i];
  std::vector<int> rListB=OrbitIntersection(TheGRP, gListB);
  for (int i=0; i<n; i++)
    rListB[i] = 1 - rListB[i];
  return rListB;
}

template<typename Tgroup>
Face OrbitIntersection(Tgroup const& GRP, Face const& gList)
{
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> LGen = GRP.GeneratorsOfGroup();
  int n = GRP.n_act();
  Face rList = gList;
  while(true) {
    size_t eSumPrev = rList.count();
    for (auto & eGen : LGen) {
      for (int i=0; i<n; i++) {
        std::size_t j = OnPoints(i, eGen);
        if (rList[i] == 0)
          rList[j]=0;
      }
    }
    size_t eSum = rList.count();
    if (eSum == eSumPrev)
      break;
  }
  return rList;
}


template<typename Tgroup>
Face OrbitUnion(Tgroup const& GRP, Face const& gList)
{
  int n = GRP.n_act();
  Face gListB(n);
  for (int i=0; i<n; i++)
    gListB[i] = 1 - gList[i];
  Face rListB=OrbitIntersection(GRP, gListB);
  for (int i=0; i<n; i++)
    rListB[i] = 1 - rListB[i];
  return rListB;
}






//
// Several building of new groups.
//

template<typename Tgroup>
Tgroup ReducedGroupAction(Tgroup const& TheGRP, Face const& eList)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  int nb=eList.count();
  if (nb == 0) {
    std::cerr << "Call of ReducedGroupAction with 0 points\n";
    throw TerminalException{1};
  }
  std::vector<int> ListPositionRev(TheGRP.n_act(), -1);
  int aRow=eList.find_first();
  std::vector<int> ListPosition(nb,-1);
  for (int iRow=0; iRow<nb; iRow++) {
    ListPositionRev[aRow]=iRow;
    ListPosition[iRow]=aRow;
    aRow=eList.find_next(aRow);
  }
  std::vector<Telt> ListGen;
  for (auto & eGen : TheGRP.GeneratorsOfGroup()) {
    std::vector<Tidx> v(nb);
    for (int i=0; i<nb; i++) {
      int eVal1=ListPosition[i];
      int eVal2=OnPoints(eVal1, eGen);
      int eVal3=ListPositionRev[eVal2];
      v[i]=eVal3;
    }
    ListGen.push_back(Telt(v));
  }
  return Tgroup(ListGen, nb);
}


template<typename Tgroup>
Tgroup ConjugateGroup(Tgroup const& TheGRP, typename Tgroup::Telt const& ePerm)
{
  using Telt = typename Tgroup::Telt;
  int n=TheGRP.n_act();
  Telt ePermR=~ePerm;
  std::vector<Telt> ListGen;
  for (auto & eGen : TheGRP.GeneratorsOfGroup()) {
    std::vector<int> v(n);
    for (int i=0; i<n; i++) {
      int eVal1=OnPoints(i, ePermR);
      int eVal2=OnPoints(eVal1, eGen);
      int eVal3=OnPoints(eVal2, ePerm);
      v[i]=eVal3;
    }
    ListGen.push_back(Telt(v));
  }
  return Tgroup(ListGen, n);
}


//
// Some enumeration code
//

template<typename Tgroup>
void GROUP_FuncInsertInSet(Tgroup const& TheGRP, Face const& eList, vectface &ListListSet)
{
  int nb=ListListSet.size();
  for (int iList=0; iList<nb; iList++) {
    bool test=TheGRP.RepresentativeAction_OnSets(eList, ListListSet[iList]).first;
    if (test)
      return;
  }
  ListListSet.push_back(eList);
}


template<typename Tgroup>
void GROUP_FuncInsertInSet_UseInv(Tgroup const& TheGRP,
				  Face const& eList,
				  std::vector<int> const& eInv, 
				  vectface & ListSet,
				  std::vector<std::vector<int>> & ListInv)
{
  int nb=ListSet.size();
  for (int iList=0; iList<nb; iList++)
    if (eInv == ListInv[iList]) {
      bool test=TheGRP.RepresentativeAction_OnSets(eList, ListSet[iList]).first;
      if (test)
	return;
    }
  ListSet.push_back(eList);
  ListInv.push_back(eInv);
}


//
// Some combinatorial algorithms using only the generators of the group.
//


template<typename Tgroup>
std::vector<int> ComputeFullOrbitPoint(Tgroup const& TheGRP, int const& ePoint)
{
  using Telt = typename Tgroup::Telt;
  int n = TheGRP.n_act();
  IntegerSubsetStorage Vorb = VSLT_InitializeStorage(n);
  IntegerSubsetStorage Vactive = VSLT_InitializeStorage(n);
  std::vector<int> eList;
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  VSLT_StoreValue(Vactive, ePoint);
  while(true) {
    if (VSLT_IsEmpty(Vactive))
      break;
    int TheFirst=VSLT_TheFirstPosition(Vactive);
    VSLT_RemoveValue(Vactive, TheFirst);
    VSLT_StoreValue(Vorb, TheFirst);
    eList.push_back(TheFirst);
    for (auto & eGen : ListGen) {
      int NewPt = OnPoints(TheFirst, eGen);
      if (!VSLT_IsItInSubset(Vorb, NewPt) && !VSLT_IsItInSubset(Vactive, NewPt))
	VSLT_StoreValue(Vactive, NewPt);
    }
  }
  return eList;
}

template<typename Tgroup>
vectface DecomposeOrbitPoint(Tgroup const& TheGRP, Face const& eList)
{
  int nbPoint=TheGRP.n_act();
  IntegerSubsetStorage Vlist = VSLT_InitializeStorage(nbPoint);
  vectface ListOrb(nbPoint);
  int len=eList.count();
  int aRow=eList.find_first();
  for (int i=0; i<len; i++) {
    VSLT_StoreValue(Vlist, aRow);
    aRow=eList.find_next(aRow);
  }
  while(true) {
    if (VSLT_IsEmpty(Vlist))
      break;
    int TheFirst=VSLT_TheFirstPosition(Vlist);
    std::vector<int> eOrb=ComputeFullOrbitPoint(TheGRP, TheFirst);
    Face vectOrb(nbPoint);
    for (auto & ePt : eOrb) {
      vectOrb[ePt]=1;
      VSLT_RemoveValue(Vlist, ePt);
    }
    ListOrb.push_back(vectOrb);
  }
  return ListOrb;
}

template<typename Tobj, typename Tgen>
std::vector<Tobj> OrbitSplittingGeneralized(std::vector<Tobj> const& PreListTotal,
					    std::vector<Tgen> const& ListGen,
					    std::function<Tobj(Tobj const&,Tgen const&)> const& TheAct)
{
  std::vector<Tobj> TheReturn;
  std::unordered_set<Tobj> ListTotal;
  for (auto eObj : PreListTotal)
    ListTotal.insert(eObj);
  while(true) {
    auto iter=ListTotal.begin();
    if (iter == ListTotal.end())
      break;
    Tobj eObj=*iter;
    TheReturn.push_back(eObj);
    std::unordered_set<Tobj> Additional;
    Additional.insert(eObj);
    ListTotal.erase(eObj);
    std::unordered_set<Tobj> SingleOrbit;
    while(true) {
      std::unordered_set<Tobj> NewElts;
      for (auto const& gObj : Additional) {
	for (auto const& eGen : ListGen) {
	  Tobj fObj=TheAct(gObj, eGen);
	  if (SingleOrbit.count(fObj) == 0 && Additional.count(fObj) == 0) {
	    if (NewElts.count(fObj) == 0) {
#ifdef DEBUG
	      if (ListTotal.count(fObj) > 0) {
		NewElts.insert(fObj);
		ListTotal.erase(fObj);
	      } else {
		std::cerr << "Orbit do not match, PANIC!!!\n";
		throw TerminalException{1};
	      }
#else
              NewElts.insert(fObj);
              ListTotal.erase(fObj);
#endif
	    }
	  }
	}
      }
      for (auto & uSet : Additional)
	SingleOrbit.insert(uSet);
      if (NewElts.size() == 0)
	break;
      Additional = std::move(NewElts);
    }
  }
  return TheReturn;
}



template<typename Telt>
Face OnFace(Face const& eSet, Telt const& eElt)
{
  int nbExt=eSet.size();
  Face fSet(nbExt);
  boost::dynamic_bitset<>::size_type pos=eSet.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    int jExt=eElt.at(pos);
    fSet[jExt]=1;
    pos = eSet.find_next(pos);
  }
  return fSet;
}



template<typename Tgroup>
vectface OrbitSplittingSet(vectface const& PreListTotal, Tgroup const& TheGRP)
{
  using Telt = typename Tgroup::Telt;
  vectface TheReturn(TheGRP.n_act());
  std::unordered_set<Face> ListTotal;
  for (auto eFace : PreListTotal)
    ListTotal.insert(eFace);
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  while(true) {
    std::unordered_set<Face>::iterator iter=ListTotal.begin();
    if (iter == ListTotal.end())
      break;
    Face eSet=*iter;
    TheReturn.push_back(eSet);
    std::unordered_set<Face> Additional{eSet};
    ListTotal.erase(eSet);
    std::unordered_set<Face> SingleOrbit;
    while(true) {
      std::unordered_set<Face> NewElts;
      for (auto const& gSet : Additional) {
	for (auto const& eGen : ListGen) {
	  Face fSet=OnFace(gSet, eGen);
	  if (SingleOrbit.count(fSet) == 0 && Additional.count(fSet) == 0) {
	    if (NewElts.count(fSet) == 0) {
#ifdef DEBUG
	      if (ListTotal.count(fSet) > 0) {
		NewElts.insert(fSet);
		ListTotal.erase(fSet);
	      } else {
		std::cerr << "Orbit do not matched, PANIC!!!\n";
		throw TerminalException{1};
	      }
#else
              NewElts.insert(fSet);
              ListTotal.erase(fSet);
#endif
	    }
	  }
	}
      }
      for (auto & uSet : Additional)
	SingleOrbit.insert(uSet);
      if (NewElts.size() == 0)
	break;
      Additional = std::move(NewElts);
    }
  }
  return TheReturn;
}


//
// Double Coset Computation
//


template<typename Tgroup>
vectface DoubleCosetDescription_Representation(Tgroup const& BigGRP, Tgroup const& SmaGRP, LocalInvInfo const& LocalInv, Face const& eList, std::ostream & os)
{
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  //  os << "Beginning of DoubleCosetDescription\n";
  std::vector<Telt> ListGen=BigGRP.GeneratorsOfGroup();
  Tgroup TheStab=BigGRP.Stabilizer_OnSets(eList);
  //  os << "BigGRP.size=" << BigGRP.size() << " TheStab.size=" << TheStab.size() << "\n";
  Tint TotalSize=BigGRP.size() / TheStab.size();
  //  os << "TotalSize=" << TotalSize << "\n";
  //
  struct Local {
    int status;
    Face eFace;
    std::vector<int> eInv;
  };
  Tint SizeGen=0;
  std::vector<Local> ListLocal;
  auto DoubleCosetInsertEntry=[&](Face const& testList) -> void {
    std::vector<int> eInv=GetLocalInvariantWeightMatrix_Enhanced<int>(LocalInv, testList);
    for (auto const& fLocal : ListLocal) {
      bool test = SmaGRP.RepresentativeAction_OnSets(fLocal.eFace, testList).first;
      if (test)
	return;
    }
    ListLocal.push_back({0,testList,eInv});
    Tgroup fStab=SmaGRP.Stabilizer_OnSets(testList);
    Tint OrbSizeSma=SmaGRP.size() / fStab.size();
    SizeGen += OrbSizeSma;
  };
  DoubleCosetInsertEntry(eList);
  while(true) {
    bool DoSomething=false;
    size_t nbLocal=ListLocal.size();
    for (size_t iLocal=0; iLocal<nbLocal; iLocal++)
      if (ListLocal[iLocal].status == 0) {
	ListLocal[iLocal].status=1;
	DoSomething=true;
	Face eFace=ListLocal[iLocal].eFace;
	for (auto const& eGen : ListGen) {
	  Face eNewList=OnFace(eFace, eGen);
	  DoubleCosetInsertEntry(eNewList);
	}
      }
    if (!DoSomething)
      break;
  }
  vectface ListListSet;
  for (auto & eRec : ListLocal)
    ListListSet.push_back(eRec.eFace);
  if (SizeGen == TotalSize)
    return ListListSet;
  os << "After Iteration loop SizeGen=" << SizeGen << " TotalSize=" << TotalSize << "\n";
  vectface PartialOrbit = std::move(ListListSet);
  auto IsPresent=[&](Face const& testList) -> bool {
    for (auto & fList : PartialOrbit)
      if (fList == testList)
	return true;
    return false;
  };
  size_t pos_start=0;
  while(true) {
    size_t n_orb = PartialOrbit.size();
    for (size_t i_orb=pos_start; i_orb<n_orb; i_orb++) {
      for (auto & eGen : ListGen) {
	Face eNewList=OnFace(PartialOrbit[i_orb], eGen);
	if (!IsPresent(eNewList)) {
	  PartialOrbit.push_back(eNewList);
	  DoubleCosetInsertEntry(eNewList);
	  if (SizeGen == TotalSize) {
	    vectface ListListFin;
	    for (auto & eRec : ListLocal)
	      ListListFin.push_back(eRec.eFace);
	    return ListListFin;
	  }
	}
      }
    }
    pos_start = n_orb;
  }
  os << "Likely not reachable stage\n";
  throw TerminalException{1};
}


template<typename T>
struct popable_vector {
private:
  std::vector<T> V;
  size_t pos;
public:
  popable_vector() : pos(0) {}
  popable_vector(std::vector<T> const& _V) : V(_V), pos(_V.size()) {}
  size_t size() const
  {
    return pos;
  }
  void push_back(T const& val)
  {
    if (V.size() == pos) {
      V.push_back(val);
      pos++;
      return;
    }
    V[pos] = val;
    pos++;
  }
  T pop()
  {
    pos--;
    return V[pos];
  }
};


template<typename Tgroup>
vectface DoubleCosetDescription_Canonic(Tgroup const& BigGRP, Tgroup const& SmaGRP, Face const& eList, std::ostream & os)
{
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  std::vector<Telt> ListGen=BigGRP.GeneratorsOfGroup();
  Tgroup TheStab=BigGRP.Stabilizer_OnSets(eList);
  Tint TotalSize=BigGRP.size() / TheStab.size();
  //
  Tint SizeGen=0;
  auto IncreaseSize=[&](Face const& eList) -> void {
    Tgroup fStab=SmaGRP.Stabilizer_OnSets(eList);
    Tint OrbSizeSma=SmaGRP.size() / fStab.size();
    SizeGen += OrbSizeSma;
  };
  std::unordered_set<Face> SetFace;
  vectface CurrList(BigGRP.n_act());
  auto DoubleCosetInsertEntry_first=[&](Face const& testList) -> void {
    Face faceCan = SmaGRP.CanonicalImage(testList);
    if (SetFace.count(faceCan) > 0)
      return;
    CurrList.push_back(faceCan);
    SetFace.insert(faceCan);
    IncreaseSize(faceCan);
  };
  auto DoubleCosetInsertEntry_second=[&](Face const& testList) -> void {
    Face faceCan = SmaGRP.CanonicalImage(testList);
    if (SetFace.count(faceCan) > 0)
      return;
    SetFace.insert(faceCan);
    IncreaseSize(faceCan);
  };
  DoubleCosetInsertEntry_first(eList);
  while(true) {
    if (CurrList.size() == 0)
      break;
    Face eFace = CurrList.pop();
    for (auto const& eGen : ListGen) {
      Face eNewList=OnFace(eFace, eGen);
      DoubleCosetInsertEntry_first(eNewList);
    }
  }
  vectface ListListSet(BigGRP.n_act());
  for (auto & eFace : SetFace)
    ListListSet.push_back(eFace);
  if (SizeGen == TotalSize)
    return ListListSet;
  os << "After Iteration loop SizeGen=" << SizeGen << " TotalSize=" << TotalSize << "\n";
  std::unordered_set<Face> PartialOrbit = SetFace;
  vectface ListListSet_pop(ListListSet);
  while(true) {
    Face eFace = ListListSet_pop.pop();
    for (auto & eGen : ListGen) {
      Face eNewList=OnFace(eFace, eGen);
      if (PartialOrbit.count(eNewList) == 0) {
        PartialOrbit.insert(eNewList);
        ListListSet_pop.push_back(eNewList);
        DoubleCosetInsertEntry_second(eNewList);
        if (SizeGen == TotalSize) {
          vectface ListListFin(BigGRP.n_act());
          for (auto & eFace : SetFace)
            ListListFin.push_back(eFace);
          return ListListFin;
        }
      }
    }
  }
  os << "Likely not reachable stage\n";
  throw TerminalException{1};
}







template<typename Tgroup>
vectface OrbitSplittingListOrbit(Tgroup const& BigGRP, Tgroup const& SmaGRP, vectface eListBig, std::ostream & os)
{
  os << "|BigGRP|=" << BigGRP.size() << " |SmaGRP|=" << SmaGRP.size() << "\n";
  if (BigGRP.size() == SmaGRP.size())
    return eListBig;
  WeightMatrix<int> WMat=WeightMatrixFromPairOrbits<int,Tgroup>(SmaGRP, os);
  LocalInvInfo LocalInv=ComputeLocalInvariantStrategy(WMat, SmaGRP, "pairinv", os);
  vectface eListSma(BigGRP.n_act());
  for (auto & eSet : eListBig) {
    //    vectface ListListSet=DoubleCosetDescription_Representation(BigGRP, SmaGRP, LocalInv, eSet, os);
    vectface ListListSet=DoubleCosetDescription_Canonic(BigGRP, SmaGRP, eSet, os);
    eListSma.append(ListListSet);
  }
  os << "OrbitSplitting |eListBig|=" << eListBig.size() << " |eListSma|=" << eListSma.size() << "\n";
  return eListSma;
}



#endif
