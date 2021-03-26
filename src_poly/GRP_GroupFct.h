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
    std::vector<int> v(n);
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
    std::vector<int> v(nb);
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
void GROUP_FuncInsertInSet(Tgroup const& TheGRP, Face const& eList, std::vector<Face> &ListListSet)
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
				  std::vector<Face> &ListSet,
				  std::vector<std::vector<int> > &ListInv)
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
std::vector<Face> DecomposeOrbitPoint(Tgroup const& TheGRP, Face const& eList)
{
  int nbPoint=TheGRP.n_act();
  IntegerSubsetStorage Vlist = VSLT_InitializeStorage(nbPoint);
  std::vector<Face> ListOrb;
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
      for (auto const& gObj : Additional)
	for (auto const& eGen : ListGen) {
	  Tobj fObj=TheAct(gObj, eGen);
	  if (SingleOrbit.find(fObj) == SingleOrbit.end() && Additional.find(fObj) == Additional.end()) {
	    if (NewElts.find(fObj) == NewElts.end()) {
	      if (ListTotal.find(fObj) != ListTotal.end()) {
		NewElts.insert(fObj);
		ListTotal.erase(fObj);
	      }
#ifdef DEBUG
	      else {
		std::cerr << "Orbit do not match, PANIC!!!\n";
		throw TerminalException{1};
	      }
#endif
	    }
	  }
	}
      for (auto & uSet : Additional)
	SingleOrbit.insert(uSet);
      if (NewElts.size() == 0)
	break;
      Additional=NewElts;
    }
  }
  return TheReturn;
}



template<typename Telt>
Face OnFace(Face const& eSet, Telt const& eElt)
{
  int nbExt=eSet.size();
  Face fSet(nbExt);
  for (int iExt=0; iExt<nbExt; iExt++)
    if (eSet[iExt] == 1) {
      int jExt=eElt.at(iExt);
      fSet[jExt]=1;
    }
  return fSet;
}



template<typename Tgroup>
std::vector<Face> OrbitSplittingSet(std::vector<Face> const& PreListTotal,
				    Tgroup const& TheGRP)
{
  using Telt = typename Tgroup::Telt;
  std::vector<Face> TheReturn;
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
      for (auto const& gSet : Additional)
	for (auto const& eGen : ListGen) {
	  Face fSet=OnFace(gSet, eGen);
	  if (SingleOrbit.find(fSet) == SingleOrbit.end() && Additional.find(fSet) == Additional.end()) {
	    if (NewElts.find(fSet) == NewElts.end()) {
	      if (ListTotal.find(fSet) != ListTotal.end()) {
		NewElts.insert(fSet);
		ListTotal.erase(fSet);
	      }
#ifdef DEBUG
	      else {
		std::cerr << "Orbit do not matched, PANIC!!!\n";
		throw TerminalException{1};
	      }
#endif
	    }
	  }
	}
      for (auto & uSet : Additional)
	SingleOrbit.insert(uSet);
      if (NewElts.size() == 0)
	break;
      Additional=NewElts;
    }
  }
  return TheReturn;
}


//
// Double Coset Computation
//


template<typename Tgroup>
std::vector<Face> DoubleCosetDescription(Tgroup const& BigGRP,
					 Tgroup const& SmaGRP,
					 LocalInvInfo const& LocalInv,
					 Face const& eList, std::ostream & os)
{
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  os << "Beginning of DoubleCosetDescription\n";
  std::vector<Telt> ListGen=BigGRP.GeneratorsOfGroup();
  Tgroup TheStab=BigGRP.Stabilizer_OnSets(eList);
  os << "BigGRP.size=" << BigGRP.size() << " TheStab.size=" << TheStab.size() << "\n";
  Tint TotalSize=BigGRP.size() / TheStab.size();
  os << "TotalSize=" << TotalSize << "\n";
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
    os << "SmaGRP.size=" << SmaGRP.size() << " fStab.size=" << fStab.size() << "\n";
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
  os << "After Iteration loop SizeGen=" << SizeGen << " TotalSize=" << TotalSize << "\n";
  std::vector<Face> ListListSet;
  for (auto & eRec : ListLocal)
    ListListSet.push_back(eRec.eFace);
  if (SizeGen == TotalSize)
    return ListListSet;
  std::vector<Face> PartialOrbit=ListListSet;
  auto IsPresent=[&](Face const& testList) -> bool {
    for (auto & fList : PartialOrbit)
      if (fList == testList)
	return true;
    return false;
  };
  while(true) {
    for (auto & eGen : ListGen) {
      size_t len=PartialOrbit.size();
      for (size_t i=0; i<len; i++) {
	Face eNewList=OnFace(PartialOrbit[i], eGen);
	if (!IsPresent(eNewList)) {
	  PartialOrbit.push_back(eNewList);
	  DoubleCosetInsertEntry(eNewList);
	  if (SizeGen == TotalSize) {
	    std::vector<Face> ListListFin;
	    for (auto & eRec : ListLocal)
	      ListListFin.push_back(eRec.eFace);
	    return ListListFin;
	  }
	}
      }
    }
  }
  os << "Likely not reachable stage\n";
  throw TerminalException{1};
}



template<typename Tgroup>
std::vector<Face> OrbitSplittingListOrbit(Tgroup const& BigGRP, Tgroup const& SmaGRP, std::vector<Face> eListBig, std::ostream & os)
{
  os << "|BigGRP|=" << BigGRP.size() << " |SmaGRP|=" << SmaGRP.size() << "\n";
  if (BigGRP.size() == SmaGRP.size())
    return eListBig;
  {
    std::ofstream os1("ORBSPLIT_BigGRP");
    std::ofstream os2("ORBSPLIT_BigGRP.gap");
    std::ofstream os3("ORBSPLIT_SmaGRP");
    std::ofstream os4("ORBSPLIT_SmaGRP.gap");
    std::ofstream os5("ORBSPLIT_ListBig");
    std::ofstream os6("ORBSPLIT_ListBig.gap");
    WriteGroup      (os1, BigGRP);
    WriteGroupGAP   (os2, BigGRP);
    WriteGroup      (os3, SmaGRP);
    WriteGroupGAP   (os4, SmaGRP);
    WriteListFace   (os5, eListBig);
    WriteListFaceGAP(os6, eListBig);
  }
  WeightMatrix<int,int> WMat=WeightMatrixFromPairOrbits<int,int,Tgroup>(SmaGRP, os);
  LocalInvInfo LocalInv=ComputeLocalInvariantStrategy(WMat, SmaGRP, "pairinv", os);
  os << "We do the algorithm\n";
  std::vector<Face> eListSma;
  size_t iter=0;
  for (auto & eSet : eListBig) {
    os << "iter=" << iter << " Before DoubleCosetDescription\n";
    std::vector<Face> ListListSet=DoubleCosetDescription(BigGRP, SmaGRP, LocalInv, eSet, os);
    os << "      |ListListSet|=" << ListListSet.size() << "\n";
    for (auto & eCos : ListListSet)
      eListSma.push_back(eCos);
    os << "      |eListSma|=" << eListSma.size() << "\n";
    iter++;
  }
  return eListSma;
}



#endif
