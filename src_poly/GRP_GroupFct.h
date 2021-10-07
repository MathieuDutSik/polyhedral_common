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
  int n_i;
  is >> n_i;
  is >> nbGen;
  Tidx n = Tidx(n_i);
  std::cerr << "n=" << n << " nbGen=" << nbGen << "\n";
  std::vector<Telt> ListGen;
  for (int iGen=0; iGen<nbGen; iGen++) {
    std::vector<Tidx> v(n);
    for (Tidx i=0; i<n; i++) {
      int eVal_i;
      is >> eVal_i;
      Tidx eVal = Tidx(eVal_i);
      if (eVal >= n) {
	std::cerr << "Error in ReadGroup function\n";
	std::cerr << "Number of elements acted on n=" << n << "\n";
	std::cerr << "But eVal=" << eVal << "\n";
	throw TerminalException{1};
      }
      v[i] = eVal;
    }
    ListGen.emplace_back(std::move(Telt(std::move(v))));
  }
  std::cerr << "We have read the generators\n";
  return Tgroup(ListGen, n);
}

template<typename Tgroup>
void WriteGroup(std::ostream &os, Tgroup const& TheGRP)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  int nbGen=ListGen.size();
  Tidx n = TheGRP.n_act();
  os << TheGRP.n_act() << " " << nbGen << "\n";
  for (auto & eGen : ListGen) {
    for (Tidx i=0; i<n; i++) {
      Tidx eVal=OnPoints(i, eGen);
      os << " " << int(eVal);
    }
    os << "\n";
  }
}


template<typename Tgroup>
void WriteGroupMakeUp(std::ostream &os, Tgroup const& TheGRP)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  os << "nb acting element=" << TheGRP.n_act() << "\n";
  std::vector<Telt> ListGen=TheGRP.GeneratorsOfGroup();
  int nbGen=0;
  Tidx n = TheGRP.n_act();
  for (auto & eGen : ListGen) {
    for (Tidx i=0; i<TheGRP.n_act(); i++) {
      Tidx eVal = OnPoints(i, eGen);
      os << " " << int(eVal);
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
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> ListGen=TheGRP.GeneratorsOfGroup();
  os << "local eListList, ListGen, GRP;\n";
  os << "eListList:=[\n";
  bool IsFirst=true;
  Tidx n = TheGRP.n_act();
  for (auto & eGen : ListGen) {
    if (!IsFirst)
      os << ",\n";
    IsFirst=false;
    os << "[";
    for (Tidx i=0; i<n; i++) {
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
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<int> eListReturn;
  Tidx n=TheGRP.n_act();
  std::vector<int> rList = gList;
  auto LGen = TheGRP.GeneratorsOfGroup();
  while(true) {
    Tidx eSumPrev = 0;
    for (Tidx i=0; i<n; i++)
      eSumPrev += rList[i];
    for (Tidx i=0; i<n; i++) {
      if (rList[i] == 0) {
        for (auto & eGen : LGen) {
          int j = OnPoints(i, eGen);
	  rList[j]=0;
        }
      }
    }
    std::size_t eSum=0;
    for (Tidx i=0; i<n; i++)
      eSum += rList[i];
    if (eSum == eSumPrev)
      break;
  }
  return rList;
}


template<typename Tgroup>
std::vector<int> OrbitUnion(Tgroup const& TheGRP, std::vector<int> const& gList)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n=TheGRP.n_act();
  std::vector<int> gListB(n);
  for (Tidx i=0; i<n; i++)
    gListB[i] = 1 - gList[i];
  std::vector<int> rListB = OrbitIntersection(TheGRP, gListB);
  for (Tidx i=0; i<n; i++)
    rListB[i] = 1 - rListB[i];
  return rListB;
}

template<typename Tgroup>
Face OrbitIntersection(Tgroup const& GRP, Face const& gList)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> LGen = GRP.GeneratorsOfGroup();
  Tidx n = GRP.n_act();
  Face rList = gList;
  while(true) {
    size_t eSumPrev = rList.count();
    for (Tidx i=0; i<n; i++) {
      if (rList[i] == 0) {
        for (auto & eGen : LGen) {
          Tidx j = OnPoints(i, eGen);
          rList[j]=0;
        }
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
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = GRP.n_act();
  Face gListB(n);
  for (Tidx i=0; i<n; i++)
    gListB[i] = 1 - gList[i];
  Face rListB=OrbitIntersection(GRP, gListB);
  for (Tidx i=0; i<n; i++)
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
  Tidx nb=eList.count();
  if (nb == 0) {
    std::cerr << "Call of ReducedGroupAction with 0 points\n";
    throw TerminalException{1};
  }
  std::vector<Tidx> ListPositionRev(TheGRP.n_act(), -1);
  boost::dynamic_bitset<>::size_type aRow=eList.find_first();
  std::vector<Tidx> ListPosition(nb);
  for (Tidx iRow=0; iRow<nb; iRow++) {
    ListPositionRev[aRow] = iRow;
    ListPosition[iRow] = Tidx(aRow);
    aRow=eList.find_next(aRow);
  }
  std::vector<Telt> ListGen;
  for (auto & eGen : TheGRP.GeneratorsOfGroup()) {
    std::vector<Tidx> v(nb);
    for (size_t i=0; i<nb; i++) {
      Tidx eVal1=ListPosition[i];
      Tidx eVal2=OnPoints(eVal1, eGen);
      Tidx eVal3=ListPositionRev[eVal2];
      v[i]=eVal3;
    }
    ListGen.emplace_back(std::move(Telt(std::move(v))));
  }
  return Tgroup(ListGen, nb);
}


template<typename Tgroup>
Tgroup ConjugateGroup(Tgroup const& TheGRP, typename Tgroup::Telt const& ePerm)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n=TheGRP.n_act();
  Telt ePermR=~ePerm;
  std::vector<Telt> ListGen;
  for (auto & eGen : TheGRP.GeneratorsOfGroup()) {
    std::vector<Tidx> v(n);
    for (Tidx i=0; i<n; i++) {
      Tidx eVal1=OnPoints(i, ePermR);
      Tidx eVal2=OnPoints(eVal1, eGen);
      Tidx eVal3=OnPoints(eVal2, ePerm);
      v[i]=eVal3;
    }
    ListGen.emplace_back(std::move(Telt(std::move(v))));
  }
  return Tgroup(ListGen, n);
}


//
// Some enumeration code
//

template<typename Tgroup>
void GROUP_FuncInsertInSet(Tgroup const& TheGRP, Face const& eList, vectface &ListListSet)
{
  for (auto & fList : ListListSet) {
    bool test=TheGRP.RepresentativeAction_OnSets(eList, fList).first;
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
  size_t nb=ListSet.size();
  for (size_t iList=0; iList<nb; iList++)
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
  boost::dynamic_bitset<>::size_type aRow=eList.find_first();
  for (int i=0; i<len; i++) {
    VSLT_StoreValue(Vlist, int(aRow));
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
  using Tidx=typename Telt::Tidx;
  size_t nbExt=eSet.size();
  Face fSet(nbExt);
  boost::dynamic_bitset<>::size_type pos=eSet.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    Tidx jExt=eElt.at(pos);
    fSet[jExt]=1;
    pos = eSet.find_next(pos);
  }
  return fSet;
}



template<typename Telt>
void OnFace_inplace(Face & fSet, Face const& eSet, Telt const& eElt)
{
  using Tidx=typename Telt::Tidx;
  fSet.reset();
  boost::dynamic_bitset<>::size_type pos=eSet.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    Tidx jExt=eElt.at(pos);
    fSet[jExt]=1;
    pos = eSet.find_next(pos);
  }
}



template<typename Telt>
vectface OrbitFace(const Face& f, const std::vector<Telt>& LGen)
{
  size_t len = f.size();
  std::unordered_set<Face> set_f;
  vectface vf(len);
  Face gw(len);
  std::vector<uint8_t> status;
  auto func_insert=[&](const Face& f) -> void {
    if (set_f.find(f) == set_f.end()) {
      vf.push_back(f);
      set_f.insert(f);
      status.push_back(0);
    }
  };
  func_insert(f);
  while(true) {
    bool is_finished=true;
    size_t len = vf.size();
    for (size_t i=0; i<len; i++) {
      if (status[i] == 0) {
        is_finished=false;
        status[i] = 1;
        Face fw = vf[i];
        for (auto & eGen : LGen) {
          OnFace_inplace(gw, fw, eGen);
          func_insert(gw);
        }
      }
    }
    if (is_finished)
      break;
  }
  return vf;
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
  Face fSet(TheGRP.n_act());
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
	  OnFace_inplace(fSet, gSet, eGen);
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


template<typename Tgroup, typename Tidx_value>
vectface DoubleCosetDescription_Representation(Tgroup const& BigGRP, Tgroup const& SmaGRP, WeightMatrix<true, int, Tidx_value> const& WMat, Face const& eList)
{
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  std::vector<Telt> ListGen=BigGRP.GeneratorsOfGroup();
  Tgroup TheStab=BigGRP.Stabilizer_OnSets(eList);
  Tint TotalSize=BigGRP.size() / TheStab.size();
  //
  struct Local {
    int status;
    Face eFace;
    size_t eInv;
  };
  Tint SizeGen=0;
  std::vector<Local> ListLocal;
  auto DoubleCosetInsertEntry=[&](Face const& testList) -> void {
    size_t eInv=GetLocalInvariantWeightMatrix(WMat, testList);
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
  std::cerr << "Likely not reachable stage\n";
  throw TerminalException{1};
}



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
  Face eFaceImg(BigGRP.n_act());
  while(true) {
    if (CurrList.size() == 0)
      break;
    Face eFace = CurrList.pop();
    for (auto const& eGen : ListGen) {
      OnFace_inplace(eFaceImg, eFace, eGen);
      DoubleCosetInsertEntry_first(eFaceImg);
    }
  }
  vectface ListListSet(BigGRP.n_act());
  for (auto & eFace : SetFace)
    ListListSet.push_back(eFace);
  if (SizeGen == TotalSize)
    return ListListSet;
  //  os << "After Iteration loop SizeGen=" << SizeGen << " TotalSize=" << TotalSize << "\n";
  std::unordered_set<Face> PartialOrbit = SetFace;
  while(true) {
    Face eFace = ListListSet.pop();
    for (auto & eGen : ListGen) {
      OnFace_inplace(eFaceImg, eFace, eGen);
      if (PartialOrbit.count(eFaceImg) == 0) {
        PartialOrbit.insert(eFaceImg);
        ListListSet.push_back(eFaceImg);
        DoubleCosetInsertEntry_second(eFaceImg);
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
vectface DoubleCosetDescription_Exhaustive(Tgroup const& BigGRP, Tgroup const& SmaGRP, Face const& eList, std::ostream & os)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = BigGRP.n_act();
  std::vector<Telt> LGenBig = BigGRP.GeneratorsOfGroup();
  //
  vectface vf(n);
  std::vector<uint8_t> status;
  std::unordered_set<Face> SetFace;
  auto f_insert=[&](const Face& f) -> void {
    if (SetFace.find(f) == SetFace.end()) {
      vf.push_back(f);
      status.push_back(0);
      SetFace.insert(f);
    }
  };
  size_t miss_val = std::numeric_limits<size_t>::max();
  auto get_undone=[&]() -> size_t {
    for (size_t i=0; i<status.size(); i++)
      if (status[i] == 0)
        return i;
    return miss_val;
  };
  f_insert(eList);
  Face eFaceImg(n);
  while (true) {
    size_t pos = get_undone();
    if (pos == miss_val)
      break;
    status[pos] = 1;
    Face f = vf[pos];
    for (auto& eGen : LGenBig) {
      OnFace_inplace(eFaceImg, f, eGen);
      f_insert(eFaceImg);
    }
  }
  return OrbitSplittingSet(vf, SmaGRP);
}





template<typename Tgroup>
vectface OrbitSplittingListOrbit_spec(Tgroup const& BigGRP, Tgroup const& SmaGRP, const vectface& eListBig, std::string const& method_split, std::ostream & os)
{
  std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
  os << "|BigGRP|=" << BigGRP.size() << " |SmaGRP|=" << SmaGRP.size() << "\n";
  using Tidx_value = uint16_t;
  WeightMatrix<true,int,Tidx_value> WMat;
  if (method_split == "repr") {
    WMat = WeightMatrixFromPairOrbits<Tgroup,Tidx_value>(SmaGRP);
  }
  vectface eListSma(BigGRP.n_act());
  for (auto & eSet : eListBig) {
    if (method_split == "repr") {
      vectface ListListSet=DoubleCosetDescription_Representation<Tgroup,Tidx_value>(BigGRP, SmaGRP, WMat, eSet);
      eListSma.append(ListListSet);
    }
    if (method_split == "canonic") {
      vectface ListListSet=DoubleCosetDescription_Canonic(BigGRP, SmaGRP, eSet, os);
      eListSma.append(ListListSet);
    }
    if (method_split == "exhaustive") {
      vectface ListListSet=DoubleCosetDescription_Exhaustive(BigGRP, SmaGRP, eSet, os);
      eListSma.append(ListListSet);
    }
  }
  std::chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  int elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
  os << "OrbitSplitting elapsed_seconds=" << elapsed_seconds << " |eListBig|=" << eListBig.size() << " |eListSma|=" << eListSma.size() << "\n";
  return eListSma;
}




template<typename Tgroup>
vectface OrbitSplittingListOrbit(Tgroup const& BigGRP, Tgroup const& SmaGRP, const vectface& eListBig, std::ostream & os)
{
  std::string method_split = "canonic";
  return OrbitSplittingListOrbit_spec(BigGRP, SmaGRP, eListBig, method_split, os);
}


// This is for the lifting of orbits coming with the dual description of perfect cones.
template<typename Tgroup>
void OrbitSplittingPerfectFacet(Tgroup const& BigGRP, Tgroup const& SmaGRP, const vectface& eListBig, std::ostream & os2, std::ostream& os3)
{
  std::cerr << "|BigGRP|=" << BigGRP.size() << " |SmaGRP|=" << SmaGRP.size() << "\n";
  size_t nb_orbit_big = eListBig.size();
  mpz_class nb_orbit_sma;
  size_t pos=0;
  for (auto & eSet : eListBig) {
    pos++;
    vectface ListListSet=DoubleCosetDescription_Canonic(BigGRP, SmaGRP, eSet, std::cerr);
    mpz_class orb_siz = ListListSet.size();
    nb_orbit_sma += orb_siz;
    for (auto & eFace : ListListSet) {
      mpz_class res = getsetasint<mpz_class>(eFace);
      os3 << res << "\n";
    }
    std::cerr << "iInc=" << pos << " / " << nb_orbit_big << " |ListInc2|=" << nb_orbit_sma << " |LInc|=" << orb_siz << "\n";
  }
  os2 << nb_orbit_sma << "\n";
}





#endif
