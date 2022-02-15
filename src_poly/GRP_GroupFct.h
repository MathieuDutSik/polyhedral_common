#ifndef INCLUDE_TEMP_GROUP_FCT_H
#define INCLUDE_TEMP_GROUP_FCT_H

#include "Temp_common.h"
#include "Basic_file.h"
#include "Basic_string.h"
#include "hash_functions.h"
#include "COMB_Stor.h"
#include "NumberTheory.h"
#include "Boost_bitset.h"



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
	std::cerr << "Error in ReadGroup function at i=" << i << "/" << n_i << "\n";
	std::cerr << "Number of elements acted on n=" << n << " iGen=" << iGen << "/" << nbGen << "\n";
	std::cerr << "But eVal=" << eVal << " eVal_i=" << eVal_i << "\n";
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

void f_print(std::vector<int> const& V, std::string const& estr)
{
  std::cerr << estr << " =";
  for (auto & val : V)
    std::cerr << " " << val;
  std::cerr << "\n";
};


template<typename Tgroup>
std::vector<int> OrbitIntersection(Tgroup const& TheGRP, std::vector<int> const& gList)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n=TheGRP.n_act();
  std::vector<int> rList = gList;
  auto f_sum=[&]() -> Tidx {
    Tidx eSum = 0;
    for (Tidx i=0; i<n; i++)
      eSum += rList[i];
    return eSum;
  };
  auto LGen = TheGRP.GeneratorsOfGroup();
  Tidx eSum = f_sum();
  f_print(rList, "input(rList)");
  while(true) {
    for (Tidx i=0; i<n; i++) {
      if (rList[i] == 0) {
        for (auto & eGen : LGen) {
          int j = OnPoints(i, eGen);
	  rList[j]=0;
        }
      }
    }
    f_print(rList, "iter(rList)");
    Tidx eSumNew = f_sum();
    if (eSum == eSumNew)
      break;
    eSum = eSumNew;
  }
  f_print(rList, "returning(rList)");
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
  f_print(gList, "OrbitUnion(gList)");
  f_print(gListB, "OrbitUnion(gListB)");
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
  size_t eSum = rList.count();
  while(true) {
    for (Tidx i=0; i<n; i++) {
      if (rList[i] == 0) {
        for (auto & eGen : LGen) {
          Tidx j = OnPoints(i, eGen);
          rList[j]=0;
        }
      }
    }
    size_t eSumNew = rList.count();
    if (eSum == eSumNew)
      break;
    eSum = eSumNew;
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

template<typename Telt>
Telt ReduceElementAction(Telt const& eElt, Face const& eList)
{
  using Tidx = typename Telt::Tidx;
  Tidx nb=eList.count();
  if (nb == 0) {
    std::cerr << "Call of ReducedGroupAction with 0 points\n";
    throw TerminalException{1};
  }
  std::vector<Tidx> ListPositionRev(eElt.size(), -1);
  boost::dynamic_bitset<>::size_type aRow=eList.find_first();
  std::vector<Tidx> ListPosition(nb);
  for (Tidx iRow=0; iRow<nb; iRow++) {
    ListPositionRev[aRow] = iRow;
    ListPosition[iRow] = Tidx(aRow);
    aRow=eList.find_next(aRow);
  }
  std::vector<Tidx> v(nb);
  for (size_t i=0; i<nb; i++) {
    Tidx eVal1=ListPosition[i];
    Tidx eVal2=OnPoints(eVal1, eElt);
    Tidx eVal3=ListPositionRev[eVal2];
    v[i]=eVal3;
  }
  return Telt(std::move(v));
}



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
  using Telt = typename Tgroup::Telt;
  for (auto & fList : ListListSet) {
    std::optional<Telt> test=TheGRP.RepresentativeAction_OnSets(eList, fList);
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
  using Telt = typename Tgroup::Telt;
  size_t nb=ListSet.size();
  for (size_t iList=0; iList<nb; iList++)
    if (eInv == ListInv[iList]) {
      std::optional<Telt> test=TheGRP.RepresentativeAction_OnSets(eList, ListSet[iList]);
      if (test)
	return;
    }
  ListSet.push_back(eList);
  ListInv.push_back(eInv);
}

//
// Some combinatorial algorithms using only the generators of the group.
//

template<typename Telt>
std::vector<int> ComputeFullOrbitPoint(const size_t& n, const std::vector<Telt>& ListGen, int const& ePoint)
{
  IntegerSubsetStorage Vorb = VSLT_InitializeStorage(n);
  IntegerSubsetStorage Vactive = VSLT_InitializeStorage(n);
  std::vector<int> eList;
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



template<typename Telt>
vectface DecomposeOrbitPoint_Kernel(const std::vector<Telt>& LGen, Face const& eList)
{
  size_t nbPoint=eList.size();
  IntegerSubsetStorage Vlist = VSLT_InitializeStorage(nbPoint);
  vectface ListOrb(nbPoint);
  size_t len=eList.count();
  boost::dynamic_bitset<>::size_type aRow=eList.find_first();
  for (size_t i=0; i<len; i++) {
    VSLT_StoreValue(Vlist, size_t(aRow));
    aRow=eList.find_next(aRow);
  }
  while(true) {
    if (VSLT_IsEmpty(Vlist))
      break;
    size_t TheFirst = VSLT_TheFirstPosition(Vlist);
    std::vector<int> eOrb=ComputeFullOrbitPoint(nbPoint, LGen, TheFirst);
    Face vectOrb(nbPoint);
    for (auto & ePt : eOrb) {
      vectOrb[ePt]=1;
      VSLT_RemoveValue(Vlist, ePt);
    }
    ListOrb.push_back(vectOrb);
  }
  return ListOrb;
}



template<typename Tgroup>
vectface DecomposeOrbitPoint(Tgroup const& TheGRP, Face const& eList)
{
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> LGen = TheGRP.GeneratorsOfGroup();
  return DecomposeOrbitPoint_Kernel(LGen, eList);
}



template<typename Tgroup>
vectface DecomposeOrbitPoint_Full(Tgroup const& TheGRP)
{
  size_t n = TheGRP.n_act();
  Face eList(n);
  for (size_t i=0; i<n; i++)
    eList[i] = 1;
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> LGen = TheGRP.GeneratorsOfGroup();
  return DecomposeOrbitPoint_Kernel(LGen, eList);
}



template<typename Telt>
vectface DecomposeOrbitPoint_KernelFull(const size_t& n, const std::vector<Telt>& LGen)
{
  Face eList(n);
  for (size_t i=0; i<n; i++)
    eList[i] = 1;
  return DecomposeOrbitPoint_Kernel(LGen, eList);
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



template<typename Telt, typename Tobj, typename Tact>
std::vector<std::pair<Tobj,Telt>> OrbitWithRepresentative(const Telt& id, std::vector<Telt> const& ListGen, Tobj const& x, Tact act)
{
  std::vector<std::pair<Tobj,Telt>> ListObj{ {x, id} };
  std::unordered_set<Tobj> SetObj{x};
  size_t curr_pos = 0;
  while(true) {
    size_t len=ListObj.size();
    if (curr_pos == len)
      break;
    for (size_t u=curr_pos; u<len; u++) {
      for (auto & eElt : ListGen) {
        Tobj eImg = act(ListObj[u].first, eElt);
        if (SetObj.count(eImg) == 0) {
          Telt NewElt = ListObj[u].second * eElt;
          ListObj.push_back( {eImg, NewElt} );
          SetObj.insert(eImg);
        }
      }
    }
    curr_pos = len;
  }
  return ListObj;
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



template<typename Tgroup, typename F>
void OrbitSplittingSet_Kernel(vectface const& PreListTotal, Tgroup const& TheGRP, F f)
{
  using Telt = typename Tgroup::Telt;
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
    f(eSet, SingleOrbit);
  }
}



template<typename Tgroup>
vectface OrbitSplittingSet(vectface const& PreListTotal, Tgroup const& TheGRP)
{
  vectface TheReturn(TheGRP.n_act());
  auto f=[&](Face const& eSet, [[maybe_unused]] std::unordered_set<Face> const& SingleOrbit) -> void {
    TheReturn.push_back(eSet);
  };
  OrbitSplittingSet_Kernel(PreListTotal, TheGRP, f);
  return TheReturn;
}



template<typename Tgroup>
vectface OrbitSplittingSet_GetMinimalOrbit(vectface const& PreListTotal, Tgroup const& TheGRP)
{
  using Tidx = typename Tgroup::Telt::Tidx;
  Tidx len = TheGRP.n_act();
  vectface TheReturn(len);
  Face TheMin;
  bool HasMin = false;
  std::cerr << "OrbitSizes =";
  auto f=[&]([[maybe_unused]] Face const& eSet, std::unordered_set<Face> const& SingleOrbit) -> void {
    //    std::cerr << "f : begin\n";
    vectface orbit(len);
    Face minF;
    bool IsFirst = true;
    std::cerr << " " << SingleOrbit.size();
    for (auto & uSet : SingleOrbit) {
      orbit.push_back(uSet);
      if (IsFirst) {
        minF = uSet;
        IsFirst = false;
      } else {
        if (uSet < minF)
          minF = uSet;
      }
    }
    //    std::cerr << "f : orbit and minF built\n";
    // Now doing the comparison with existing data
    auto set_return=[&]() -> void {
      TheReturn = std::move(orbit);
      TheMin = minF;
    };
    if (!HasMin) {
      set_return();
      HasMin = true;
    } else {
      if (orbit.size() < TheReturn.size()) {
        set_return();
      } else {
        if (orbit.size() == TheReturn.size()) {
          if (minF < TheMin) {
            set_return();
          }
        }
      }
    }
    //    std::cerr << "f : end\n";
  };
  OrbitSplittingSet_Kernel(PreListTotal, TheGRP, f);
  std::cerr << " -- ";
  return TheReturn;
}



// Test if f1 is a subset of f2
bool is_subset(Face const& f1, Face const& f2)
{
  boost::dynamic_bitset<>::size_type pos=f1.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    if (f2[pos] == 0)
      return false;
    pos = f1.find_next(pos);
  }
  return true;
}


/*
  This is a combinatorial algorithm for finding all group elements g such that  set2 \subset g . set1
  We could probably do better with double coset decompositions.
 */
template<typename Tgroup>
std::vector<std::pair<Face, typename Tgroup::Telt>> FindContainingOrbit(Tgroup const& GRP_ext, Face const& set1, Face const& set2)
{
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tgroup stab1 = GRP_ext.Stabilizer_OnSets(set1);
  Tgroup stab2 = GRP_ext.Stabilizer_OnSets(set2);
#ifdef DEBUG_GROUP
  std::cerr << "|stab1=" << stab1.size() << " |stab2|=" << stab2.size() << "\n";
#endif
  std::vector<Telt> LGen = GRP_ext.GeneratorsOfGroup();
  Tidx n = GRP_ext.n_act();
  Telt id(n);
  auto f_act=[](const Face& x, const Telt& u) -> Face {
    return OnFace(x, u);
  };
  std::vector<std::pair<Face,Telt>> LPair = OrbitWithRepresentative(id, LGen, set1, f_act);
#ifdef DEBUG_GROUP
  std::cerr << "|LPair|=" << LPair.size() << "\n";
#endif
  std::unordered_map<Face,Telt> map_face_elt;
  vectface list_face(n);
  for (auto & ePair : LPair) {
    if (is_subset(set2, ePair.first)) {
      map_face_elt[ePair.first] = ePair.second;
      list_face.push_back(ePair.first);
    }
  }
#ifdef DEBUG_GROUP
  std::cerr << "|list_face|=" << list_face.size() << "\n";
#endif
  vectface vfRepr = OrbitSplittingSet(list_face, stab2);
  std::vector<std::pair<Face, typename Tgroup::Telt>> list_ret;
  for (auto & eRepr : vfRepr) {
    Telt eElt = map_face_elt[eRepr];
    list_ret.push_back( {eRepr, eElt} );
  }
  return list_ret;
}


#endif
