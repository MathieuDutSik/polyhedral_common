#ifndef INCLUDE_TEMP_GROUP_FCT_H
#define INCLUDE_TEMP_GROUP_FCT_H

#include "Temp_common.h"
#include "Basic_file.h"
#include "Basic_string.h"
#include "hash_functions.h"
#include "COMB_Stor.h"
#include "NumberTheory.h"
#include "Boost_bitset.h"
#include <permlib/permlib_api.h>


//
// Template general code for input output of groups
//

template<typename Tgroup>
Tgroup ReadGroup(std::istream &is)
{
  using Telt = Tgroup::Telt;
  if (!is.good()) {
    std::cerr << "ReadGroup operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  int nbGen;
  permlib::dom_int n;
  is >> n;
  is >> nbGen;
  std::cerr << "n=" << n << " nbGen=" << nbGen << "\n";
  std::vector<Telt> ListGen;
  for (int iGen=0; iGen<nbGen; iGen++) {
    std::vector<permlib::dom_int> v(n);
    for (std::size_t i=0; i<n; i++) {
      permlib::dom_int eVal;
      is >> eVal;
      if (eVal > n-1) {
	std::cerr << "Error in ReadGroup function\n";
	std::cerr << "Number of elements acted on n=" << n << "\n";
	std::cerr << "But eVal=" << eVal << "\n";
	throw TerminalException{1};
      }
      v[i]=eVal;
    }
    ListGen.push_back(permlib::Permutation(v));
  }
  return Tgroup(ListGen, n);
}

template<typename Tgroup>
void WriteGroup(std::ostream &os, Tgroup const& TheGRP)
{
  using Telt = Tgroup::Telt;
  using Tint = Tgroup::Tint;
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  int nbGen=ListGen.size();
  os << TheGRP.n << " " << nbGen << "\n";
  for (auto & eGen : ListGen) {
    for (std::size_t i=0; i<TheGRP.n; i++) {
      int eVal=OnPoints(i, eGen);
      os << " " << eVal;
    }
    os << "\n";
  }
}


template<typename Tgroup>
void WriteGroupMakeUp(std::ostream &os, Tgroup const& TheGRP)
{
  using Telt = Tgroup::Telt;
  os << "nb acting element=" << TheGRP.n_act() << "\n";
  std::vector<Telt> ListGen=TheGRP.GeneratorsOfGroup();
  int nbGen=0;
  for (auto & eGen : ListGen) {
    for (std::size_t i=0; i<TheGRP.n; i++) {
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
  using Telt = Tgroup::Telt;
  std::vector<Telt> ListGen=TheGRP.GeneratorsOfGroup();
  os << "local eListList, ListGen, GRP;\n";
  os << "eListList:=[\n";
  bool IsFirst=true;
  for (auto & eGen : ListGen) {
    if (!IsFirst)
      os << ",\n";
    IsFirst=false;
    os << "[";
    for (std::size_t i=0; i<TheGRP.n; i++) {
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
    for (std::size_t i=0; i<n; i++)
      eSumPrev += rList[i];
    for (auto & eGen : LGen) {
      for (std::size_t i=0; i<n; i++) {
	int j = OnPoints(i, eGen);
	if (rList[i] == 0)
	  rList[j]=0;
      }
    }
    std::size_t eSum=0;
    for (std::size_t i=0; i<n; i++)
      eSum += rList[i];
    if (eSum == eSumPrev)
      break;
  }
  return rList;
}


template<typename Tgroup>
std::vector<int> OrbitUnion(Tgroup const& TheGRP, std::vector<int> const& gList)
{
  int n=TheGRP.n;
  std::vector<int> gListB(n);
  for (int i=0; i<n; i++)
    gListB[i] = 1 - gList[i];
  std::vector<int> rListB=OrbitIntersection(TheGRP, gListB);
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
  using Telt = Tgroup::Telt;
  int nb=eList.count();
  if (nb == 0) {
    std::cerr << "Call of ReducedGroupAction with 0 points\n";
    throw TerminalException{1};
  }
  std::vector<int> ListPositionRev(TheGRP.n, -1);
  int aRow=eList.find_first();
  std::vector<int> ListPosition(nb,-1);
  for (int iRow=0; iRow<nb; iRow++) {
    ListPositionRev[aRow]=iRow;
    ListPosition[iRow]=aRow;
    aRow=eList.find_next(aRow);
  }
  std::vector<Telt> ListGen;
  for (auto & eGen : TheGRP.GeneratorsOfGroup()) {
    std::vector<permlib::dom_int> v(nb);
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
  using Telt = Tgroup::Telt;
  int n=TheGRP.n;
  Telt ePermR=~ePerm;
  std::vector<Telt> ListGen;
  for (auto & eGen : TheGRP.GeneratorsOfGroup()) {
    std::vector<permlib::dom_int> v(n);
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
// Specific implementation with permlib
//




typedef std::shared_ptr<permlib::PermutationGroup> PermutationGroupPtr;
typedef boost::dynamic_bitset<> DsetList;


/*
std::vector<permlib::dom_int> GetBaseGroup(TheGroupFormat const& eGRP)
{
  std::vector<permlib::dom_int> eList;
  for (auto & eTrans : eGRP.group->U) {
    permlib::dom_int eElt=eTrans.element();
    eList.push_back(eElt);
  }
  return eList;
}
*/





std::set<int> GetSetFrom_DB(Face const& eList)
{
  int nb=eList.count();
  std::set<int> eSet;
  int aRow=eList.find_first();
  for (int i=0; i<nb; i++) {
    eSet.insert(aRow);
    aRow=eList.find_next(aRow);
  }
  return eSet;
}

Face Face_EquivSimplification(Face const& eFace)
{
  int siz=eFace.size();
  int cnt=eFace.count();
  if (2*cnt > siz) {
    Face retFace(siz);
    for (int i=0; i<siz; i++)
      retFace[i]=1 - eFace[i];
    return retFace;
  }
  else {
    return eFace;
  }
}





std::pair<bool,permlib::Permutation> PERMLIB_TestEquivalenceGeneralNaked(int const& n, PermutationGroupPtr Tconst& group, Face const& eList1, Face const& eList2, int const& eMethod)
{
  permlib::Permutation::ptr mappingElement;
  int nb1=eList1.count();
  int nb2=eList2.count();
  if (nb1 != nb2)
    return {false, {}};
  //
  Face NewList1=Face_EquivSimplification(eList1);
  Face NewList2=Face_EquivSimplification(eList2);
  std::set<int> eSet1=GetSetFrom_DB(NewList1);
  std::set<int> eSet2=GetSetFrom_DB(NewList2);
  if (eMethod == 0)
    mappingElement = permlib::setImage_classic(*group, eSet1.begin(), eSet1.end(), eSet2.begin(), eSet2.end());
  if (eMethod == 1)
    mappingElement = permlib::setImage_partition(*group, eSet1.begin(), eSet1.end(), eSet2.begin(), eSet2.end());
  if (mappingElement) {
    return {true, *mappingElement};
  }
  return {false, {}};
}

std::pair<bool,permlib::Permutation> PERMLIB_TestEquivalenceGeneral(int const& n, PermutationGroupPtr const& group, Face const& eList1, Face const& eList2)
{
  mpz_class MaxSize=10000;
  int eMethod = 1;
  if (TheGRP.size < MaxSize)
    eMethod=0;
  return PERMLIB_TestEquivalenceGeneralNaked(n, group, eList1, eList2, eMethod);
}



bool PERMLIB_TestEquivalence(int const& n, PermutationGroupPtr const& group, Face const& eList1, Face const& eList2)
{
  std::pair<bool,permlib::Permutation> eRes=PERMLIB_TestEquivalenceGeneral(n, group, eList1, eList2);
  return eRes.first;
}







Face PERMLIB_Canonicalization(int const& n, PermutationGroupPtr const& group, Face const& eList)
{
  DsetList eListI(n), eListO(n);
  int siz=eList.count();
  int aRow=eList.find_first();
  for (int i=0; i<siz; i++) {
    eListI[aRow]=1;
    aRow=eList.find_next(aRow);
  }
  eListO=smallestSetImage(n, *group, eListI);
  Face TheRet;
  for (int i=0; i<n; i++)
    if (eListO[i] == 1)
      TheRet[i]=1;
#ifdef DEBUG_GROUP
  bool test=PERMLIB_TestEquivalence(n, *group, eList, TheRet);
  if (!test) {
    std::cerr << "We have major debugging to do\n";
    throw TerminalException{1};
  }
#endif
  return TheRet;
}



PermutationGroupPtr PERMLIB_GetStabilizer_general(PermutationGroupPtr const& group, Face const& eList, int const& opt)
{
  Face NewList=Face_EquivSimplification(eList);
  std::set<int> eSet=GetSetFrom_DB(NewList);
  if (opt == 0)
    return setStabilizer_classic(*TheGRP.group, eSet.begin(), eSet.end());
  if (opt == 1)
    return setStabilizer_partition(*TheGRP.group, eSet.begin(), eSet.end());
}





PermutationGroupPtr PERMLIB_GetStabilizer(PermutationGroupPtr const& group, Face const& eList)
{
  mpz_class MaxSize=10000;
  if (TheGRP.size < MaxSize) {
    return PERMLIB_GetStabilizer_general(group, eList, 0);
  }
  else {
    return PERMLIB_GetStabilizer_general(group, eList, 1);
  }
}





template<typename Tint_inp>
struct TheGroupFormat {
private:
  int n;
  Tint_inp size;
  PermutationGroupPtr group;
  TheGroupFormat(int const& _n, Tint_inp const& _size, PermutationGroupPtr const& _group) : n(_n), size(_size), group(_group)
  {
  }
public:
  using Tint = Tint_inp;
  using Telt = permlib::Permutation;
  TheGroupFormat(std::vector<permlib::Permutation> const& ListPerm, int const& n_inp)
  {
    std::vector<permlib::Permutation::ptr> generatorList;
    for (auto & eGen : ListPerm) {
      std::vector<permlib::dom_int> v(n);
      for (std::size_t i=0; i<n; i++)
        v[i]=eGen.at(i);
      generatorList.push_back(permlib::Permutation::ptr(new permlib::Permutation(v)));
    }
    n = n_inp;
    group = construct(n, generatorList.begin(), generatorList.end());
    size = TheGroupRet.group->order<mpz_class>();
  }
  Face CanonicalImage(Face const& eFace) const
  {
    return PERMLIB_Canonicalization(n, group, eFace);
  }
  TheGroupFormat Stabilizer_OnSets(Face const& f) const
  {
    PermutationGroupPtr group_stab = PERMLIB_GetStabilizer(group, f);
    Tint_inp size_stab = group_stab->order<Tint_inp>();
    return TheGroupFormat(n, size_stab, group_stab);
  }
  std::pair<bool,permlib::Permutation> RepresentativeAction_OnSets(Face const& f1, Face const& f2) const
  {
    return PERMLIB_TestEquivalenceGeneral(n, group, f1, f2);
  }
  std::vector<Telt> GeneratorsOfGroup() const
  {
    // copy operation, but that is what it is.
    std::vector<Telt> LGen;
    for (auto & eval : group->S) {
      LGen.push_back(*eval);
    }
    return LGen;
  }
  Tint size() const
  {
    return size;
  }
  int n_act() const
  {
    return n;
  }
};


namespace std {
  template<>
  struct hash<permlib::Permutation>
  {
    std::size_t operator()(permlib::Permutation const& eElt) const
    {
      size_t len=eElt.size();
      std::vector<permlib::dom_int> V(len);
      for (size_t i=0; i<len; i++)
        V[i] = eElt.at(i);
      return std::hash<std::vector<permlib::dom_int>>()(V);
    }
  };
}


Face eEltImage(Face const& eSet, permlib::Permutation const& eElt)
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





void WriteVectorInt(std::ostream &os, std::vector<int> const& OneInc)
{
  int i, siz;
  siz=OneInc.size();
  for (i=0; i<siz; i++)
    os << " " << OneInc[i];
  os << "\n";
}


void GROUP_FuncInsertInSet_UseInv(TheGroupFormat const& TheGRP,
				  Face const& eList,
				  std::vector<int> const& eInv, 
				  std::vector<Face> &ListSet,
				  std::vector<std::vector<int> > &ListInv)
{
  int nb=ListSet.size();
  for (int iList=0; iList<nb; iList++)
    if (eInv == ListInv[iList]) {
      bool test=TestEquivalence(TheGRP, eList, ListSet[iList]);
      if (test)
	return;
    }
  ListSet.push_back(eList);
  ListInv.push_back(eInv);
}





void GROUP_FuncInsertInSet(TheGroupFormat const& TheGRP, Face const& eList, std::vector<Face> &ListListSet)
{
  int iList, nb;
  nb=ListListSet.size();
  for (iList=0; iList<nb; iList++) {
    bool test=TestEquivalence(TheGRP, eList, ListListSet[iList]);
    if (test)
      return;
  }
  ListListSet.push_back(eList);
}





























std::list<int> ComputeFullOrbitPoint(TheGroupFormat const& TheGRP, int const& ePoint)
{
  IntegerSubsetStorage *Vorb;
  IntegerSubsetStorage *Vactive;
  Vorb = new IntegerSubsetStorage;
  Vactive = new IntegerSubsetStorage;
  std::list<int> eList;
  std::list<permlib::Permutation::ptr> ListGen;
  int TheFirst;
  permlib::dom_int n=TheGRP.n;
  ListGen=TheGRP.group->S;
  VSLT_InitializeStorage(Vorb, n);
  VSLT_InitializeStorage(Vactive, n);
  VSLT_StoreValue(Vactive, ePoint);
  while(true) {
    if (VSLT_IsEmpty(Vactive) == 1)
      break;
    TheFirst=VSLT_TheFirstPosition(Vactive);
    VSLT_RemoveValue(Vactive, TheFirst);
    VSLT_StoreValue(Vorb, TheFirst);
    eList.push_back(TheFirst);
    for (auto & eGen : ListGen) {
      permlib::dom_int NewPt=eGen->at(TheFirst);
      if (!VSLT_IsItInSubset(Vorb, NewPt) && !VSLT_IsItInSubset(Vactive, NewPt))
	VSLT_StoreValue(Vactive, NewPt);
    }
  }
  VSLT_FreeStorage(Vorb);
  VSLT_FreeStorage(Vactive);
  delete Vorb;
  delete Vactive;
  return eList;
}

std::vector<Face> DecomposeOrbitPoint(TheGroupFormat const& TheGRP, Face const& eList)
{
  IntegerSubsetStorage *Vlist;
  Vlist = new IntegerSubsetStorage;
  std::vector<Face> ListOrb;
  int len, i, TheFirst;
  len=eList.count();
  int nbPoint=TheGRP.n;
  VSLT_InitializeStorage(Vlist, nbPoint);
  int aRow=eList.find_first();
  for (i=0; i<len; i++) {
    VSLT_StoreValue(Vlist, aRow);
    aRow=eList.find_next(aRow);
  }
  while(true) {
    if (VSLT_IsEmpty(Vlist) == 1)
      break;
    TheFirst=VSLT_TheFirstPosition(Vlist);
    std::list<int> eOrb=ComputeFullOrbitPoint(TheGRP, TheFirst);
    Face vectOrb(nbPoint);
    for (auto & ePt : eOrb) {
      vectOrb[ePt]=1;
      VSLT_RemoveValue(Vlist, ePt);
    }
    ListOrb.push_back(vectOrb);
  }
  VSLT_FreeStorage(Vlist);
  delete Vlist;
  return ListOrb;
}

struct PairEltPerm {
  permlib::dom_int eElt;
  permlib::Permutation ePerm;
};


struct MyFormTransversal {
  permlib::dom_int eElt;
  std::vector<PairEltPerm> ListPair;
};


permlib::Permutation IdentityPermutation(int const& n)
{
  std::vector<permlib::dom_int> v(n);
  for (int i=0; i<n; i++)
    v[i]=i;
  return permlib::Permutation(v);
}


MyFormTransversal GetListPermutation(PermutationGroupPtr TheGRP,
				     permlib::SchreierTreeTransversal<permlib::Permutation> const& eTrans)
{
  permlib::dom_int n=TheGRP->n;
  permlib::dom_int eElt=eTrans.element();
  std::unordered_set<permlib::Permutation, std::hash<permlib::Permutation>> ListPermWork;
  for (std::shared_ptr<permlib::Permutation> & p : eTrans.GetMtransversal() ) {
    if (p) {
      permlib::Permutation ePerm=*p;
      ListPermWork.insert(ePerm);
    }
  }
  std::unordered_set<permlib::dom_int> SetOrbit;
  std::vector<PairEltPerm> ListPair;
  std::vector<bool> StatusDone;
  std::function<void(permlib::dom_int,permlib::Permutation)> fInsert=[&](permlib::dom_int const& eVal, permlib::Permutation const& ePerm) -> void {
    std::unordered_set<permlib::dom_int>::iterator iter=SetOrbit.find(eVal);
    if (iter == SetOrbit.end()) {
      SetOrbit.insert(eVal);
      PairEltPerm ePair{eVal, ePerm};
      ListPair.push_back(ePair);
      StatusDone.push_back(false);
    }
  };
  permlib::Permutation ePerm=IdentityPermutation(n);
  fInsert(eElt, ePerm);
  while(true) {
    bool IsFinished=true;
    int len=ListPair.size();
    for (int i=0; i<len; i++)
      if (!StatusDone[i]) {
	StatusDone[i]=true;
	IsFinished=false;
	for (auto & fPerm : ListPermWork) {
	  permlib::dom_int fVal=fPerm.at(ListPair[i].eElt);
	  permlib::Permutation eProd=ListPair[i].ePerm*fPerm;
	  permlib::dom_int gVal=eProd.at(eElt);
	  if (gVal != fVal) {
	    std::cerr << "Inconsistency here on the permutation product\n";
	    throw TerminalException{1};
	  }
	  fInsert(fVal, eProd);
	}
      }
    if (IsFinished)
      break;
  }
  return {eElt, ListPair};
}

struct IteratorGrp {
  int n;
  std::vector<size_t> ListPos;
  std::vector<MyFormTransversal> ListTrans;
};


IteratorGrp GetInitialIterator(TheGroupFormat const& eGRP)
{
  int n=eGRP.n;
  std::vector<MyFormTransversal> ListTrans;
  for (auto & eTrans : eGRP.group->U) {
    MyFormTransversal TheTrans=GetListPermutation(eGRP.group, eTrans);
    ListTrans.push_back(TheTrans);
  }
  int nbClass=ListTrans.size();
  std::vector<size_t> ListPos(nbClass, 0);
  return {n, ListPos, ListTrans};
}

struct OrbitMinimumArr {
  int n;
  mpz_class GRPsize;
  std::vector<MyFormTransversal> ListTrans;
};

OrbitMinimumArr GetInitialMinimumArray(TheGroupFormat const& eGRP)
{
  //  std::cerr << "GetInitialMinimumArray |GRP|=" << eGRP.size << "\n";
  IteratorGrp eIter=GetInitialIterator(eGRP);
  return {eIter.n, eGRP.size, eIter.ListTrans};
}

struct ResultMinimum {
  Face eMin;
  mpz_class OrbitSize;
};

ResultMinimum GetMinimumRecord(OrbitMinimumArr const& ArrMin, Face const& eFace)
{
  int nbTrans=ArrMin.ListTrans.size();
  int n=ArrMin.n;
  std::vector<Face> ListFace(nbTrans+1,eFace);
  std::vector<size_t> ListPos(nbTrans,0);
  Face FaceMin=eFace;
  mpz_class nbAtt=0;
  //  std::cerr << "|GRP|=" << ArrMin.GRPsize << "\n";
  auto Increment=[&]() -> int {
    for (int i=0; i<nbTrans; i++)
      if (ListPos[i] < ArrMin.ListTrans[i].ListPair.size() -1) {
	ListPos[i]++;
	for (int j=0; j<i; j++)
	  ListPos[j]=0;
	for (int j=i; j>=0; j--) {
	  int ePos=ListPos[j];
	  for (int iCol=0; iCol<n; iCol++) {
	    int jCol=ArrMin.ListTrans[j].ListPair[ePos].ePerm.at(iCol);
	    ListFace[j][jCol]=ListFace[j+1][iCol];
	  }
	}
	return 0;
      }
    return -1;
  };
  while(true) {
    if (ListFace[0] < FaceMin) {
      FaceMin=ListFace[0];
      nbAtt=0;
    }
    if (ListFace[0] == FaceMin) 
      nbAtt++;
    int test=Increment();
    if (test == -1)
      break;
  }
  //  std::cerr << "nbAtt=" << nbAtt << "\n";
  mpz_class eOrbitSize=ArrMin.GRPsize / nbAtt;
  return {FaceMin, eOrbitSize};
}






permlib::Permutation GetPermutation(IteratorGrp const& eIter)
{
  int n=eIter.n;
  permlib::Permutation ePerm=IdentityPermutation(n);
  int nbClass=eIter.ListPos.size();
  for (int i=0; i<nbClass; i++) {
    size_t ePos=eIter.ListPos[i];
    ePerm = eIter.ListTrans[i].ListPair[ePos].ePerm*ePerm;
  }
  return ePerm;
}

bool IsFinalIterator(IteratorGrp const& eIter)
{
  int nbClass=eIter.ListPos.size();
  for (int i=0; i<nbClass; i++) {
    size_t siz=eIter.ListTrans[i].ListPair.size();
    if (eIter.ListPos[i] < siz-1)
      return false;
  }
  return true;
}

int IteratorIncrease(IteratorGrp & eIter)
{
  int nbClass=eIter.ListPos.size();
  for (int i=0; i<nbClass; i++)
    if (eIter.ListPos[i] < eIter.ListTrans[i].ListPair.size() -1) {
      eIter.ListPos[i]++;
      for (int j=0; j<i; j++)
	eIter.ListPos[j]=0;
      return 0;
    }
  return -1;
}


bool TestBelongingInGroup(IteratorGrp const& eIter, permlib::Permutation const& ePerm)
{
  permlib::Permutation eWork=ePerm;
  int nbClass=eIter.ListTrans.size();
  std::function<bool(permlib::dom_int,int)> fUpdate=[&](permlib::dom_int const& eElt,int const& i) -> bool {
    permlib::dom_int eImg=eWork.at(eElt);
    int len=eIter.ListTrans[i].ListPair.size();
    for (int j=0; j<len; j++)
      if (eIter.ListTrans[i].ListPair[j].eElt == eImg) {
	eWork = eWork * (~eIter.ListTrans[i].ListPair[j].ePerm);
	return true;
      }
    return false;
  };
  for (int i=0; i<nbClass; i++) {
    permlib::dom_int eElt=eIter.ListTrans[i].eElt;
    bool test1=fUpdate(eElt, i);
    if (!test1)
      return false;
  }
  return eWork.isIdentity();
}




bool IsSubgroup(TheGroupFormat const& g1, TheGroupFormat const& g2)
{
  IteratorGrp eIter=GetInitialIterator(g1);
  for (auto & eGen : g2.group->S) {
    bool test=TestBelongingInGroup(eIter, *eGen);
    if (!test)
      return false;
  }
  return true;
}



template<typename T>
permlib::Permutation SortingPerm(std::vector<T> const & ListV)
{
  struct PairData {
    std::size_t i;
    T x;
  };
  std::size_t len=ListV.size();
  //  std::cerr << "len=" << len << "\n";
  std::vector<PairData> ListPair(len);
  for (std::size_t i=0; i<len; i++) {
    PairData ePair{i, ListV[i]};
    ListPair[i]=ePair;
  }
  sort(ListPair.begin(), ListPair.end(),
       [](PairData const & x1, PairData const& x2) -> bool {
	 if (x1.x < x2.x)
	   return true;
	 if (x2.x < x1.x)
	   return false;
	 return x1.i< x2.i;
       });
  std::vector<permlib::dom_int> v(len);
  for (std::size_t i=0; i<len; i++)
    v[i]=ListPair[i].i;
  return permlib::Permutation(v);
}





std::vector<int> PermutationOrbit(permlib::Permutation const& ePerm)
{
  //  std::cerr << "  Beginning of PermutationOrbit\n";
  int siz=ePerm.size();
  //  for (int i=0; i<siz; i++) {
  //    std::cerr << "i=" << i << " img=" << ePerm.at(i) << "\n";
  //  }
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
    //    std::cerr << "iPoint=" << iPoint << " idxOrbit=" << idxOrbit << "\n";
    if (iPoint == -1)
      break;
    int iPointWork=iPoint;
    while(true) {
      StatusOrbit[iPointWork]=idxOrbit;
      iPointWork=ePerm.at(iPointWork);
      //      std::cerr << "  iPointWork=" << iPointWork << "\n";
      if (iPointWork == iPoint)
	break;
    }
    idxOrbit++;
  }
  //  std::cerr << "  End of PermutationOrbit\n";
  return StatusOrbit;
}


void WritePermutationGAP(std::ostream&os, permlib::Permutation const& ePerm)
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
      TheList[i]=ePerm.at(TheList[i-1]);
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







std::vector<Face> OrbitSplittingSet(std::vector<Face> const& PreListTotal,
				    TheGroupFormat const& TheGRP)
{
  std::vector<Face> TheReturn;
  std::unordered_set<Face> ListTotal;
  for (auto eFace : PreListTotal)
    ListTotal.insert(eFace);
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
	for (auto const& eGen : TheGRP.group->S) {
	  Face fSet=eEltImage(gSet, *eGen);
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



std::vector<Face> DoubleCosetDescription(TheGroupFormat const& BigGRP,
					 TheGroupFormat const& SmaGRP,
					 LocalInvInfo const& LocalInv,
					 Face const& eList, std::ostream & os)
{
  os << "Beginning of DoubleCosetDescription\n";
  std::list<permlib::Permutation::ptr> ListGen=BigGRP.group->S;
  TheGroupFormat TheStab=GetStabilizer(BigGRP, eList);
  os << "BigGRP.size=" << BigGRP.size << " TheStab.size=" << TheStab.size << "\n";
  mpz_class TotalSize=BigGRP.size / TheStab.size;
  os << "TotalSize=" << TotalSize << "\n";
  //
  struct Local {
    int status;
    Face eFace;
    std::vector<int> eInv;
  };
  mpz_class SizeGen=0;
  std::vector<Local> ListLocal;
  auto DoubleCosetInsertEntry=[&](Face const& testList) -> void {
    std::vector<int> eInv=GetLocalInvariantWeightMatrix_Enhanced<int>(LocalInv, testList);
    for (auto const& fLocal : ListLocal) {
      bool testCL=TestEquivalenceGeneralNaked(SmaGRP, fLocal.eFace, testList, 0).TheReply;
      bool testPA=TestEquivalenceGeneralNaked(SmaGRP, fLocal.eFace, testList, 1).TheReply;
      bool test=TestEquivalence(SmaGRP, fLocal.eFace, testList);
      os << "fLocal.eFace=\n";
      WriteFace(os, fLocal.eFace);
      os << "    testList=\n";
      WriteFace(os, testList);
      os << "fLocal.eFace=" << fLocal.eFace << "\n";
      os << "    testList=" << testList << "  testCL=" << testCL << " testPA=" << testPA << "\n";
      if (test)
	return;
    }
    ListLocal.push_back({0,testList,eInv});
    TheGroupFormat fStab=GetStabilizer(SmaGRP, testList);
    os << "SmaGRP.size=" << SmaGRP.size << " fStab.size=" << fStab.size << "\n";
    mpz_class OrbSizeSma=SmaGRP.size / fStab.size;
    SizeGen += OrbSizeSma;
    os << "Now SizeGen=" << SizeGen << " OrbSizeSma=" << OrbSizeSma << " |ListLocal|=" << ListLocal.size() << "\n";
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
	  Face eNewList=eEltImage(eFace, *eGen);
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
	Face eNewList=eEltImage(PartialOrbit[i], *eGen);
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



std::vector<Face> OrbitSplittingListOrbit(TheGroupFormat const& BigGRP, TheGroupFormat const& SmaGRP, std::vector<Face> eListBig, std::ostream & os)
{
  os << "|BigGRP|=" << BigGRP.size << " |SmaGRP|=" << SmaGRP.size << "\n";
  if (BigGRP.size == SmaGRP.size)
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
  WeightMatrix<int,int> WMat=WeightMatrixFromPairOrbits<int,int>(SmaGRP, os);
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
