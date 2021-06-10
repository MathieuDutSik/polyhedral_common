#ifndef INCLUDE_PERMLIB_SPECIFIC_H
#define INCLUDE_PERMLIB_SPECIFIC_H


#include "Boost_bitset.h"
#include "hash_functions.h"
#include <permlib/permlib_api.h>


int OnPoints(int const& i, permlib::Permutation const& elt)
{
  return elt.at(i);
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





std::pair<bool,permlib::Permutation> PERMLIB_TestEquivalenceGeneralNaked(int const& n, PermutationGroupPtr const& group, Face const& eList1, Face const& eList2, int const& eMethod)
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

template<typename Tint>
std::pair<bool,permlib::Permutation> PERMLIB_TestEquivalenceGeneral(int const& n, PermutationGroupPtr const& group, Tint const& grp_size, Face const& eList1, Face const& eList2)
{
  Tint MaxSize=10000;
  int eMethod = 1;
  if (grp_size < MaxSize)
    eMethod=0;
  return PERMLIB_TestEquivalenceGeneralNaked(n, group, eList1, eList2, eMethod);
}



template<typename Tint>
bool PERMLIB_TestEquivalence(int const& n, PermutationGroupPtr const& group, Tint const& grp_size, Face const& eList1, Face const& eList2)
{
  return PERMLIB_TestEquivalenceGeneral(n, group, grp_size, eList1, eList2).first;
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
  eListO=smallestSetImage(*group, eListI);
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
    return permlib::setStabilizer_classic(*group, eSet.begin(), eSet.end());
  else
    return permlib::setStabilizer_partition(*group, eSet.begin(), eSet.end());
}





template<typename Tint>
PermutationGroupPtr PERMLIB_GetStabilizer(PermutationGroupPtr const& group, Tint const& grp_size, Face const& eList)
{
  Tint MaxSize=10000;
  if (grp_size < MaxSize) {
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
  Tint_inp e_size;
  PermutationGroupPtr group;
  TheGroupFormat(int const& _n, Tint_inp const& _size, PermutationGroupPtr const& _group) : n(_n), e_size(_size), group(_group)
  {
  }
public:
  using Tint = Tint_inp;
  using Telt = permlib::Permutation;
  TheGroupFormat(std::vector<permlib::Permutation> const& ListPerm, int const& n_inp) : n(n_inp)
  {
    std::vector<permlib::Permutation::ptr> generatorList;
    for (auto & eGen : ListPerm) {
      std::vector<int> v(n);
      for (int i=0; i<n; i++)
        v[i]=eGen.at(i);
      generatorList.push_back(permlib::Permutation::ptr(new permlib::Permutation(v)));
    }
    group = construct(n, generatorList.begin(), generatorList.end());
    e_size = group->order<Tint>();
  }
  TheGroupFormat(int const& n_inp) : TheGroupFormat({}, n) 
  {
  }
  TheGroupFormat() : TheGroupFormat(0)
  {
  }
  Face CanonicalImage(Face const& eFace) const
  {
    return PERMLIB_Canonicalization(n, group, eFace);
  }
  TheGroupFormat<Tint> Stabilizer_OnSets(Face const& f) const
  {
    PermutationGroupPtr group_stab = PERMLIB_GetStabilizer(group, e_size, f);
    Tint_inp size_stab = group_stab->order<Tint_inp>();
    return TheGroupFormat(n, size_stab, group_stab);
  }
  std::pair<bool,permlib::Permutation> RepresentativeAction_OnSets(Face const& f1, Face const& f2) const
  {
    return PERMLIB_TestEquivalenceGeneral(n, group, e_size, f1, f2);
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
    return e_size;
  }
  int n_act() const
  {
    return n;
  }
  PermutationGroupPtr get_group() const // USe only by some function specific to permlib.
  {
    return group;
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





void WriteVectorInt(std::ostream &os, std::vector<int> const& OneInc)
{
  int i, siz;
  siz=OneInc.size();
  for (i=0; i<siz; i++)
    os << " " << OneInc[i];
  os << "\n";
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


template<typename Tint>
IteratorGrp GetInitialIterator(TheGroupFormat<Tint> const& eGRP)
{
  int n=eGRP.n_act();
  std::vector<MyFormTransversal> ListTrans;
  for (auto & eTrans : eGRP.get_group()->U) {
    MyFormTransversal TheTrans=GetListPermutation(eGRP.get_group(), eTrans);
    ListTrans.push_back(TheTrans);
  }
  int nbClass=ListTrans.size();
  std::vector<size_t> ListPos(nbClass, 0);
  return {n, ListPos, ListTrans};
}

template<typename Tint>
struct OrbitMinimumArr {
  int n;
  Tint GRPsize;
  std::vector<MyFormTransversal> ListTrans;
};

template<typename Tint>
OrbitMinimumArr<Tint> GetInitialMinimumArray(TheGroupFormat<Tint> const& eGRP)
{
  //  std::cerr << "GetInitialMinimumArray |GRP|=" << eGRP.size << "\n";
  IteratorGrp eIter=GetInitialIterator(eGRP);
  return {eIter.n, eGRP.size(), eIter.ListTrans};
}

template<typename Tint>
struct ResultMinimum {
  Face eMin;
  Tint OrbitSize;
};

template<typename Tint>
ResultMinimum<Tint> GetMinimumRecord(OrbitMinimumArr<Tint> const& ArrMin, Face const& eFace)
{
  int nbTrans=ArrMin.ListTrans.size();
  int n=ArrMin.n;
  std::vector<Face> ListFace(nbTrans+1,eFace);
  std::vector<size_t> ListPos(nbTrans,0);
  Face FaceMin=eFace;
  Tint nbAtt=0;
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
  Tint eOrbitSize=ArrMin.GRPsize / nbAtt;
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





template<typename Tint>
bool IsSubgroup(TheGroupFormat<Tint> const& g1, TheGroupFormat<Tint> const& g2)
{
  IteratorGrp eIter=GetInitialIterator(g1);
  for (auto & eGen : g2.group->S) {
    bool test=TestBelongingInGroup(eIter, *eGen);
    if (!test)
      return false;
  }
  return true;
}

























#endif
