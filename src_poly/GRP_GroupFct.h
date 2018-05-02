#ifndef TEMP_GROUP_FUNCTIONS
#define TEMP_GROUP_FUNCTIONS

#include "Temp_common.h"
#include "Basic_file.h"
#include "Basic_string.h"
#include "COMB_Stor.h"
#include "NumberTheory.h"
#include "Face_bitset.h"
#include <permlib/permlib_api.h>

typedef std::shared_ptr<permlib::PermutationGroup> PermutationGroupPtr;
typedef boost::dynamic_bitset<> DsetList;

struct TheGroupFormat {
  permlib::dom_int n;
  mpz_class size;
  PermutationGroupPtr group;
};

Face eEltImage(Face const& eSet, permlib::Permutation const& eElt)
{
  int nbExt=eSet.size();
  Face fSet(nbExt);
  //  std::cerr << "eEltImage nbExt=" << nbExt << "\n";
  for (int iExt=0; iExt<nbExt; iExt++)
    if (eSet[iExt] == 1) {
      int jExt=eElt.at(iExt);
      //      std::cerr << "  iExt=" << iExt << " jExt=" << jExt << "\n";
      fSet[jExt]=1;
    }
  return fSet;
}


TheGroupFormat GetPermutationGroup(permlib::dom_int const& n, std::vector<permlib::Permutation> const& ListPerm)
{
  std::list<permlib::Permutation::ptr> generatorList;
  for (auto & eGen : ListPerm) {
    std::vector<permlib::dom_int> v(n);
    for (std::size_t i=0; i<n; i++)
      v[i]=eGen.at(i);
    generatorList.push_back(permlib::Permutation::ptr(new permlib::Permutation(v)));
  }
  TheGroupFormat TheGroupRet;
  TheGroupRet.n=n;
  TheGroupRet.group=construct(n, generatorList.begin(), generatorList.end());
  TheGroupRet.size=TheGroupRet.group->order<mpz_class>();
  return TheGroupRet;
}



TheGroupFormat ReadGroup(std::istream &is)
{
  if (!is.good()) {
    std::cerr << "ReadGroup operation failed because stream is not valid\n";
    throw TerminalException{1};
  }

  int nbGen;
  permlib::dom_int n;
  is >> n;
  is >> nbGen;
  std::cerr << "n=" << n << " nbGen=" << nbGen << "\n";
  std::vector<permlib::Permutation> ListGen;
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
  return GetPermutationGroup(n, ListGen);
}



void WriteGroup(std::ostream &os, TheGroupFormat const& TheGRP)
{
  std::list<permlib::Permutation::ptr> ListGen=TheGRP.group->S;
  int nbGen=ListGen.size();
  os << TheGRP.n << " " << nbGen << "\n";
  for (auto & eGen : ListGen) {
    for (std::size_t i=0; i<TheGRP.n; i++) {
      permlib::dom_int eVal=eGen->at(i);
      os << " " << eVal;
    }
    os << "\n";
  }
}

void WriteGroupMakeUp(std::ostream &os, TheGroupFormat const& TheGRP)
{
  os << "nb acting element=" << TheGRP.n << "\n";
  std::list<permlib::Permutation::ptr> ListGen=TheGRP.group->S;
  int nbGen=0;
  for (auto & eGen : ListGen) {
    for (std::size_t i=0; i<TheGRP.n; i++) {
      permlib::dom_int eVal=eGen->at(i);
      os << " " << eVal;
    }
    os << "\n";
    nbGen++;
  }
  os << "nbGen=" << nbGen << "\n";
  os << "size=" << TheGRP.group->order<mpz_class>() << "\n";
}

void WriteGroupGAP(std::ostream &os, TheGroupFormat const& TheGRP)
{
  std::list<permlib::Permutation::ptr> ListGen=TheGRP.group->S;
  os << "local eListList, ListGen, GRP;\n";
  os << "eListList:=[\n";
  bool IsFirst=true;
  for (auto & eGen : ListGen) {
    if (!IsFirst)
      os << ",\n";
    IsFirst=false;
    os << "[";
    for (std::size_t i=0; i<TheGRP.n; i++) {
      permlib::dom_int eVal=1+eGen->at(i);
      if (i>0)
	os << ",";
      os << eVal;
    }
    os << "]";
  }
  os << "];\n";
  os << "ListGen:=List(eListList, PermList);\n";
  os << "GRP:=Group(ListGen);\n";
  os << "SetSize(GRP, " << TheGRP.size << ");\n";
  os << "return GRP;\n";
}



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

std::vector<permlib::dom_int> GetBaseGroup(TheGroupFormat const& eGRP)
{
  std::vector<permlib::dom_int> eList;
  for (auto & eTrans : eGRP.group->U) {
    permlib::dom_int eElt=eTrans.element();
    eList.push_back(eElt);
  }
  return eList;
}

Face Face_EquivSimplification(Face const& eFace)
{
  int siz=eFace.size();
  int cnt=eFace.count();
  if (2*cnt > siz) {
    Face retFace(siz);
    for (int i=0; i<siz; i++) {
      int eVal;
      if (eFace[i] == 1)
	eVal=0;
      else
	eVal=1;
      retFace[i]=eVal;
    }
    return retFace;
  }
  else {
    return eFace;
  }
}





EquivTest<permlib::Permutation> TestEquivalenceGeneralNaked(TheGroupFormat const& TheGRP, Face const& eList1, Face const& eList2, int const& eMethod)
{
  //  std::cerr << "Begin of TestEquivalenceGeneral\n";
  permlib::Permutation::ptr mappingElement;
  int nb1=eList1.count();
  int nb2=eList2.count();
  if (nb1 != nb2)
    return {false, {}};
  //
  //  std::ofstream os(eFileRES);
  Face NewList1=Face_EquivSimplification(eList1);
  Face NewList2=Face_EquivSimplification(eList2);
  std::set<int> eSet1=GetSetFrom_DB(NewList1);
  std::set<int> eSet2=GetSetFrom_DB(NewList2);
  //  std::cerr << "Before call to setImage\n";
  if (eMethod == 0)
    mappingElement = permlib::setImage_classic(*TheGRP.group, eSet1.begin(), eSet1.end(), eSet2.begin(), eSet2.end());
  if (eMethod == 1)
    mappingElement = permlib::setImage_partition(*TheGRP.group, eSet1.begin(), eSet1.end(), eSet2.begin(), eSet2.end());
  //  std::cerr << "After call to setImage\n";
  if (mappingElement) {
    //    os << "IsEquivalence\n";
    //    os << "mappingElement=" << *mappingElement << "\n";
    return {true, *mappingElement};
  }
  //  os << "NOT equivalent\n";
  return {false, {}};
}

EquivTest<permlib::Permutation> TestEquivalenceGeneral(TheGroupFormat const& TheGRP, Face const& eList1, Face const& eList2)
{
  mpz_class MaxSize=10000;
  int eMethod;
  if (TheGRP.size < MaxSize) {
    eMethod=0;
  }
  else {
    eMethod=1;
  }
  return TestEquivalenceGeneralNaked(TheGRP, eList1, eList2, eMethod);
}



bool TestEquivalence(TheGroupFormat const& TheGRP, Face const& eList1, Face const& eList2)
{
  EquivTest<permlib::Permutation> eRes=TestEquivalenceGeneral(TheGRP, eList1, eList2);
  return eRes.TheReply;
}



void WriteVectorInt(std::ostream &os, std::vector<int> const& OneInc)
{
  int i, siz;
  siz=OneInc.size();
  for (i=0; i<siz; i++)
    os << " " << OneInc[i];
  os << "\n";
}


Face GROUP_Canonicalization(TheGroupFormat const& TheGRP, Face const& eList)
{
  DsetList eListI(TheGRP.n), eListO(TheGRP.n);
  int n=TheGRP.n;
  int siz=eList.count();
  int aRow=eList.find_first();
  for (int i=0; i<siz; i++) {
    eListI[aRow]=1;
    aRow=eList.find_next(aRow);
  }
  eListO=smallestSetImage(*TheGRP.group, eListI);
  Face TheRet;
  for (int i=0; i<n; i++)
    if (eListO[i] == 1)
      TheRet[i]=1;
  bool test=TestEquivalence(TheGRP, eList, TheRet);
  if (!test) {
    std::cerr << "We have major debugging to do\n";
    throw TerminalException{1};
  }
  return TheRet;
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










TheGroupFormat PERMLIB_GetStabilizer_general(TheGroupFormat const& TheGRP, Face const& eList, int const& opt)
{
  Face NewList=Face_EquivSimplification(eList);
  std::set<int> eSet=GetSetFrom_DB(NewList);
  TheGroupFormat TheStab;
  TheStab.n=TheGRP.n;
  if (opt == 0)
    TheStab.group=setStabilizer_classic(*TheGRP.group, eSet.begin(), eSet.end());
  if (opt == 1)
    TheStab.group=setStabilizer_partition(*TheGRP.group, eSet.begin(), eSet.end());
  TheStab.size=TheStab.group->order<mpz_class>();
  return TheStab;
}





TheGroupFormat GetStabilizer(TheGroupFormat const& TheGRP, Face const& eList)
{
  mpz_class MaxSize=10000;
  if (TheGRP.size < MaxSize) {
    return PERMLIB_GetStabilizer_general(TheGRP, eList,0);
  }
  else {
    return PERMLIB_GetStabilizer_general(TheGRP, eList,1);
  }
}












std::vector<int> OrbitIntersection(TheGroupFormat const& TheGRP, std::vector<int> const& gList)
{
  std::vector<int> eListReturn;
  permlib::dom_int n=TheGRP.n;
  std::vector<int> rList = gList;
  while(true) {
    std::size_t eSumPrev=0;
    for (std::size_t i=0; i<n; i++)
      eSumPrev += rList[i];
    for (auto & eGen : TheGRP.group->S) {
      for (std::size_t i=0; i<n; i++) {
	std::size_t j=eGen->at(i);
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


std::vector<int> OrbitUnion(TheGroupFormat const& TheGRP, std::vector<int> const& gList)
{
  int n=TheGRP.n;
  std::vector<int> gListB(n);
  for (int i=0; i<n; i++)
    gListB[i]=1-gList[i];
  std::vector<int> rListB=OrbitIntersection(TheGRP, gListB);
  for (int i=0; i<n; i++)
    rListB[i]=1-rListB[i];
  return rListB;
}



TheGroupFormat ReducedGroupAction(TheGroupFormat const& TheGRP, Face const& eList)
{
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
  std::vector<permlib::Permutation> ListGen;
  for (auto & eGen : TheGRP.group->S) {
    std::vector<permlib::dom_int> v(nb);
    for (int i=0; i<nb; i++) {
      int eVal1=ListPosition[i];
      int eVal2=eGen->at(eVal1);
      int eVal3=ListPositionRev[eVal2];
      v[i]=eVal3;
    }
    ListGen.push_back(permlib::Permutation(v));
  }
  return GetPermutationGroup(nb, ListGen);
}

TheGroupFormat ConjugateGroup(TheGroupFormat const& TheGRP, permlib::Permutation const& ePerm)
{
  int n=TheGRP.n;
  permlib::Permutation ePermR=~ePerm;
  std::vector<permlib::Permutation> ListGen;
  for (auto & eGen : TheGRP.group->S) {
    std::vector<permlib::dom_int> v(n);
    for (int i=0; i<n; i++) {
      int eVal1=ePermR.at(i);
      int eVal2=eGen->at(eVal1);
      int eVal3=ePerm.at(eVal2);
      v[i]=eVal3;
    }
    ListGen.push_back(permlib::Permutation(v));
  }
  return GetPermutationGroup(n, ListGen);
}








std::list<int> ComputeFullOrbitPoint(TheGroupFormat const& TheGRP, int const& ePoint)
{
  IntegerSubsetStorage *Vorb;
  IntegerSubsetStorage *Vactive;
  std::list<int> eList;
  std::list<permlib::Permutation::ptr> ListGen;
  int TheFirst;
  permlib::dom_int n=TheGRP.n;
  ListGen=TheGRP.group->S;
  Vorb=(IntegerSubsetStorage*)malloc(sizeof(IntegerSubsetStorage));
  Vactive=(IntegerSubsetStorage*)malloc(sizeof(IntegerSubsetStorage));
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
  free(Vorb);
  free(Vactive);
  return eList;
}

std::vector<Face> DecomposeOrbitPoint(TheGroupFormat const& TheGRP, Face const& eList)
{
  IntegerSubsetStorage *Vlist;
  std::vector<Face> ListOrb;
  int len, i, TheFirst;
  Vlist=(IntegerSubsetStorage*)malloc(sizeof(IntegerSubsetStorage));
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
  free(Vlist);
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
  std::set<permlib::Permutation> ListPermWork;
  for (std::shared_ptr<permlib::Permutation> & p : eTrans.GetMtransversal() ) {
    if (p) {
      permlib::Permutation ePerm=*p;
      ListPermWork.insert(ePerm);
    }
  }
  std::set<permlib::dom_int> SetOrbit;
  std::vector<PairEltPerm> ListPair;
  std::vector<bool> StatusDone;
  std::function<void(permlib::dom_int,permlib::Permutation)> fInsert=[&](permlib::dom_int const& eVal, permlib::Permutation const& ePerm) -> void {
    std::set<permlib::dom_int>::iterator iter=SetOrbit.find(eVal);
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
  //  std::cerr << "Beginning of WritePermutationGAP\n";
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
  //  std::cerr << "End of WritePermutationGAP\n";
}








#endif
