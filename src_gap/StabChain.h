#ifndef DEFINE_STAB_CHAIN
#define DEFINE_STAB_CHAIN

#define DEBUG_GROUP

/*
  Version of the code used. It must be from Summer 2016.
 */


/*
  Questions on the GAP code:
  ---Why QuickInverseRepresentative is not used more if it is the same but maybe quicker
     than InverseRepresentative

  Not needed for the time being:
  ---MembershipTestKnownBase needed for PCGS
  ---MinimalElementCosetStabChain for coset representatives.
  ---ListStabChain since we are using a different design in our code.
  ---OrbitStabChain which seems used by no code at all.
  ---StabChainForcePoint seems not to be used
  ---knownBase seems strange. We should be able to work completely without it,
     since I never set up the knownBase in the first place.
  ---

  For what it seems:
  ---The transversal entry are permutation (when assigned) and all included in the
     S.labels
  ---The S.labels are all identical from the top to bottom.
  ---Significantly, this seems also the case of S.identity.
  ---The attribute "relativeOrders" is related to pcgs and not needed here.
*/

#include<vector>


#include "PermGroup.h"



namespace gap {


  
template<typename Telt>
int GetLabelIndex(std::vector<Telt> & labels, Telt const& u)
{
  int nbLabel=labels.size();
  for (int iLabel=0; iLabel<nbLabel; iLabel++)
    if (labels[iLabel] == u)
      return iLabel;
  labels.push_back(u);
  return nbLabel;
}

template<typename Telt>
int GetLabelIndex_const(std::vector<Telt> const& labels, Telt const& u)
{
  int nbLabel=labels.size();
  for (int iLabel=0; iLabel<nbLabel; iLabel++)
    if (labels[iLabel] == u)
      return iLabel;
  return -1;
}


template<typename Telt>
struct StabLevel {
  //  std::vector<Telt> generators;  // We should use genlabels
  std::vector<int> transversal;
  std::vector<int> orbit;
  std::vector<int> genlabels; // Not used in Random algorithm
  Face cycles;
  // entry below are specific to the random algorithms:
  std::vector<Telt> treegen;
  std::vector<Telt> treegeninv;
  std::vector<Telt> aux;
  int treedepth;
  int diam;
};
// other possible entries:
// transimages, genimages, labelimages, idimage



 
// The labels are put on top since they are all identical.
// The stabilizers can be added in any way:
// ---At the top of the chain
// ---Removed in the middle
// Therefore we need a data set that allows us to do that and
// native std::vector are not adequate.
template<typename Telt>
struct StabChain {
  int n;
  Telt identity;
  bool UseCycle;
  std::vector<Telt> labels;
  std::vector<StabLevel<Telt>> stabilizer;
};

template<typename Telt>
std::ostream& operator<<(std::ostream& os, StabChain<Telt> const& Stot)
{
  os << "n=" << Stot.n << " identity=" << GapStyleString(Stot.identity) << "\n";
  os << "labels=[";
  bool IsFirst=true;
  for (auto & eLabel : Stot.labels) {
    if (!IsFirst)
      os << " , ";
    IsFirst=false;
    os << GapStyleString(eLabel);
  }
  os << "]\n";
  int nbLev=Stot.stabilizer.size();
  for (int iLev=0; iLev<nbLev; iLev++) {
    os << "iLev=" << iLev << "\n";
    os << "  transversal =";
    for (auto & eVal : Stot.stabilizer[iLev].transversal)
      os << " " << eVal;
    os << "\n";
    os << "  orbit=";
    for (auto & eVal : Stot.stabilizer[iLev].orbit)
      os << " " << eVal;
    os << "\n";
    os << "  genlabels=";
    for (auto & eVal : Stot.stabilizer[iLev].genlabels)
      os << " " << eVal;
    os << "\n";
  }
  return os;
}



 
template<typename Telt>
StabChain<Telt> RestrictedStabChain(StabChain<Telt> const& Stot, int const& eLev)
{
  int nbLevel=Stot.stabilizer.size();
  std::vector<StabLevel<Telt>> stabilizerRed;
  for (int uLev=eLev; uLev<nbLevel; uLev++)
    stabilizerRed.push_back(Stot.stabilizer[uLev]);
  return {Stot.n, Stot.identity, Stot.UseCycle, Stot.labels, stabilizerRed};
}


template<typename Telt>
StabLevel<Telt> EmptyStabLevel()
{
  std::vector<int> transversal;
  std::vector<int> orbit;
  std::vector<int> genlabels;
  Face cycles;
  std::vector<Telt> treegen;
  std::vector<Telt> treegeninv;
  std::vector<Telt> aux;
  int treedepth = 0;
  int diam = 0;
  return {transversal, orbit, genlabels, cycles, treegen, treegeninv, aux, treedepth, diam};
}


 
template<typename Telt>
StabChain<Telt> EmptyStabChain(int const& n)
{
  Telt id(n);
  return {n, id, false, {id}, {EmptyStabLevel<Telt>()} };
}


template<typename Telt>
StabChain<Telt> EmptyStabChainPlusNode(int const& n, int const& bas)
{
  Telt id(n);
  StabLevel<Telt> eLevel = EmptyStabLevel<Telt>();
  eLevel.orbit.push_back(bas);
  return {n, id, false, {id}, {eLevel} };
}


 
 
template<typename Telt>
int BasePoint(StabChain<Telt> const& Stot, int const& eLev)
{
  int siz=Stot.stabilizer.size();
  if (eLev >= siz)
    return -1;
  return Stot.stabilizer[eLev].orbit[0];
}

/*
template<typename Telt>
void RemoveOneStabilizer(StabChain<Telt> & S, int const& idx)
{
  int nbInd=S.ListIndices.size();
  std::vector<int> NewListIndices(nbInd-1);
  int pos=0;
  for (int iInd=0; iInd<nbInd; iInd++) {
    if (iInd != idx) {
      NewListIndices[pos]=S.ListIndices[iInd];
      pos++;
    }
  }
  S.ListIndices=NewListIndices;
}
*/


/* 
template<typename Telt>
void InsertStabilizer(StabChain<Telt> & S, int const& ThePos)
{
  int nbInd=S.ListIndices.size();
  int nbStab=S.stabilizer.size();
  if (nbStab <= nbInd) {
    S.stabilizer.push_back({});
    nbStab++;
  }
  std::vector<int> AttainedIndex(nbStab,0);
  for (int iInd=0; iInd<nbInd; iInd++) {
    int pos=S.ListIndices[iInd];
    AttainedIndex[pos]=1;
  }
  int FirstFree=-1;
  for (int iStab=0; iStab<nbStab; iStab++)
    if (FirstFree == -1 && AttainedIndex[iStab] == 0)
      FirstFree=iStab;
  std::vector<int> NewListIndices(nbInd+1);
  int pos=0;
  for (int iInd=0; iInd<nbInd; iInd++) {
    if (iInd == ThePos) {
      NewListIndices[pos]=FirstFree;
      pos++;
    }
    NewListIndices[pos]=ListIndices[iInd];
    pos++;
  }
  S.ListIndices=NewListIndices;
}
*/

 
 

   // almost certainly wrong code below.
template<typename Telt>
void RemoveStabChain(StabChain<Telt> & Stot)
{
  Telt TheId=Stot.identity;
  Stot.stabilizer.clear();
  Stot.stabilizer[0].genlabels.clear();
}


template<typename Telt>
bool IsInBasicOrbit(StabChain<Telt> const& eStab, int const& lev, int const& pnt)
{
  int eVal=eStab.stabilizer[lev].transversal[pnt];
  if (eVal == -1)
    return false;
  return true;
}

// This correspond to the code of
// pnt^g in GAP code.
template<typename Telt>
int PowAct(int const& pnt, Telt const& g)
{
  return g.at(pnt);
}


template<typename Telt>
Telt InverseRepresentative(StabChain<Telt> const& Stot, int const& eLev, int const& pnt)
{
  int bpt=Stot.stabilizer[eLev].orbit[0];
  Telt rep=Stot.identity;
  int pntw=pnt;
  //  std::cerr << "bpt=" << bpt << " pnt=" << pnt << "\n";
  while(pntw != bpt) {
    int idx=Stot.stabilizer[eLev].transversal[pntw];
    //    std::cerr << "idx=" << idx << "\n";
    Telt te=Stot.labels[idx];
    //    std::cerr << "te = " << GapStyleString(te) << "\n";
    //    std::cerr << "Before pntw=" << pntw << "\n";
    pntw=PowAct(pntw, te);
    //    std::cerr << " After pntw=" << pntw << "\n";
    rep = rep * te;
  }
  return rep;
}


template<typename Telt>
std::vector<Telt> InverseRepresentativeWord(StabChain<Telt> const& Stot, int const& eLev, int const& pnt)
{
  int bpt=Stot.stabilizer[eLev].orbit[0];
  std::vector<Telt> word;
  int pntw=pnt;
  while(pntw != bpt) {
    int idx=Stot.stabilizer[eLev].transversal[pntw];
    Telt te=Stot.labels[idx];
    pntw=PowAct(pntw, te);
    word.push_back(te);
  }
  return word;
}


template<typename Telt>
Telt SiftedPermutation(StabChain<Telt> const& Stot, int const& eLev, Telt const& g)
{
  int len=Stot.stabilizer.size()-1;
  Telt gW=g;
  for (int iLev=eLev; iLev<len; iLev++) {
    if (gW.isIdentity())
      return gW;
    int bpt=Stot.stabilizer[iLev].orbit[0];
    int img=PowAct(bpt, gW);
    if (Stot.stabilizer[iLev].transversal[img] == -1)
      return gW;
    while(true) {
      if (img == bpt)
	break;
      int idx=Stot.stabilizer[iLev].transversal[img];
      gW = gW * Stot.labels[idx];
      img = PowAct(bpt, gW);
    }
  }
  return gW;
}


template<typename Telt>
std::vector<int> BaseStabChain(StabChain<Telt> const& Stot)
{
  int len=Stot.stabilizer.size();
  std::vector<int> base(len);
  for (int iLev=0; iLev<len; iLev++) {
    int eVal=Stot.stabilizer[iLev].orbit[0];
    base[iLev]=eVal;
  }
  return base;
}



template<typename Telt, typename Tint>
Tint SizeStabChain(StabChain<Telt> const& Stot)
{
  Tint size=1;
  int len=Stot.stabilizer.size();
  for (int iLev=0; iLev<len; iLev++) {
    int siz=Stot.stabilizer[iLev].orbit.size();
    Tint siz_i = siz;
    size *= siz_i;
  }
  return size;
}
 
template<typename Telt>
std::vector<Telt> StrongGeneratorsStabChain(StabChain<Telt> const& Stot, int const& TheLev)
{
  std::set<Telt> sgs_set;
  int len=Stot.stabilizer.size();
  for (int iLev=TheLev; iLev<len; iLev++) {
    int siz=Stot.stabilizer[iLev].genlabels.size();
    if (siz == 0)
      break;
    for (auto & pos : Stot.stabilizer[iLev].genlabels)
      sgs_set.insert(Stot.labels[pos]);
  }
  std::vector<Telt> sgs(sgs_set.begin(), sgs_set.end());
  return sgs;
}


template<typename Telt>
std::vector<int> IndicesStabChain(StabChain<Telt> const& Stot)
{
  int len=Stot.stabilizer.size();
  std::vector<int> ind(len);
  for (int iLev=0; iLev<len; iLev++) {
    int siz=Stot.stabilizer[iLev].orbit.size();
    ind[iLev]=siz;
  }
  return ind;
}



template<typename Telt>
Telt LargestElementStabChain(StabChain<Telt> const& Stot)
{
  Telt rep=Stot.stabilizer[0].identity;
  int len=Stot.stabilizer.size();
  for (int iLev=0; iLev<len; iLev++) {
    if (Stot.stabilizer[iLev].genlabels.size() == 0)
      break;
    int pnt=Stot.stabilizer[iLev].orbit[0];
    int min=0;
    int val=0;
    for (auto & i : Stot.stabilizer[iLev].orbit) {
      int img=PowAct(i, rep);
      if (img > val) {
	min=i;
	val=img;
      }
    }
    while(true) {
      if (pnt == min)
	break;
      int idx=Stot.stabilizer[iLev].transversal[min];
      Telt gen=Stot.labels[idx];
      rep=LeftQuotient(gen,rep);
      min=PowAct(min, gen);
    }
  }
  return rep;
}

 

template<typename Telt>
std::vector<Telt> ElementsStabChain(StabChain<Telt> const& Stot)
{
  std::vector<Telt> elms;
  auto LevelIncrease=[&](StabLevel<Telt> const& eLev) -> void {
    std::vector<Telt> NewElms;
    for (auto & pnt : eLev.orbit) {
      Telt rep=eLev.identity;
      while (PowAct(eLev.orbit[0], rep) != pnt) {
	int jpt=SlashAct(pnt, rep);
	int idx=Stot.stabilizer[eLev].transversal[jpt];
	rep=LeftQuotient(Stot.labels[idx], rep);
      }
      for (auto & eStb : elms)
	NewElms.push_back(eStb * rep);
    }
    elms=NewElms;
  };
  elms={Stot.stabilizer[0].identity};
  int len=Stot.stabilizer.size();
  for (int iLev=0; iLev<len; iLev++) {
    int jLev=len-1-iLev;
    LevelIncrease(Stot.stabilizer[jLev]);
  }
  return elms;
}


// is base is empty then this just replaces the IsBound(options.base)
template<typename Tint>
struct StabChainOptions {
  std::vector<int> base;
  std::vector<int> knownBase;
  int random;
  bool reduced;
  Tint size;
  Tint limit;
};


template<typename Tint>
StabChainOptions<Tint> GetStandardOptions()
{
  std::vector<int> base;
  std::vector<int> knownBase;
  int random = 1000;
  bool reduced=true;
  Tint size=0;
  Tint limit=0;
  return {base, knownBase, random, reduced, size, limit};
}


template<typename Telt>
bool IsTrivial(std::vector<Telt> const& G)
{
  for (auto & eElt : G)
    if (!eElt.isIdentity())
      return false;
  return true;
}

template<tpename Telt>
std::vector<int> MovedPoints(StabChain<Telt> const& S)
{
  std::set<int> LIdx;
  for (auto & eChain : S.stabilizer)
    for (auto & eIdx : eChain.genlabels)
      LIdx.insert(eIdx);
  int n=S.n;
  auto IsMoved=[&](int const& ePt) -> bool {
    for (auto & eIdx : LIdx) {
      if (S.labels[eIdx].at(ePt) != ePt)
	return false;
    }
    return true;
  };
  std::vector<int> LMoved;
  for (int i=0; i<n; i++)
    if (IsMoved(i))
      LMoved.push_back(i);
  return LMoved;
}



template<typename Telt>
int LargestMovedPoint(std::vector<Telt> const& LGen)
{
  int n=LGen[0].size();
  std::vector<int> Status(n, 1);
  for (auto & eGen : LGen) {
    for (int u=0; u<n; u++) {
      int v=eGen.at(u);
      if (u != v)
	Status[u]=0;
    }
  }
  int eMov=0;
  for (int u=0; u<n; u++) {
    if (Status[u] == 0)
      eMov=u;
  }
  eMov++;
  return eMov;
}

// Most likely very buggy
template<typename Telt>
void InsertTrivialStabilizer(StabChain<Telt> & Stot, int const& eLev, int const& pnt)
{
  StabLevel<Telt> eStab;
  Stot.stabilizer.emplace(Stot.stabilizer.begin() + eLev + 1, eStab);
  Stot.stabilizer[eLev+1] = Stot.stabilizer[eLev];
  InitializeSchreierTree(Stot, eLev, pnt);
}



template<typename Telt>
StabChain<Telt> StabChainBaseStrongGenerators(std::vector<int> const& base, std::vector<Telt> const& sgs)
{
  int n=sgs[0].size();
  StabChain<Telt> Stot = EmptyStabChain<Telt>(n);
  int nbGen=sgs.size();
  Face status(nbGen);
  for (int i=0; i<nbGen; i++)
    status[i] = 1;
  int basSiz=base.size();
  for (int iBas=0; iBas<basSiz; iBas++) {
    int pnt=base[iBas];
    std::vector<Telt> sgsFilt;
    for (int i=0; i<nbGen; i++)
      if (status[i] == 1)
	sgsFilt.push_back(sgs[i]);
    InsertTrivialStabilizer(Stot, iBas, pnt);
    AddGeneratorsExtendSchreierTree(Stot, iBas, sgsFilt);
    for (int i=0; i<nbGen; i++)
      if (status[i] == 1 && PowAct(pnt, sgs[i]) != pnt)
	status[i]=0;
  }
  return Stot;
}


template<typename Telt>
void AddGeneratorsExtendSchreierTree(StabChain<Telt> & Stot, int const& eLev, std::vector<Telt> const& newgens)
{
  int nbLabel=Stot.labels.size();
  std::vector<int> ListAtt(nbLabel);
  for (int i=0; i<nbLabel; i++)
    ListAtt[i]=i;
  Face old=BlistList(ListAtt, Stot.stabilizer[eLev].genlabels);
  old[0]=true;
  std::cerr << "old =";
  for (int i=0; i<int(old.size()); i++)
    std::cerr << " " << old[i];
  std::cerr << "\n";
  std::cerr << "Before genlabels =";
  for (auto & eVal : Stot.stabilizer[eLev].genlabels)
    std::cerr << " " << eVal;
  std::cerr << "\n";
  Face ald=old;
  for (auto & gen : newgens) {
    int pos = PositionVect(Stot.labels, gen);
    if (pos == -1) {
      Stot.labels.push_back(gen);
      old.push_back(false);
      ald.push_back(true);
      int posG=Stot.labels.size() - 1;
      Stot.stabilizer[eLev].genlabels.push_back(posG);
    }
    else {
      if (!ald[pos])
	Stot.stabilizer[eLev].genlabels.push_back(pos);
    }
  }
  std::cerr << "Before test old =";
  for (int i=0; i<int(old.size()); i++)
    std::cerr << " " << old[i];
  std::cerr << "\n";
  std::cerr << "transversal =";
  for (int i=0; i<int(Stot.stabilizer[eLev].transversal.size()); i++)
    std::cerr << " " << Stot.stabilizer[eLev].transversal[i];
  std::cerr << "\n";
  std::cerr << "After genlabels =";
  for (auto & eVal : Stot.stabilizer[eLev].genlabels)
    std::cerr << " " << eVal;
  std::cerr << "\n";
  int len = Stot.stabilizer[eLev].orbit.size();
  int i=0;
  if (Stot.UseCycle) {
    std::cerr << "Before cycles =";
    for (int u=0; u<int(Stot.stabilizer[eLev].cycles.size()); u++)
      std::cerr << " " << Stot.stabilizer[eLev].cycles[u];
    std::cerr << "\n";
    while (i < int(Stot.stabilizer[eLev].orbit.size())) {
      std::cerr << "i=" << i << "\n";
      for (int& j : Stot.stabilizer[eLev].genlabels) {
	std::cerr << "  j=" << j << "\n";
	if (i > len-1 || old[j] == 0) {
	  int img=SlashAct(Stot.stabilizer[eLev].orbit[i], Stot.labels[j]);
	  std::cerr << "    After the test img=" << img << "\n";
	  if (Stot.stabilizer[eLev].transversal[img] != -1) {
	    Stot.stabilizer[eLev].cycles[i]=true;
	  }
	  else {
	    Stot.stabilizer[eLev].transversal[img]=j;
	    Stot.stabilizer[eLev].orbit.push_back(img);
	    Stot.stabilizer[eLev].cycles.push_back(false);
	  }
	}
      }
      i++;
    }
    std::cerr << "After cycles =";
    for (int u=0; u<int(Stot.stabilizer[eLev].cycles.size()); u++)
      std::cerr << " " << Stot.stabilizer[eLev].cycles[u];
    std::cerr << "\n";
  }
  else {
    while (i < int(Stot.stabilizer[eLev].orbit.size())) {
      for (int& j : Stot.stabilizer[eLev].genlabels) {
	if (i > len || old[j] == 0) {
	  int img=SlashAct(Stot.stabilizer[eLev].orbit[i], Stot.labels[j]);
	  if (Stot.stabilizer[eLev].transversal[img] == -1) {
	    Stot.stabilizer[eLev].transversal[img]=j;
	    Stot.stabilizer[eLev].orbit.push_back(img);
	  }
	}
      }
      i++;
    }
  }
}

template<typename Telt>
void ChooseNextBasePoint(StabChain<Telt> & Stot, int const& eLev, std::vector<int> const& base, std::vector<Telt> const& newgens)
{
  std::cerr << "base =";
  for (auto & eVal : base)
    std::cerr << " " << eVal;
  std::cerr << "\n";
  auto IsFullyStable=[&](int const& eBas) -> bool {
    for (auto & eGen : newgens) {
      if (PowAct(eBas, eGen) != eBas)
	return false;
    }
    return true;
  };
  int i = 0;
  int len=base.size();
  while(true) {
    if (i == len)
      break;
    int eBas=base[i];
    if (!IsFullyStable(eBas))
      break;
    i++;
  }
  int pnt;
  if (i < len)
    pnt=base[i];
  else
    pnt=SmallestMovedPoint(newgens);
  int bpt, pos;
  std::cerr << "eLev=" << eLev << " |Stot.stabilizer|=" << Stot.stabilizer.size() << "\n";
  if (Stot.stabilizer[eLev].orbit.size() > 0) {
    bpt = Stot.stabilizer[eLev].orbit[0];
    pos = PositionVect(base, bpt);
  }
  else {
    bpt = Stot.n + 444; // value in GAP is infinity
    pos = -1;
  }
  std::cerr << "BPT/POS bpt=" << bpt << " pos=" << pos << "\n";
  if ((pos != -1 && i < pos) || (pos == -1 && i<int(base.size())) || (pos == -1 && pnt < bpt)) {
    std::cerr << "InsertTrivialStabilizer pnt=" << pnt << " bpt=" << bpt << " pos=" << pos << "\n";
    InsertTrivialStabilizer(Stot, eLev, pnt);
    if (Stot.UseCycle) {
      Face eFace(1);
      eFace[0] = 0;
      Stot.stabilizer[eLev].cycles = eFace;
      std::cerr << "   Initializing cycles\n";
    }
  }
}



template<typename Telt, typename Tint>
void StabChainStrong(StabChain<Telt> & Stot, int const& eLev, std::vector<Telt> const& newgens, StabChainOptions<Tint> const& options)
{
  std::cerr << "   StabChainStrong call at eLev=" << eLev << " |newgens|=" << newgens.size() << "\n";
  ChooseNextBasePoint(Stot, eLev, options.base, newgens);
  std::cerr << "|Stot.stabilizer|=" << Stot.stabilizer.size() << " eLev=" << eLev << "\n";
  std::cerr << "|Stot.stabilizer[eLev].orbit|=" << Stot.stabilizer[eLev].orbit.size() << "\n";
  
  int pnt = Stot.stabilizer[eLev].orbit[0];
  int len = Stot.stabilizer[eLev].orbit.size();
  int old = Stot.stabilizer[eLev].genlabels.size();
  AddGeneratorsExtendSchreierTree(Stot, eLev, newgens);
  
  //# If a new generator fixes the base point, put it into the stabilizer.
  for (auto & eGen : newgens)
    if (eGen.isIdentity() == false && PowAct(pnt, eGen) == pnt) {
      std::cerr << "   1: Calling StabChainStrong with eGen=" << GapStyleString(eGen) << "\n";
      StabChainStrong(Stot, eLev+1, {eGen}, options);
    }
    
  // # Compute the Schreier generators (seems to work better backwards).
  std::vector<int> pnts = ClosedInterval(0, Stot.stabilizer[eLev].orbit.size());
  if (Stot.UseCycle)
    pnts=ListBlist(pnts, Stot.stabilizer[eLev].cycles);
  std::cerr << "   pnts =";
  for (auto & eVal : pnts)
    std::cerr << " " << eVal;
  std::cerr << " Usecycle=" << Stot.UseCycle;
  if (Stot.UseCycle) {
    std::cerr << " cycles=";
    for (int i=0; i<int(Stot.stabilizer[eLev].orbit.size()); i++)
      std::cerr << " " << Stot.stabilizer[eLev].cycles[i];
  }
  std::cerr << "\n";
  int gen1=0;
  for (int& i : Reversed(pnts)) {
    int p=Stot.stabilizer[eLev].orbit[i];
    Telt rep=InverseRepresentative(Stot, eLev, p );
    if (i < len)
      gen1=old;
    for (int & j : ClosedInterval(gen1, Stot.stabilizer[eLev].genlabels.size())) {
      Telt g = Stot.labels[ Stot.stabilizer[eLev].genlabels[j] ];
      if (Stot.stabilizer[eLev].transversal[ SlashAct(p, g) ] != Stot.stabilizer[eLev].genlabels[j]) {
        Telt sch = SiftedPermutation(Stot, eLev, Inverse(g*rep));
	if (!sch.isIdentity()) {
	  std::cerr << "   2: Calling StabChainStrong with sch=" << GapStyleString(sch) << "\n";
	  StabChainStrong(Stot, eLev+1, {sch}, options );
	}
      }
    }
  }
}


template<typename Telt>
bool StabChainForcePoint(StabChain<Telt> & Stot, int const& eLev, int const& pnt)
{
  if (Stot.stabilizer[eLev].transversal[pnt] == -1) {
    if (IsFixedStabilizer(Stot, eLev, pnt )) {
      InsertTrivialStabilizer(Stot, eLev, pnt);
    }
    else {
      if (!StabChainForcePoint(Stot, eLev+1, pnt) || !StabChainSwap(Stot, eLev))
	return false;
    }
  }
  return true;
}


template<typename Telt>
std::vector<Telt> GetListGenerators(StabChain<Telt> const& Stot, int const& eLev)
{
  std::vector<Telt> LGens;
  for (auto & posGen : Stot.stabilizer[eLev].genlabels)
    LGens.push_back(Stot.labels[posGen]);
  return LGens;
}

template<typename Telt>
bool StabChainSwap(StabChain<Telt> & Stot, int const& eLev)
{
  std::cerr << "   Running StabChainSwap\n";
  int n=Stot.n;
  int a = Stot.stabilizer[eLev].orbit[0];
  int b = Stot.stabilizer[eLev+1].orbit[0];
  //
  std::vector<Telt> LGens = GetListGenerators(Stot, eLev);
  //
  StabChain<Telt> Ttot = EmptyStabChainPlusNode<Telt>(n, b);
  AddGeneratorsExtendSchreierTree(Ttot, 0, LGens);
  //
  StabChain<Telt> Tstab = EmptyStabChainPlusNode<Telt>(n, a);
  int nbLev=Stot.stabilizer.size();
  if (eLev+2 < nbLev) {
    std::vector<Telt> LGensB = GetListGenerators(Stot, eLev+2);
    AddGeneratorsExtendSchreierTree(Tstab, 0, LGensB);
  }
  //
  int ind = 0;
  int len = Stot.stabilizer[eLev].orbit.size() * Stot.stabilizer[eLev+1].orbit.size() / Ttot.stabilizer[0].orbit.size();
  while (int(Tstab.stabilizer[0].orbit.size()) < len) {
    int pnt;
    while(true) {
      ind++;
      if (ind > int(Stot.stabilizer[eLev].orbit.size()))
	return false;
      pnt = Stot.stabilizer[eLev].orbit[ind];
      if (Tstab.stabilizer[0].transversal[pnt] == -1)
	break;
    }
    int img = b;
    int i = pnt;
    while (i != a) {
      int posGen=Stot.stabilizer[eLev].transversal[i];
      img = PowAct(img, Stot.labels[posGen]);
      i = PowAct(i, Stot.labels[posGen]);
    }
    if (Stot.stabilizer[eLev+1].transversal[img] != -1) {
      Telt gen = Stot.identity;
      while (PowAct(pnt, gen) != a) {
	int posGen=Stot.stabilizer[eLev].transversal[PowAct(pnt, gen)];
	gen = gen * Stot.labels[posGen];
      }
      while (PowAct(b, gen) != b) {
	int posGen=Stot.stabilizer[eLev+1].transversal[PowAct(pnt, gen)];
	gen = gen * Stot.labels[posGen];
      }
      AddGeneratorsExtendSchreierTree(Tstab, 0, {gen});
    }
  }
  auto MappingIndex=[&](StabChain<Telt> const& Wtot, int const& idx) -> int {
    if (idx == -1)
      return -1;
    Telt eElt=Wtot.labels[idx];
    return GetLabelIndex(Stot.labels, eElt);
  };
  auto MapAtLevel=[&](StabChain<Telt> const& Wtot, int const& TargetLev) -> void {
    Stot.stabilizer[TargetLev].genlabels.clear();
    for (int const& posGen : Wtot.stabilizer[0].genlabels) {
      int posGenMap=MappingIndex(Wtot, posGen);
      Stot.stabilizer[TargetLev].genlabels.push_back(posGenMap);
    }
    Stot.stabilizer[TargetLev].orbit = Wtot.stabilizer[0].orbit;
    for (int u=0; u<n; u++) {
      int idx=Wtot.stabilizer[0].transversal[u];
      int idxMap=MappingIndex(Wtot, idx);
      Stot.stabilizer[TargetLev].transversal[u] = idxMap;
    }
  };
  MapAtLevel(Ttot, eLev);
  if (Tstab.stabilizer[0].orbit.size() == 1)
    Stot.stabilizer.erase(Stot.stabilizer.begin() + eLev + 1);
  else
    MapAtLevel(Tstab, eLev+1);
  return true;
}





// maybe use std::map<T, T> instead
template<typename T>
T LabsLims(T const& lab, std::function<T(T const&)> const& hom, std::vector<T> & labs, std::vector<T> & lims)
{
  int pos=PositionVect(labs, lab);
  if (pos == -1) {
    int pos=labs.size();
    labs.push_back(lab);
    T img=hom(lab);
    lims.push_back(img);
  }
  return lims[pos];
}

// We significantly change the functionality 
template<typename Telt>
void ConjugateStabChain(StabChain<Telt> & Stot, int const& TheLev, Telt const& cnj)
{
  int n=Stot.n;
  int nbLev=Stot.stabilizer.size();
  for (int uLev=TheLev; uLev<nbLev; uLev++) {
    
  }
}

static const int int

// value of reduced
//  reduced = -1 corresponds to reduced = -1 in GAP code
//  reduced = 0  corresponds to reduced = false in GAP code
//  reduced = 1  corresponds to reduced = true in GAP code
template<typename Telt>
bool ChangeStabChain(StabChain<Telt> & Stot, int const& TheLev, std::vector<int> const& base, int const& reduced)
{
  Telt cnj = Stot.identity;
  std::vector<int> newBase;
  int i=0;
  int eLev=TheLev;
  int basSiz=base.size();
  while (eLev < int(Stot.stabilizer.size())-1 || i < basSiz) {
    int old=BasePoint(Stot, eLev);
    if (Stot.stabilizer[eLev].genlabels.size() == 0 && (reduced == int_true || i >= basSiz)) {
      RemoveStabChain(Stot);
      i = basSiz;
    }
    else if (i < basSiz) {
      int newpnt = SlashAct(base[i], cnj);
      i++;
      if (reduced == int_reducedm1) {
	newBase.push_back(newpnt);
	if (newpnt != old) {
	  if (IsFixedStabilizer(Stot, eLev, newpnt)) {
	    InsertTrivialStabilizer(Stot, eLev, newpnt);
	  }
#ifdef DEBUG_GROUP
	  else {
	    std::cerr << "<base> must be an extension of base of <G>\n";
	    throw TerminalException{1};
	  }
#endif
	}
	eLev++;
      }
      else if (reduced == int_false || !IsFixedStabilizer(Stot, eLev, newpnt )) {
	if (eLev < int(Stot.stabilizer.size())-1) {
	  if (!StabChainForcePoint(Stot, eLev, newpnt))
	    return false;
	  cnj = LeftQuotient(InverseRepresentative(Stot, eLev, newpnt), cnj);
	}
	else {
	  InsertTrivialStabilizer(Stot, eLev, newpnt);
	}
	newBase.push_back(Stot.stabilizer[eLev].orbit[0]);
	eLev++;
      }
    }
    else if (PositionVect(newBase, old) != -1 || (reduced == int_true && Stot.stabilizer[eLev].orbit.size() == 1)) {
      int nbStab=Stot.stabilizer.size();
      for (int u=eLev; u<nbStab-1; u++)
	Stot.stabilizer[u] = Stot.stabilizer[u+1];
      Stot.stabilizer.pop_back();
    }
    else {
      newBase.push_back(old);
      eLev++;
    }
  }
  if (!cnj.isIdentity())
    ConjugateStabChain(Stot, TheLev, cnj);
  return true;
}

template<typename Telt>
bool ExtendStabChain(StabChain<Telt> & Stot, int const & TheLev, std::vector<int> const& base)
{
  return ChangeStabChain(Stot, TheLev, base, int_reducedm1);
}


template<typename Telt>
bool ReduceStabChain(StabChain<Telt> & Stot, int const& TheLev)
{
  return ChangeStabChain(Stot, TheLev, {}, int_true);
}


template<typename Telt>
void InitializeSchreierTree(StabChain<Telt> & Stot, int const& eLev, int const& pnt)
{
  int n=Stot.n;
  //
  std::vector<int> transversal(n, -1);
  transversal[pnt] = 0;
  //
  Stot.stabilizer[eLev].orbit = {pnt};
  Stot.stabilizer[eLev].transversal = transversal;
}
    

template<typename Telt>
bool TestEqualityAtLevel(StabChain<Telt> const& L, StabChain<Telt> const& R, int const& lev)
{
  int nbLevL=L.stabilizer.size();
  int nbLevR=R.stabilizer.size();
  if (nbLevL != nbLevR)
    return false;
  for (int eLev=lev; eLev<nbLevL; eLev++) {
    if (L.stabilizer[eLev].orbit != R.stabilizer[eLev].orbit)
      return false;
    if (L.stabilizer[eLev].transversal.size() != R.stabilizer[eLev].transversal.size())
      return false;
    ine lenL=L.stabilizer[eLev].transversal.size();
    if (int u=0; u<lenL; u++) {
      if (L.stabilizer[eLev].transversal[u] == -1 && R.stabilizer[eLev].transversal[u] != -1)
	return false;
      if (L.stabilizer[eLev].transversal[u] != -1 && R.stabilizer[eLev].transversal[u] == -1)
	return false;
      if (L.stabilizer[eLev].transversal[u] != -1) {
	int idxL=L.stabilizer[eLev].transversal[u];
	int idxR=R.stabilizer[eLev].transversal[u];
	Telt permL=L.labels[idxL];
	Telt permR=L.labels[idxL];
	if (permL != permR)
	  return false;
      }
    }
  }
  return true;
}


template<typename Telt>
void SetStabChainFromLevel(StabChain<Telt> & R, StabChain<Telt> const& L, int const& lev)
{
  int nbLevL=L.stabilizer.size();
  int nbLevR=R.stabilizer.size();
  if (nbLevL != nbLevR) {
    std::cerr << "We should have nbLevL = nbLevR. Maybe wrong code here\n";
    std::cerr << "nbLEvL=" << nbLevL << " nbLevR=" << nbLevR << "\n";
    throw TerminalException{1};
  }
  int n=L.n;
  for (int eLev=lev; eLev<nbLevL; eLev++) {
    R.stabilizer[eLev].orbit = L.stabilizer[eLev].orbit;
    for (int i=0; i<n; i++) {
      int idx=L.stabilizer[eLev].transversal[i];
      int posPerm;
      if (idx == -1) {
	posPerm=-1;
      }
      else {
	Telt ePerm=L.labels[idx];
	posPerm = GetLabelIndex(R.labels, ePerm);
      }
      R.stabilizer[eLev].transversal[i] = posPerm;
    }
    R.stabilizer[eLev].genlabels.clear();
    for (int eVal : L.stabilizer[eLev].genlabels) {
      Telt ePerm=L.labels[eVal];
      int pos = GetLabelIndex(R.labels, ePerm);
      R.stabilizer[eLev].genlabels.push_back(pos);
    }
  }
}




template<typename Telt>
bool IsFixedStabilizer(StabChain<Telt> const& Stot, int const& eLev, int const& pnt)
{
  for (auto & posGen : Stot.stabilizer[eLev].genlabels) {
    if (pnt != PowAct(pnt, Stot.labels[posGen]))
      return false;
  }
  return true;
}


template<typename Telt>
Telt MinimalElementCosetStabChain(StabChain<Telt> const& Stot, Telt const& g)
{
  int nbLev=Stot.stabilizer.size();
  Telt gRet=g;
  for (int eLev=0; eLev<nbLev; eLev++) {
    if (Stot.stabilizer[eLev].genlabels.size() == 0)
      return gRet;
    int pMin=Stot.n + 1;
    for (auto & i : Stot.stabilizer[eLev].orbit) {
      int a=PowAct(i,gRet);
      if (a < pMin)
	pMin=a;
    }
    int bp=Stot.stabilizer[eLev].orbit[0];
    int pp=SlashAct(pMin, gRet);
    while (bp != pp) {
      int pos=Stot.stabilizer[eLev].transversal[pp];
      gRet=LeftQuotient(Stot.labels[pos], gRet);
      pp = SlashAct(pMin, gRet);
    }
  }
  return gRet;
}






template<typename Telt>
StabChain<Telt> HomomorphismMapping(StabChain<Telt> const& Stot, std::function<Telt(Telt const&)> const& f)
{
  Telt idMap = f(Stot.identity);
  int nMap=idMap.size();
  auto fVector =[&](std::vector<Telt> const& V) -> std::vector<Telt> {
    std::vector<Telt> Vret;
    for (auto & eElt : V)
      Vret.push_back(f(eElt));
    return Vret;
  };
  std::vector<Telt> labelsMap = fVector(Stot.labels);
  std::vector<StabLevel<Telt>> stabilizerMap;
  for (auto & eLevel : Stot.stabilizer) {
    StabLevel<Telt> eLevelMap{eLevel.transversal, eLevel.orbit, eLevel.genlabels, eLevel.cycles, fVector(eLevel.treegen), fVector(eLevel.treegeninv), fVector(eLevel.aux), eLevel.treedepth, eLevel.diam};
    stabilizerMap.push_back(eLevelMap);
  }
  return {nMap, idMap, Stot.UseCycle, labelsMap, stabilizerMap};
}




}


#endif
