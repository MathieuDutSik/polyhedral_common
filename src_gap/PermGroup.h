#ifndef DEFINE_INCLUDE_PERM_GROUP
#define DEFINE_INCLUDE_PERM_GROUP

#include "Face_basic.h"


template<typename Telt>
int SmallestMovedPoint(std::vector<Telt> const& LGen)
{
  if (LGen.size() == 0)
    return -1;
  int n=LGen[0].size();
  for (int u=0; u<n; u++) {
    bool IsOK=false;
    for (auto & eGen : LGen)
      if (eGen.at(u) != u)
	IsOK=true;
    if (IsOK)
      return u;
  }
  return -1;
}


template<typename Telt>
std::vector<int> OrbitPerms(std::vector<Telt> const& gens, int const& d)
{
  std::vector<int> orb;
  int siz=gens[0].size();
  Face eFace(siz);
  auto InsertValue=[&](int const& val) -> void {
    orb.push_back(val);
    eFace[val]=1;
  };
  InsertValue(d);
  int posDone=0;
  while(true) {
    int posTot=orb.size();
    if (posTot == posDone)
      break;
    for (int u=posDone; u<posTot; u++) {
      int pnt=orb[u];
      for (auto & eGen : gens) {
	int img=PowAct(pnt, eGen);
	if (eFace[img] == 0)
	  InsertValue(img);
      }
    }
    posDone=posTot;
  }
  return orb;
}


template<typename Telt>
std::vector<std::vector<int>> OrbitsPerms(std::vector<Telt> const& gens, std::vector<int> const& D)
{
  std::vector<std::vector<int>> orbs;
  int max=gens[0].size();
  Face dom(max);
  for (auto & eV : D)
    dom[eV]=1;
  while(true) {
    if (dom.count() == 0)
      break;
    std::vector<int> orb;
    auto insert=[&](int const& eV) -> void {
      orb.push_back(eV);
      dom[eV]=0;
    };
    int fst=dom.find_first();
    insert(fst);
    int posDone=0;
    while(true) {
      int posTot=orb.size();
      if (posTot == posDone)
	break;
      for (int u=posDone; u<posTot; u++) {
	int pnt=orb[u];
	for (auto & eGen : gens) {
	  int img=PowAct(pnt, eGen);
	  if (dom[img] == 1)
	    insert(img);
	}
      }
      posDone=posTot;
    }
    orbs.push_back(orb);
  }
  return orbs;
}

template<typename Telt>
int SmallestMovedPointsPerms(std::vector<Telt> const& gens)
{
  int siz=0;
  for (int i=0; i<siz; i++) {
    for (auto & eGen : gens)
      if (PowAct(i, eGen) != i)
	return i;
  }
  return -1;
}

template<typename Telt>
int LargestMovedPointsPerms(std::vector<Telt> const& gens)
{
  int siz=0;
  for (int i=0; i<siz; i++) {
    int j=siz - 1 - i;
    for (auto & eGen : gens)
      if (PowAct(j, eGen) != j)
	return j;
  }
  return -1;
}



template<typename Telt>
std::vector<int> MovedPointsPerms(std::vector<Telt> const& gens)
{
  std::vector<int> ListMoved;
  int siz=gens[0].size();
  for (int i=0; i<siz; i++) {
    auto IsMoved=[&](int const& ePt) -> bool {
      for (auto & eGen : gens)
	if (PowAct(ePt, eGen) != ePt)
	  return true;
      return false;
    };
    if (IsMoved(i))
      ListMoved.push_back(i);
  }
  return ListMoved;
}




template<typename Telt>
Telt RestrictedPermNC(Telt const& x, std::vector<int> const& listRes)
{
  int n=x.size();
  std::vector<int> MapRev(n, -1);
  int nbRes=listRes.size();
  for (int iRes=0; iRes<nbRes; iRes++) {
    int ePt=listRes[iRes];
    MapRev[ePt] = iRes;
  }
  std::vector<int> eList(nbRes);
  for (int iRes=0; iRes<nbRes; iRes++) {
    int ePt=listRes[iRes];
    int ePtImg=x.at(ePt);
    int iResImg=MapRev[ePtImg];
    eList[iRes]=iResImg;
  }
  Telt eRet(eList);
  return eRet;
}

template<typename Telt, typename Tobj>
std::vector<Tobj> Orbit(std::vector<Telt> const& ListGen, Tobj const& x, std::function<Tobj(Tobj const&,Telt const&)> const& act)
{
  std::vector<Tobj> ListObj{x};
  std::vector<int> ListStat(0);
  std::set<Tobj> SetObj{x};
  while(true) {
    bool IsFinished=true;
    int len=ListObj.size();
    for (int u=0; u<len; u++) {
      if (ListStat[u] == 0) {
	IsFinished=false;
	ListStat[u]=1;
	for (auto & eElt : ListGen) {
	  Tobj eImg = act(ListObj[u], eElt);
	  if (SetObj.find(eImg) == SetObj.end()) {
	    ListObj.push_back(eImg);
	    ListStat.push_back(0);
	    SetObj.insert(eImg);
	  }
	}
      }
    }
    if (IsFinished)
      break;
  }
  return ListObj;
}

/*
template<typename Telt>
std::vector<int> OrbitPerms(std::vector<Telt> const& ListGen, int const & ePt)
{
  auto act=[](Telt const& u, int const& eVal) -> int {
    return u.at(eVal);
  };
  return Orbit(ListGen, ePt, act);
}
*/


template<typename Telt>
int CycleLength(Telt const& u, int const& x)
{
  int CycleLen=0;
  int xFirst=x;
  int xWork=x;
  while(true) {
    int xImg=u.at(xWork);
    CycleLen++;
    if (xImg == xFirst)
      break;
    xWork = xImg;
  }
  return CycleLen;
}


template<typename Telt>
std::vector<int> OnSets(std::vector<int> const& V, Telt const& u)
{
  std::vector<int> Vret;
  for (auto & ePt : V)
    Vret.push_back(u.at(ePt));
  sort(Vret.begin(), Vret.end());
  return Vret;
}

template<typename Telt>
Telt PowerGroupElement(Telt const& u, int const& n)
{
  if (n <= 0) {
    std::cerr << "We should have n >= 1\n";
    throw TerminalException{1};
  }
  Telt pow = u;
  for (int i=1; i<n; i++)
    pow = pow * u;
  return pow;
}




#endif
