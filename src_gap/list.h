#ifndef SRC_GAP_LIST_H
#define SRC_GAP_LIST_H

#include <vector>
#include <Boost_bitset.h>

template<typename T>
std::vector<T> Reversed(std::vector<T> const& eList)
{
  int len=eList.size();
  std::vector<T> retList(len);
  for (int u=0; u<len; u++) {
    int pos=len - 1 - u;
    retList[pos] = eList[u];
  }
  return retList;
}



std::vector<int> ClosedInterval(int const& begin, int const& end)
{
  std::vector<int> clInt;
  for (int u=begin; u<end; u++)
    clInt.push_back(u);
  return clInt;
}


Face BlistList(std::vector<int> const& list, std::vector<int> const& sub)
{
  int len=list.size();
  Face ret(len);
  for (auto & eVal : sub) {
    int pos=PositionVect(list, eVal);
    ret[pos]=true;
  }
  return ret;
}


int PositionNthTrueBlist(Face const& blist, int const& hpos)
{
  boost::dynamic_bitset<>::size_type b=blist.find_first();
  int epos=0;
  while (true) {
    if (b == boost::dynamic_bitset<>::npos)
      return -1;
    if (epos == hpos)
      return b;
    // We increase
    b=blist.find_next(b);
    epos++;
  }
}



std::vector<int> ListBlist(std::vector<int> const& list, Face const& blist)
{
  std::vector<int> ret;
  int len=list.size();
  for (int i=0; i<len; i++) {
    if (blist[i])
      ret.push_back(list[i]);
  }
  return ret;
}

int SizeBlist(Face const& blist)
{
  return blist.count();
}

bool IsSubsetBlist(Face const& a, Face const& b)
{
  int siz=a.size();
  for (int i=0; i<siz; i++) {
    if (b[i] == 1 && a[i] == 0)
      return false;
  }
  return true;
}



void UniteBlist(Face & a, Face const& b)
{
  int siz=a.size();
  for (int i=0; i<siz; i++) {
    if (b[i] == 1)
      a[i]=1;
  }
}



void IntersectBlist(Face & a, Face const& b)
{
  int siz=a.size();
  for (int i=0; i<siz; i++) {
    int val=0;
    if (a[i] == 1 && b[i] == 1)
      val=1;
    a[i]=val;
  }
}


void SubtractBlist(Face & a, Face const& b)
{
  int siz=a.size();
  for (int i=0; i<siz; i++) {
    if (b[i] == 1)
      a[i]=0;
  }
}



#endif
