#ifndef DEFINE_COMB_VECTORS_H
#define DEFINE_COMB_VECTORS_H


#include <vector>
#include <set>


template<typename T>
std::vector<T> UnionVect(std::vector<T> const& V1, std::vector<T> const& V2)
{
  std::set<T> eSet;
  for (auto & eVal : V1)
    eSet.insert(eVal);
  for (auto & eVal : V2)
    eSet.insert(eVal);
  std::vector<T> eV;
  for (auto & eVal : eSet)
    eV.push_back(eVal);
  return eV;
}

// The difference V1 - V2
template<typename T>
std::vector<T> DifferenceVect(std::vector<T> const& V1, std::vector<T> const& V2)
{
  std::set<T> eSet2;
  for (auto & eVal : V2)
    eSet2.insert(eVal);
  std::vector<T> eV;
  for (auto & eVal : V1) {
    typename std::set<T>::iterator iter=eSet2.find(eVal);
    if (iter == eSet2.end())
      eV.push_back(eVal);
  }
  return eV;
}

template<typename T>
void AppendVect(std::vector<T> & V1, std::vector<T> const& V2)
{
  for (auto & eVal : V2)
    V1.push_back(eVal);
}




template<typename T>
bool IsSubset(std::vector<T> const& S1, std::vector<T> const& S2)
{
  for (auto & eVal : S2) {
    int pos=PositionVect(S1, eVal);
    if (pos == -1)
      return false;
  }
  return true;
}


template<typename T>
std::vector<T> VectorAsSet(std::vector<T> const& V)
{
  std::set<T> eSet;
  for (auto & eVal : V)
    eSet.insert(eVal);
  std::vector<T> eV;
  for (auto & eVal : eSet)
    eV.push_back(eVal);
  return eV;
}


template<typename T>
std::vector<T> IntersectionVect(std::vector<T> const& V1, std::vector<T> const& V2)
{
  std::set<T> eSet2;
  for (auto & eVal : V2)
    eSet2.insert(eVal);
  std::vector<T> eV;
  for (auto & eVal : V1) {
    typename std::set<T>::iterator iter=eSet2.find(eVal);
    if (iter != eSet2.end())
      eV.push_back(eVal);
  }
  return eV;
}






#endif
