#ifndef INCLUDE_PERM_FCT_H
#define INCLUDE_PERM_FCT_H



template<typename T>
std::vector<int> SortingPerm(std::vector<T> const & ListV)
{
  struct PairData {
    std::size_t i;
    T x;
  };
  std::size_t len=ListV.size();
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
  std::vector<int> v(len);
  for (std::size_t i=0; i<len; i++)
    v[i]=ListPair[i].i;
  return v;
}









#endif
