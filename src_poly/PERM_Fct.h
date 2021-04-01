#ifndef INCLUDE_PERM_FCT_H
#define INCLUDE_PERM_FCT_H



template<typename T, typename Tidx>
std::vector<Tidx> SortingPerm(std::vector<T> const & ListV)
{
  struct PairData {
    Tidx i;
    T x;
  };
  Tidx len=ListV.size();
  std::vector<PairData> ListPair(len);
  for (Tidx i=0; i<len; i++) {
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
  std::vector<Tidx> v(len);
  for (Tidx i=0; i<len; i++)
    v[i]=ListPair[i].i;
  return v;
}









#endif
