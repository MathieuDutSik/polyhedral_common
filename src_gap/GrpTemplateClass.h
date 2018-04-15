#ifndef GRP_TEMPLATE_CLASS_H
#define GRP_TEMPLATE_CLASS_H



template<typename T>
struct {
  RemoveOneEntry(int const& idx)
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
  InsertOneEntry(T const& eObj, int const& ThePos)
  {
    int nbInd=ListIndices.size();
    int nbObj=stabilizer.size();
    if (nbObj <= nbInd) {
      ListObj.push_back(eObj);
      nbStab++;
    }
    std::vector<int> AttainedIndex(nbObj,0);
    for (auto & pos : ListIndices)
      AttainedIndex[pos]=1;
    int FirstFree=-1;
    for (int iObj=0; iObj<nbObj; iObj++)
      if (FirstFree == -1 && AttainedIndex[iObj] == 0)
	FirstFree = iObj;
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
    ListIndices=NewListIndices;
  }
  clear()
  {
    ListIndices.clear();
  }
  T& get(int const& idx)
  {
    int eIdx=ListIndices[idx];
    return ListObj[eIdx];
  }
  T get(int const& idx) const
  {
    int eIdx=ListIndices[idx];
    return ListObj[eIdx];
  }
private:
  std::vector<T> ListObj;
  std::vector<int> ListIndices;
}


#endif
