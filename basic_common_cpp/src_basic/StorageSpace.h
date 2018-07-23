#ifndef STORAGESPACE_H
#define STORAGESPACE_H

#include <vector>
#include <utility>


template<typename T>
struct StorageSpaceLastN {
public:
  // no copy
  StorageSpaceLastN(const StorageSpaceLastN<T> &) = delete;

  // no assign
  StorageSpaceLastN& operator=(const StorageSpaceLastN<T> &) = delete;

  // no move
  StorageSpaceLastN(StorageSpaceLastN<T> &&) = delete;

  StorageSpaceLastN()
  {
    N = 0;
    lastpos = -1;
    nbCorrect = 0;
  }
  
  StorageSpaceLastN(int const& _N)
  {
    N = _N;
    lastpos = -1;
    nbCorrect = 0;
    ListVal.resize(N);
  }
  std::vector<T> RetrieveLastRelevantValues() const
  {
    std::vector<T> ListRet(nbCorrect);
    int ePos=lastpos;
    for (int i=0; i<nbCorrect; i++) {
      ListRet[i]=ListVal[ePos];
      ePos--;
      if (ePos == -1)
	ePos = N-1;
    }
    return ListRet;
  }
  void Resize(int const& NewN) 
  {
    std::vector<T> NewListVal(NewN);
    int SizInit = std::min(nbCorrect, std::min(NewN, N));
    int pos = lastpos;
    int NewLastpos = SizInit-1;
    for (int i=0; i<SizInit; i++) {
      int Npos = SizInit - 1 - i;
      NewListVal[Npos] = ListVal[pos];
      pos--;
      if (pos == -1)
	pos = N-1;
    }
    ListVal = NewListVal;
    N = NewN;
    lastpos = NewLastpos;
    nbCorrect = SizInit;
  }
  void InsertOneValue(T const& eVal)
  {
    if (N == 0)
      return;
    nbCorrect++;
    if (nbCorrect > N)
      nbCorrect = N;
    lastpos++;
    if (lastpos == N)
      lastpos = 0;
    ListVal[lastpos] = eVal;
  }
  int GetN() const
  {
    return N;
  }
  int GetNbCorrect() const
  {
    return nbCorrect;
  }
  T GetLast() const
  {
    return ListVal[lastpos];
  }
  T GetPrevious(int const& pos) const
  {
    int ThePos=lastpos - pos;
    if (ThePos < 0)
      ThePos += N;
    return ListVal[ThePos];
  }
private:
  std::vector<T> ListVal;
  int N;
  int lastpos;
  int nbCorrect;
};


#endif
