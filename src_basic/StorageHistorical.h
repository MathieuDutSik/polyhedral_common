#ifndef STORAGEHISTORICAL_H
#define STORAGEHISTORICAL_H

#include <vector>
#include <utility>




template<typename T>
struct StorageHistorical {
  struct FullRelInfo {
    bool status;
    unsigned long long eTime;
    T eVal;
  };
public:
  // no copy
  StorageHistorical(const StorageSpaceLastN<T> &) = delete;

  // no assign
  StorageHistorical& operator=(const StorageHistorical<T> &) = delete;

  // no move
  StorageHistorical(StorageHistorical<T> &&) = delete;

  StorageHistorical()
  {
    ThePeriod = -1;
    siz = 0;
  }
  
  StorageHistorical(int const& _ThePeriod)
  {
    ThePeriod = _ThePeriod;
  }
  
  void SetPeriod(int const& _ThePeriod)
  {
    ThePeriod = _ThePeriod;
  }
  
  std::vector<T> RetrieveLastRelevantValues() const
  {
    std::vector<T> ListVal;
    for (auto & eFull : ListFull)
      if (eFull.status)
	ListVal.push_back(eFull.eVal);
    return ListVal;
  }
  
  void InsertOneValue(unsigned long long const& eTime, T const& eVal)
  {
    if (ThePeriod < 0) {
      return;
    }
    bool HasFoundPlace=false;
    for (int i=0; i<siz; i++) {
      int deltaTime=int(eTime - ListFull[i].eTime);
      if (deltaTime > ThePeriod)
	ListFull[i].status = false;
      if (!HasFoundPlace && !ListFull[i].status) {
	ListFull[i] = {true, eTime, eVal};
	HasFoundPlace = true;
      }
    }
    if (!HasFoundPlace) {
      ListFull.push_back({true, eTime, eVal});
      siz++;
    }
  }
private:
  std::vector<FullRelInfo> ListFull;
  int siz;
  int ThePeriod;
};


#endif
