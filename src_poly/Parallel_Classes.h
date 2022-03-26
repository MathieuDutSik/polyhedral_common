#ifndef SRC_POLY_PARALLEL_CLASSES_H_
#define SRC_POLY_PARALLEL_CLASSES_H_

#include "ThreadManagement.h"
#include <utility>
#include <string>
#include <vector>
#include <set>

//
// Fundamental struct data types
//

template <typename Tequiv> struct EquivInfo {
  int iEntry;
  Tequiv TheEquiv;
  int nbEntryRelevant;
};

template <typename T> struct invariant_info {};

template <typename T> struct equiv_info {};

template <typename T> struct DataBank_ResultQuery {
  typedef typename equiv_info<T>::equiv_type Tequiv;
  bool test;
  T eEnt;
  Tequiv TheEquiv;
};

template <typename T> struct FctsDataBank {
  typedef typename equiv_info<T>::equiv_type Tequiv;
  std::function<std::optional<Tequiv>(T const &, T const &)> FctEquiv;
  std::function<int(T const &)> FctSize;
};

template <typename T> struct PairT_Tinv {
  typedef typename invariant_info<T>::invariant_type Tinv;
  T x;
  Tinv xInv;
};

// DataBank database for storing known results
//
// This can be for example:
// --- k-dimensional faces of a polytope as computed from the (k+1)-dimensional
// faces
// --- configuration of rank k shortest vectors from the configuration of rank
// k-1
// --- database of known dual descriptions.
//
// Functionality
// --- Use of invariant for faster check
// --- We may or may not have additional information in the data (like dual
// description)
//     but this is up to you in the design of your type T.
// --- There is not entry save for whether we did the check or not.
//
// We maintain T separate from Tinv instead of merging them because we want to
// consider the situation where not all data are in memory
//
// For example list of known double descriptions.
template <typename T> struct DataBank {
  typedef typename equiv_info<T>::equiv_type Tequiv;
  typedef typename invariant_info<T>::invariant_type Tinv;

public:
  DataBank() {
    std::cerr << "We should never init the DataBank with default constructor\n";
    throw TerminalException{1};
  }
  DataBank(bool const &eSave, bool const &eMemory, std::string const &ePrefix,
           FctsDataBank<T> const &recFct)
      : FctEquiv(recFct.FctEquiv), FctSize(recFct.FctSize) {
    MinSize = -1;
    IsSaving = eSave;
    FullDataInMemory = eMemory;
    if (!FullDataInMemory && !IsSaving) {
      std::cerr << "We have FullDataInMemory = false\n";
      std::cerr << "and IsSaving = false\n";
      throw TerminalException{1};
    }
    ThePrefix = ePrefix;
    nbEntry = 0;
    if (IsSaving) {
      bool test = FILE_CheckPrefix(ePrefix);
      if (!test) {
        std::cerr << "ePrefix=" << ePrefix << "\n";
        std::cerr << "which is not correct for DataBank\n";
        throw TerminalException{1};
      }
      while (true) {
        std::string FileSaveDAT =
            ThePrefix + "BankEntryDAT_" + IntToString(nbEntry + 1);
        std::string FileSaveINV =
            ThePrefix + "BankEntryINV_" + IntToString(nbEntry + 1);
        if (IsExistingFile(FileSaveDAT) && IsExistingFile(FileSaveINV)) {
          std::ifstream isT(FileSaveDAT);
          T eEnt;
          isT >> eEnt;
          int eSize = FctSize(eEnt);
          if (MinSize == -1 || eSize < MinSize)
            MinSize = eSize;
          if (FullDataInMemory)
            ListDat.emplace_back(std::move(eEnt));
          //
          std::ifstream isI(FileSaveINV);
          Tinv eInv;
          isI >> eInv;
          ListInv.push_back(eInv);
          //
          nbEntry++;
        } else {
          break;
        }
      }
    }
  }

  // no copy
  DataBank(const DataBank &) = delete;

  // no assign
  DataBank &operator=(const DataBank &) = delete;

  // no move
  DataBank(DataBank &&) = delete;

  EquivInfo<Tequiv> IsPresentNoLock(int const &iEntryStart, T const &eEnt,
                                    Tinv const &eInv) const {
    for (int iEntry = iEntryStart; iEntry < nbEntry; iEntry++) {
      Tinv fInv = ListInv[iEntry];
      if (fInv == eInv) {
        T fEnt = GetEntry(iEntry);
        std::optional<Tequiv> eEquiv = FctEquiv(fEnt, eEnt);
        if (eEquiv)
          return {iEntry, *eEquiv, nbEntry};
      }
    }
    return {-1, {}, nbEntry};
  }
  DataBank_ResultQuery<T> ProcessRequest(T const &eEnt, Tinv const &eInv,
                                         std::ostream &os) {
    os << "DataBank.ProcessRequest, step 2\n";
    EquivInfo<Tequiv> eEquiv = IsPresentNoLock(0, eEnt, eInv);
    os << "DataBank.ProcessRequest, step 3\n";
    if (eEquiv.iEntry != -1) {
      T fEnt = GetEntry(eEquiv.iEntry);
      return {true, std::move(fEnt), eEquiv.TheEquiv};
    }
    os << "DataBank.ProcessRequest, step 4\n";
    return {false, {}, {}};
  }
  T GetEntry(int const &i) const {
    if (FullDataInMemory)
      return ListDat[i];
    std::string FileSaveDAT =
        ThePrefix + "BankEntryDAT_" + IntToString(nbEntry + 1);
    std::ifstream isT(FileSaveDAT);
    T eEnt;
    isT >> eEnt;
    return eEnt;
  }
  int GetNbEntry() const { return nbEntry; }
  int GetMinSize() const { return MinSize; }
  std::vector<T> GetAllEntry() const {
    std::vector<T> ListEntry(nbEntry);
    for (int iEntry = 0; iEntry < nbEntry; iEntry++)
      ListEntry[iEntry] = GetEntry(iEntry);
    return ListEntry;
  }
  void InsertEntry(T const &eEnt, Tinv const &eInv) {
    EquivInfo<Tequiv> eEquiv = IsPresentNoLock(0, eEnt, eInv);
    if (eEquiv.iEntry == -1) {
      std::lock_guard<std::mutex> lk(mul);
      EquivInfo<Tequiv> fEquiv =
          IsPresentNoLock(eEquiv.nbEntryRelevant, eEnt, eInv);
      if (fEquiv.iEntry != -1)
        return;
      int eSize = FctSize(eEnt);
      if (IsSaving) {
        std::string FileSaveDAT =
            ThePrefix + "BankEntryDAT_" + IntToString(nbEntry);
        std::ofstream osT(FileSaveDAT);
        osT << eEnt;
        //
        std::string FileSaveINV =
            ThePrefix + "BankEntryINV_" + IntToString(nbEntry);
        std::ofstream osI(FileSaveINV);
        osI << eInv;
      }
      if (eSize < MinSize)
        MinSize = eSize;
      if (FullDataInMemory)
        ListDat.emplace_back(std::move(eEnt));
      ListInv.push_back(eInv);
      nbEntry = ListInv.size();
    }
  }

private:
  std::function<std::optional<Tequiv>(T const &, T const &)> FctEquiv;
  std::function<int(T const &)> FctSize;
  bool IsSaving;
  bool FullDataInMemory;
  std::string ThePrefix;
  std::mutex mul;
  tbb::concurrent_vector<T> ListDat;
  tbb::concurrent_vector<Tinv> ListInv;
  int nbEntry;
  int MinSize;
};

//
// Balinski business
//

struct TrivialBalinski {
  bool final = true;
  bool IsComplete = false;
};

std::ostream &operator<<(std::ostream &os, TrivialBalinski const &obj) {
  os << "TrivialBalinski : final =" << obj.final
     << " IsComplete=" << obj.IsComplete << "\n";
  return os;
}

template <typename T> struct balinski_info {
  typedef TrivialBalinski balinski_type;
};

// We want to allow the storing of data on disk.
// ---This implies that we must make the comparison function on
//    invariants, since this is the only thing that is in memory
// ---The right structure for the invariant is multimap.
//    Problematic is having a parallel implementation in tbb.
// ---The testing of the completeness via Balinski like theorems.
//    This is especially complex. It requires that a function be sent
//    to the code when initializing it. But which one?
//    We want the function to be such that:
//    ---It returns true/false only
//    ---It is fast.
//    ---It ends after conclusion is reached which should be fast in most
//       cases. So, a single while loop.
//    From those requirements, it looks like there is no other way than
//    doing it all complicated with a Tbalinski type and a function
//    that upgrades the Tbalinski and adds to it the final status.
//    Final status should be directly readable on the type, which should be a
//    simple struct.
// ---The invariant and equivalence functions have to be provided as input.
//    We cannot use a template enable_if formalism since things may vary too
//    much
//
// Since the comparison function must use the data for doing the comparison.
// we have no choice for a first version to have all data in memory.
// And this implies that the object type themselves will also contain
// So, we may have Tinv
//
// This is for Graph Traversal algorithm.
// Again the same Hack, we do not want to define constructor ()
// So, we have to allow this.
template <typename T> struct NewEnumerationWork {
  typedef typename invariant_info<T>::invariant_type Tinv;
  typedef typename balinski_info<T>::balinski_type Tbalinski;
  typedef typename equiv_info<T>::equiv_type Tequiv;

public:
  // no copy
  NewEnumerationWork(const NewEnumerationWork &) = delete;

  // no assign
  NewEnumerationWork &operator=(const NewEnumerationWork &) = delete;

  // no move
  NewEnumerationWork(NewEnumerationWork &&) = delete;

  // default constructor
  NewEnumerationWork() = delete;

  // non-default constructor
  NewEnumerationWork(
      bool const &eSave, bool const &eMemory, std::string const &ePrefix,
      std::function<bool(Tinv const &, Tinv const &)> const &eCompFCT,
      std::function<void(Tbalinski &, T const &, Tinv const &,
                         std::ostream &)> const &eUP,
      std::function<std::optional<Tequiv>(T const &, T const &)> const &fEquiv,
      std::ostream &os)
      : ListPendingCases(
            std::set<int, std::function<bool(int const &, int const &)>>(
                [&](int const &i, int const &j) -> bool {
                  if (eCompFCT(ListInv[i], ListInv[j]))
                    return true;
                  if (eCompFCT(ListInv[j], ListInv[i]))
                    return false;
                  return i < j;
                })),
        UpgradeBalinskiStat(eUP), TestEquivalence(fEquiv) {
    IsSaving = eSave;
    FullDataInMemory = eMemory;
    ThePrefix = ePrefix;
    nbEntry = 0;
    nbEntryTreated = 0;
    if (!FullDataInMemory && !IsSaving) {
      std::cerr << "We have FullDataInMemory == false\n";
      std::cerr << "and IsSaving == false\n";
      std::cerr << "This simply cannot work\n";
      throw TerminalException{1};
    }
    // Template functionology
    //    auto comp=[&](int const& i, int const& j) -> bool {return
    //    CompFCT(ListInv[i],ListInv[j]);}; ListPendingCases=std::set<int,
    //    decltype(comp)> (comp); UpgradeBalinskiStat=eUP;
    //
    if (IsSaving) {
      bool test = FILE_CheckPrefix(ePrefix);
      if (!test) {
        std::cerr << "ePrefix=" << ePrefix << "\n";
        std::cerr << "which is not correct for NewEnumerationWork\n";
        throw TerminalException{1};
      }
      while (true) {
        std::string FileSaveDAT =
            ThePrefix + "OrbitDAT_" + IntToString(nbEntry);
        std::string FileSaveINV =
            ThePrefix + "OrbitINV_" + IntToString(nbEntry);
        std::string FileSaveSTA =
            ThePrefix + "OrbitSTA_" + IntToString(nbEntry);
        if (IsExistingFile(FileSaveDAT) && IsExistingFile(FileSaveINV) &&
            IsExistingFile(FileSaveSTA)) {
          if (FullDataInMemory) {
            std::ifstream isT(FileSaveDAT);
            T eEnt;
            isT >> eEnt;
            ListDat.push_back(eEnt);
          }
          //
          std::ifstream isI(FileSaveINV);
          Tinv eInv;
          isI >> eInv;
          ListInv.push_back(eInv);
          //
          std::ifstream isS(FileSaveSTA);
          int eSta;
          isS >> eSta;
          ListSta.push_back(eSta);
          nbEntryTreated += eSta;
          if (eSta == 0)
            ListPendingCases.insert(nbEntry);
          //
          nbEntry++;
        } else {
          break;
        }
      }
    }
    UpdateCompleteStatusNoLock(os);
  }
  ~NewEnumerationWork() = default;
  void FuncClear() {
    if (IsSaving) {
      for (int iEntry = 0; iEntry < nbEntry; iEntry++) {
        std::string FileSaveDAT = ThePrefix + "OrbitDAT_" + IntToString(iEntry);
        RemoveFileIfExist(FileSaveDAT);
        //
        std::string FileSaveINV = ThePrefix + "OrbitINV_" + IntToString(iEntry);
        RemoveFileIfExist(FileSaveINV);
        //
        std::string FileSaveSTA = ThePrefix + "OrbitSTA_" + IntToString(iEntry);
        RemoveFileIfExist(FileSaveSTA);
      }
    }
  }

  EquivInfo<Tequiv> IsPresentNoLock(int const &iEntryStart,
                                    PairT_Tinv<T> const &eRec,
                                    std::ostream &os) {
    //    os << "Start IsPresentNoLock\n";
    //    os << "eInv=" << eInv << "\n";
    //    os << "nbEntry=" << nbEntry << "\n";
    for (int iEntry = iEntryStart; iEntry < nbEntry; iEntry++) {
      if (ListInv[iEntry] == eRec.xInv) {
        //	os << "IsPresentNoLock, ListInv[iEntry]=" << ListInv[iEntry] <<
        //"\n";
        T fEnt = GetRepresentative(iEntry);
        std::optional<Tequiv> eEquiv = TestEquivalence(fEnt, eRec.x);
        if (eEquiv) {
          //	  os << "Before exit IsPresentNoLock. Find Isomorphism\n";
          return {iEntry, *eEquiv, nbEntry};
        }
      }
    }
    //    os << "Before exit IsPresentNoLock. Did NOT Find Isomorphism\n";
    return {-1, {}, nbEntry};
  }

  DataBank_ResultQuery<T> ProcessRequest(PairT_Tinv<T> const &eRec,
                                         std::ostream &os) {
    EquivInfo<Tequiv> eEquiv = IsPresentNoLock(0, eRec, os);
    if (eEquiv.iEntry != -1) {
      T fEntry = GetRepresentative(eEquiv.iEntry);
      return {true, fEntry, eEquiv.TheEquiv};
    }
    return {false, {}, {}};
  }

  int InsertEntryLock(int const &eStart, PairT_Tinv<T> const &eRec,
                      std::ostream &os) {
    std::lock_guard<std::mutex> lk(mul);
    EquivInfo<Tequiv> fEquiv = IsPresentNoLock(eStart, eRec, os);
    if (fEquiv.iEntry != -1)
      return fEquiv.iEntry;
    if (FullDataInMemory)
      ListDat.push_back(eRec.x);
    ListInv.push_back(eRec.xInv);
    nbEntry = ListInv.size();
    int ThePos = nbEntry - 1;
    ListPendingCases.insert(ThePos);
    int eSta = 0;
    ListSta.push_back(eSta);
    //    os << "|ListPendingCases|=" << ListPendingCases.size() << " Now
    //    |ListDat|=" << ListDat.size() << "\n";
    int RetVal = -nbEntry;
    return RetVal;
  }

  int GetNbEntry() const { return nbEntry; }

  // positive value means that the orbit was already present (and value is
  // number of orbit) -1 means that the orbit is new
  int InsertEntry(PairT_Tinv<T> const &eRec, std::ostream &os) {
    EquivInfo<Tequiv> eEquiv = IsPresentNoLock(0, eRec, os);
    if (eEquiv.iEntry == -1) {
      int iEntry = InsertEntryLock(eEquiv.nbEntryRelevant, eRec, os);
      if (iEntry >= 0)
        return iEntry;
      int ThePos = -1 - iEntry;
      if (IsSaving) {
        std::string FileSaveDAT = ThePrefix + "OrbitDAT_" + IntToString(ThePos);
        std::ofstream osT(FileSaveDAT);
        osT << eRec.x;
        //
        std::string FileSaveINV = ThePrefix + "OrbitINV_" + IntToString(ThePos);
        std::ofstream osI(FileSaveINV);
        osI << eRec.xInv;
        //
        int eSta = 0;
        std::string FileSaveSTA = ThePrefix + "OrbitSTA_" + IntToString(ThePos);
        std::ofstream osS(FileSaveSTA);
        osS << eSta;
      }
      return -1;
    }
    return eEquiv.iEntry;
  }
  // The lock is actually put at some other place
  void UpdateCompleteStatusNoLock(std::ostream &os) {
    //    os << "Begin UpdateCompleteStatusNoLock\n";
    if (nbEntryTreated == 0) {
      os << "Leaving because nbEntryTreated=0\n";
      IsComplete = false;
      return;
    }
    // The Balinski iterate over the cases.
    // This is expensive, but we should conclude in not too many steps.
    //    os << "|ListPendingCases|=" << ListPendingCases.size() << "
    //    |ListDat|=" << ListDat.size() << "\n";
    Tbalinski eStat;
    for (auto &iEntry : ListPendingCases) {
      T eEnt = GetRepresentative(iEntry);
      Tinv eInv = ListInv[iEntry];
      //      os << "Before UpgradeBalinskiStat\n";
      UpgradeBalinskiStat(eStat, eEnt, eInv, os);
      //      os << "After UpgradeBalinskiStat\n";
      if (eStat.final) {
        IsComplete = eStat.IsComplete;
        os << "Leaving with IsComplete=" << IsComplete << "\n";
        return;
      }
    }
    os << eStat;
    os << "Leaving with IsComplete=true\n";
    IsComplete = true;
  }

  bool GetCompleteStatus() const {
    //    std::lock_guard<std::mutex> lk(mul);
    return IsComplete;
  }

  int GetNonTreatedOrbit(std::ostream &os) {
    std::lock_guard<std::mutex> lk(mul);
    os << "|ListPendingCases|=" << ListPendingCases.size()
       << " |ListDat|=" << ListDat.size()
       << " nbEntryTreated=" << nbEntryTreated << "\n";
    int eEntryFound = -1;
    auto iter = ListPendingCases.begin();
    while (true) {
      if (iter == ListPendingCases.end())
        break;
      int eEntry = *iter;
      std::set<int>::iterator iterB = ListUnderConsideration.find(eEntry);
      if (iterB == ListUnderConsideration.end()) {
        eEntryFound = eEntry;
        break;
      }
      iter++;
    }
    if (eEntryFound == -1)
      return -1;
    size_t reply = ListPendingCases.erase(eEntryFound);
    if (reply == 0) {
      std::cerr << "We should have found eEntryFound = " << eEntryFound << "\n";
      std::cerr << "Bug in the erase setting\n";
      throw TerminalException{1};
    }
    ListUnderConsideration.insert(eEntryFound);
    return eEntryFound;
  }

  void SetEntryAsDone(int const &iEntry, std::ostream &os) {
    int eSta = 1;
    if (IsSaving) {
      std::string FileSave = ThePrefix + "OrbitSTA_" + IntToString(iEntry);
      std::ofstream osSave(FileSave);
      osSave << eSta;
    }
    //
    std::lock_guard<std::mutex> lk(mul);
    auto iter = ListUnderConsideration.find(iEntry);
    if (iter == ListUnderConsideration.end()) {
      std::cerr << "iEntry=" << iEntry << "\n";
      std::cerr << "ListUnderConsideration does not have iEntry\n";
      throw TerminalException{1};
    }
    size_t reply = ListUnderConsideration.erase(iEntry);
    if (reply == 0) {
      std::cerr << "Error in erase iEntry=" << iEntry
                << " from ListUnderConsideration\n";
      std::cerr << "Bug in the erase setting\n";
      throw TerminalException{1};
    }
    ListSta[iEntry] = eSta;
    nbEntryTreated++;
    UpdateCompleteStatusNoLock(os);
  }

  std::vector<int> GetListPendingCases() const {
    std::lock_guard<std::mutex> lk(mul);
    std::vector<int> eVect;
    for (auto &eCase : ListPendingCases)
      eVect.push_back(eCase);
    return eVect;
  }

  Tinv GetInvariant(int const &i) const { return ListInv[i]; }

  T GetRepresentative(int const &i) const {
    if (FullDataInMemory) {
      return ListDat[i];
    } else {
      std::string FileSaveDAT = ThePrefix + "OrbitDAT_" + IntToString(i);
      if (!IsExistingFile(FileSaveDAT)) {
        std::cerr << "Error: trying to access orbit nr: " << i << "\n";
        std::cerr << "But file = " << FileSaveDAT << "\n";
        std::cerr << "is not existent\n";
        throw TerminalException{1};
      }
      std::ifstream isT(FileSaveDAT);
      T eEnt;
      isT >> eEnt;
      return eEnt;
    }
  }

  bool IsStuck() const {
    //    std::lock_guard<std::mutex> lk(mul);
    if (ListUnderConsideration.size() == 0 && ListPendingCases.size() == 0)
      return false;
    if (!IsComplete && ListPendingCases.size() == 0)
      return true;
    return false;
  }

  std::vector<T> GetListRepr() const {
    std::vector<T> ListRepr(nbEntry);
    for (int iEntry = 0; iEntry < nbEntry; iEntry++)
      ListRepr[iEntry] = GetRepresentative(iEntry);
    return ListRepr;
  }

private:
  std::mutex mul;
  bool IsSaving;
  std::string ThePrefix;
  bool FullDataInMemory;
  tbb::concurrent_vector<T> ListDat;
  tbb::concurrent_vector<Tinv> ListInv;
  std::vector<int> ListSta;
  int nbEntry;
  int nbEntryTreated;
  std::set<int, std::function<bool(int const &, int const &)>> ListPendingCases;
  std::set<int> ListUnderConsideration;
  bool IsComplete;
  std::function<void(Tbalinski &, T const &, Tinv const &, std::ostream &)>
      UpgradeBalinskiStat;
  std::function<std::optional<Tequiv>(T const &, T const &)> TestEquivalence;
};

#endif //  SRC_POLY_PARALLEL_CLASSES_H_
