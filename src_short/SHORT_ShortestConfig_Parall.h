// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_SHORT_SHORT_SHORTESTCONFIG_PARALL_H_
#define SRC_SHORT_SHORT_SHORTESTCONFIG_PARALL_H_

// clang-format off
#include "SHORT_ShortestConfig.h"
#include "Parallel_Classes.h"
#include <string>
#include <vector>
// clang-format on

template <typename T, typename Tint, typename Tgroup>
std::vector<MyMatrix<Tint>>
SHORT_SimplicialEnumeration(std::vector<MyMatrix<Tint>> const &ListSHVinp,
                            int const &NPROC, std::string const &TheMethod, std::ostream& os) {
  int nbEnt = ListSHVinp.size();
  int res_np = nbEnt % NPROC;
  int q_np = (nbEnt - res_np) / NPROC;
  std::vector<int> ListPos(NPROC + 1, 0);
  std::vector<int> ListSiz(NPROC, q_np);
  for (int i = 0; i < res_np; i++)
    ListSiz[i]++;
  for (int i = 0; i < NPROC; i++)
    ListPos[i + 1] = ListPos[i] + ListSiz[i];
  //
  bool eSave = false;
  bool eMemory = true;
  std::string ePrefix = "/irrelevant";
  std::function<std::optional<MyMatrix<Tint>>(SHVshortest<T, Tint> const &,
                                              SHVshortest<T, Tint> const &)>
      fEquiv =
          [&](SHVshortest<T, Tint> const &M1,
             SHVshortest<T, Tint> const &M2) -> std::optional<MyMatrix<Tint>> {
    return SHORT_TestEquivalence<T, Tint, Tgroup>(M1.SHV, M2.SHV, os);
  };
  std::function<int(SHVshortest<T, Tint> const &)> fSize =
      []([[maybe_unused]] SHVshortest<T, Tint> const &M) -> int { return 0; };
  FctsDataBank<SHVshortest<T, Tint>> recFct{fEquiv, fSize};
  DataBank<SHVshortest<T, Tint>> TheBank(eSave, eMemory, ePrefix, recFct);
  std::atomic<int> NbDone(0);
  std::condition_variable cv;
  std::mutex mtx_cv;
  auto TreatEntries = [&](int idx) -> void {
    for (int i = ListPos[idx]; i < ListPos[idx + 1]; i++) {
      std::vector<MyMatrix<Tint>> ListSpann =
          SHORT_SpannSimplicial<T, Tint, Tgroup>(ListSHVinp[i], ListSHVinp,
                                                 TheMethod, os);
      for (auto &eSpann : ListSpann) {
        SHVshortest<T, Tint> eEnt{eSpann};
        SHVinvariant<T, Tint> eInv = SHORT_Invariant<T, Tint>(eSpann);
        TheBank.InsertEntry(eEnt, eInv);
      }
    }
    NbDone++;
    cv.notify_one();
  };
  std::vector<std::thread> ListThr(NPROC);
  for (int iProc = 0; iProc < NPROC; iProc++)
    ListThr[iProc] = std::thread(TreatEntries, iProc);
  for (int iProc = 0; iProc < NPROC; iProc++)
    ListThr[iProc].detach();
  std::unique_lock<std::mutex> lk(mtx_cv);
  cv.wait(lk, [&] { return NbDone == NPROC; });
  //
  // Collecting output
  //
  int nbEntry = TheBank.GetNbEntry();
  std::vector<MyMatrix<Tint>> ListRet(nbEntry);
  for (int i = 0; i < nbEntry; i++)
    ListRet[i] = TheBank.GetEntry(i).SHV;
  return ListRet;
}





// clang-format off
#endif  // SRC_SHORT_SHORT_SHORTESTCONFIG_PARALL_H_
// clang-format on
