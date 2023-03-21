// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_GRP_DOUBLECOSET_H_
#define SRC_POLY_GRP_DOUBLECOSET_H_

#include "GRP_GroupFct.h"
#include "WeightMatrix.h"
#include <limits>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

//
// Double Coset Computation
//

template <typename Tgroup, typename Tidx_value>
vectface DoubleCosetDescription_Representation(
    Tgroup const &BigGRP, Tgroup const &SmaGRP,
    WeightMatrix<true, int, Tidx_value> const &WMat, Face const &eList) {
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  std::vector<Telt> ListGen = BigGRP.GeneratorsOfGroup();
  Tgroup TheStab = BigGRP.Stabilizer_OnSets(eList);
  Tint TotalSize = BigGRP.size() / TheStab.size();
  //
  struct Local {
    int status;
    Face eFace;
    size_t eInv;
  };
  Tint SizeGen = 0;
  std::vector<Local> ListLocal;
  auto DoubleCosetInsertEntry = [&](Face const &testList) -> void {
    size_t eInv = GetLocalInvariantWeightMatrix(WMat, testList);
    for (auto const &fLocal : ListLocal) {
      std::optional<Telt> test =
          SmaGRP.RepresentativeAction_OnSets(fLocal.eFace, testList);
      if (test)
        return;
    }
    ListLocal.push_back({0, testList, eInv});
    Tgroup fStab = SmaGRP.Stabilizer_OnSets(testList);
    Tint OrbSizeSma = SmaGRP.size() / fStab.size();
    SizeGen += OrbSizeSma;
  };
  DoubleCosetInsertEntry(eList);
  while (true) {
    bool DoSomething = false;
    size_t nbLocal = ListLocal.size();
    for (size_t iLocal = 0; iLocal < nbLocal; iLocal++)
      if (ListLocal[iLocal].status == 0) {
        ListLocal[iLocal].status = 1;
        DoSomething = true;
        Face eFace = ListLocal[iLocal].eFace;
        for (auto const &eGen : ListGen) {
          Face eNewList = OnFace(eFace, eGen);
          DoubleCosetInsertEntry(eNewList);
        }
      }
    if (!DoSomething)
      break;
  }
  vectface ListListSet;
  for (auto &eRec : ListLocal)
    ListListSet.push_back(eRec.eFace);
  if (SizeGen == TotalSize)
    return ListListSet;
  vectface PartialOrbit = std::move(ListListSet);
  auto IsPresent = [&](Face const &testList) -> bool {
    for (auto &fList : PartialOrbit)
      if (fList == testList)
        return true;
    return false;
  };
  size_t pos_start = 0;
  while (true) {
    size_t n_orb = PartialOrbit.size();
    for (size_t i_orb = pos_start; i_orb < n_orb; i_orb++) {
      for (auto &eGen : ListGen) {
        Face eNewList = OnFace(PartialOrbit[i_orb], eGen);
        if (!IsPresent(eNewList)) {
          PartialOrbit.push_back(eNewList);
          DoubleCosetInsertEntry(eNewList);
          if (SizeGen == TotalSize) {
            vectface ListListFin;
            for (auto &eRec : ListLocal)
              ListListFin.push_back(eRec.eFace);
            return ListListFin;
          }
        }
      }
    }
    pos_start = n_orb;
  }
  std::cerr << "Likely not reachable stage\n";
  throw TerminalException{1};
}

template <typename Tgroup>
vectface DoubleCosetDescription_Canonic(Tgroup const &BigGRP,
                                        Tgroup const &SmaGRP, Face const &eList,
                                        std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  std::vector<Telt> ListGen = BigGRP.GeneratorsOfGroup();
  Tgroup TheStab = BigGRP.Stabilizer_OnSets(eList);
  Tint TotalSize = BigGRP.size() / TheStab.size();
  //
  Tint SizeGen = 0;
  auto IncreaseSize = [&](Face const &eList) -> void {
    Tgroup fStab = SmaGRP.Stabilizer_OnSets(eList);
    Tint OrbSizeSma = SmaGRP.size() / fStab.size();
    SizeGen += OrbSizeSma;
  };
  std::unordered_set<Face> SetFace;
  vectface CurrList(BigGRP.n_act());
  auto DoubleCosetInsertEntry_first = [&](Face const &testList) -> void {
    Face faceCan = SmaGRP.CanonicalImage(testList);
    if (SetFace.count(faceCan) > 0)
      return;
    CurrList.push_back(faceCan);
    SetFace.insert(faceCan);
    IncreaseSize(faceCan);
  };
  auto DoubleCosetInsertEntry_second = [&](Face const &testList) -> void {
    Face faceCan = SmaGRP.CanonicalImage(testList);
    if (SetFace.count(faceCan) > 0)
      return;
    SetFace.insert(faceCan);
    IncreaseSize(faceCan);
  };
  DoubleCosetInsertEntry_first(eList);
  Face eFaceImg(BigGRP.n_act());
  while (true) {
    if (CurrList.size() == 0)
      break;
    Face eFace = CurrList.pop();
    for (auto const &eGen : ListGen) {
      OnFace_inplace(eFaceImg, eFace, eGen);
      DoubleCosetInsertEntry_first(eFaceImg);
    }
  }
  vectface ListListSet(BigGRP.n_act());
  for (auto &eFace : SetFace)
    ListListSet.push_back(eFace);
  if (SizeGen == TotalSize)
    return ListListSet;
  //  os << "After Iteration loop SizeGen=" << SizeGen << " TotalSize=" <<
  //  TotalSize << "\n";
  std::unordered_set<Face> PartialOrbit = SetFace;
  while (true) {
    Face eFace = ListListSet.pop();
    for (auto &eGen : ListGen) {
      OnFace_inplace(eFaceImg, eFace, eGen);
      if (PartialOrbit.count(eFaceImg) == 0) {
        PartialOrbit.insert(eFaceImg);
        ListListSet.push_back(eFaceImg);
        DoubleCosetInsertEntry_second(eFaceImg);
        if (SizeGen == TotalSize) {
          vectface ListListFin(BigGRP.n_act());
          for (auto &eFace : SetFace)
            ListListFin.push_back(eFace);
          return ListListFin;
        }
      }
    }
  }
  os << "Likely not reachable stage\n";
  throw TerminalException{1};
}

template <typename Tgroup>
vectface DoubleCosetDescription_Exhaustive(Tgroup const &BigGRP,
                                           Tgroup const &SmaGRP,
                                           Face const &eList) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = BigGRP.n_act();
  std::vector<Telt> LGenBig = BigGRP.GeneratorsOfGroup();
  //
  vectface vf(n);
  std::vector<uint8_t> status;
  std::unordered_set<Face> SetFace;
  auto f_insert = [&](const Face &f) -> void {
    if (SetFace.find(f) == SetFace.end()) {
      vf.push_back(f);
      status.push_back(0);
      SetFace.insert(f);
    }
  };
  size_t miss_val = std::numeric_limits<size_t>::max();
  auto get_undone = [&]() -> size_t {
    for (size_t i = 0; i < status.size(); i++)
      if (status[i] == 0)
        return i;
    return miss_val;
  };
  f_insert(eList);
  Face eFaceImg(n);
  while (true) {
    size_t pos = get_undone();
    if (pos == miss_val)
      break;
    status[pos] = 1;
    Face f = vf[pos];
    for (auto &eGen : LGenBig) {
      OnFace_inplace(eFaceImg, f, eGen);
      f_insert(eFaceImg);
    }
  }
  return OrbitSplittingSet(vf, SmaGRP);
}

template <typename Tgroup>
vectface OrbitSplittingListOrbit_spec(Tgroup const &BigGRP,
                                      Tgroup const &SmaGRP,
                                      const vectface &eListBig,
                                      std::string const &method_split,
                                      std::ostream &os) {
#ifdef TIMINGS
  MicrosecondTime time;
#endif
  os << "|BigGRP|=" << BigGRP.size() << " |SmaGRP|=" << SmaGRP.size() << "\n";
  using Tidx_value = uint16_t;
  WeightMatrix<true, int, Tidx_value> WMat;
  if (method_split == "repr") {
    WMat = WeightMatrixFromPairOrbits<Tgroup, Tidx_value>(SmaGRP);
  }
  vectface eListSma(BigGRP.n_act());
  for (auto &eSet : eListBig) {
    if (method_split == "repr") {
      vectface ListListSet =
          DoubleCosetDescription_Representation<Tgroup, Tidx_value>(
              BigGRP, SmaGRP, WMat, eSet);
      eListSma.append(ListListSet);
    }
    if (method_split == "canonic") {
      vectface ListListSet =
          DoubleCosetDescription_Canonic(BigGRP, SmaGRP, eSet, os);
      eListSma.append(ListListSet);
    }
    if (method_split == "exhaustive") {
      vectface ListListSet =
          DoubleCosetDescription_Exhaustive(BigGRP, SmaGRP, eSet);
      eListSma.append(ListListSet);
    }
  }
#ifdef TIMINGS
  os << "OrbitSplitting elapsed_microseconds=" << time
     << " |eListBig|=" << eListBig.size() << " |eListSma|=" << eListSma.size()
     << "\n";
#endif
  return eListSma;
}

template <typename Tgroup>
void PrintDoubleCosetCasesTestProblem(Tgroup const &BigGRP, Tgroup const &SmaGRP,
                                      const vectface &eListBig) {
  std::string strSizeSma = std::to_string(SmaGRP.size());
  std::string strSizeBig = std::to_string(BigGRP.size());
  std::string strLen = std::to_string(eListBig.size());
  std::string FileOut = "DoubleCoset_" + strSizeBig + "_" + strSizeSma + "_" + strLen;
  std::ofstream os(FileOut);
  WriteGroup(os, BigGRP);
  WriteGroup(os, SmaGRP);
  WriteListFace(os, eListBig);
}


template <typename Tgroup>
vectface OrbitSplittingListOrbit(Tgroup const &BigGRP, Tgroup const &SmaGRP,
                                 const vectface &eListBig, std::ostream &os) {
  std::string method_split = "canonic";
#ifdef PRINT_DOUBLE_COSETS_TEST_PROBLEM
  PrintDoubleCosetCasesTestProblem(BigGRP, SmaGRP, eListBig);
#endif
  return OrbitSplittingListOrbit_spec(BigGRP, SmaGRP, eListBig, method_split,
                                      os);
}

// This is for the lifting of orbits coming with the dual description of perfect
// cones.
template <typename Tgroup>
void OrbitSplittingPerfectFacet(Tgroup const &BigGRP, Tgroup const &SmaGRP,
                                const vectface &eListBig, std::ostream &os2,
                                std::ostream &os3) {
  std::cerr << "|BigGRP|=" << BigGRP.size() << " |SmaGRP|=" << SmaGRP.size()
            << "\n";
  size_t nb_orbit_big = eListBig.size();
  mpz_class nb_orbit_sma;
  size_t pos = 0;
  for (auto &eSet : eListBig) {
    pos++;
    vectface ListListSet =
        DoubleCosetDescription_Canonic(BigGRP, SmaGRP, eSet, std::cerr);
    mpz_class orb_siz = ListListSet.size();
    nb_orbit_sma += orb_siz;
    for (auto &eFace : ListListSet) {
      mpz_class res = getsetasint<mpz_class>(eFace);
      os3 << res << "\n";
    }
    std::cerr << "iInc=" << pos << " / " << nb_orbit_big
              << " |ListInc2|=" << nb_orbit_sma << " |LInc|=" << orb_siz
              << "\n";
  }
  os2 << nb_orbit_sma << "\n";
}

// clang-format off
#endif  // SRC_POLY_GRP_DOUBLECOSET_H_
// clang-format on
