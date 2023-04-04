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
    std::vector<typename Tgroup::Telt> const &BigGens, Tgroup const &SmaGRP,
    WeightMatrix<true, int, Tidx_value> const &WMat,
    Face const &eList, typename Tgroup::Tint const& TotalSize,
    [[maybe_unused]] std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
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
        for (auto const &eGen : BigGens) {
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
      for (auto &eGen : BigGens) {
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
vectface DoubleCosetDescription_Canonic(std::vector<typename Tgroup::Telt> const &BigGens, Tgroup const &SmaGRP,
                                        Face const &eList, typename Tgroup::Tint const& TotalSize,
                                        std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  //
  Tint SizeGen = 0;
  std::unordered_set<Face> SetFace;
  Tidx n = eList.size();
  vectface CurrList(n);
  auto DoubleCosetInsertEntry_first = [&](Face const &testList) -> void {
    std::pair<Face,Tint> pairCan = SmaGRP.OptCanonicalImageOrbitSize(testList);
    if (SetFace.count(pairCan.first) > 0)
      return;
    CurrList.push_back(pairCan.first);
    SetFace.insert(pairCan.first);
    SizeGen += pairCan.second;
  };
  auto DoubleCosetInsertEntry_second = [&](Face const &testList) -> void {
    std::pair<Face,Tint> pairCan = SmaGRP.OptCanonicalImageOrbitSize(testList);
    if (SetFace.count(pairCan.first) > 0)
      return;
    SetFace.insert(pairCan.first);
    SizeGen += pairCan.second;
  };
  auto get_list_list_set=[&]() -> vectface {
    vectface ListListSet(n);
    for (auto &eFace : SetFace)
      ListListSet.push_back(eFace);
    return ListListSet;
  };
  DoubleCosetInsertEntry_first(eList);
  Face eFaceImg(n);
  while (true) {
    if (SizeGen == TotalSize)
      return get_list_list_set();
    if (CurrList.size() == 0)
      break;
    Face eFace = CurrList.pop();
    for (auto const &eGen : BigGens) {
      OnFace_inplace(eFaceImg, eFace, eGen);
      DoubleCosetInsertEntry_first(eFaceImg);
    }
  }
  vectface ListListSet = get_list_list_set();
  //  os << "After Iteration loop SizeGen=" << SizeGen << " TotalSize=" <<
  //  TotalSize << "\n";
  std::unordered_set<Face> PartialOrbit = SetFace;
  while (true) {
    Face eFace = ListListSet.pop();
    for (auto &eGen : BigGens) {
      OnFace_inplace(eFaceImg, eFace, eGen);
      if (PartialOrbit.count(eFaceImg) == 0) {
        PartialOrbit.insert(eFaceImg);
        ListListSet.push_back(eFaceImg);
        DoubleCosetInsertEntry_second(eFaceImg);
        if (SizeGen == TotalSize) {
          return get_list_list_set();
        }
      }
    }
  }
  os << "Likely not reachable stage\n";
  throw TerminalException{1};
}

template <typename Tgroup>
vectface DoubleCosetDescription_Exhaustive(std::vector<typename Tgroup::Telt> const &BigGens, Tgroup const &SmaGRP,
                                           Face const &eList, typename Tgroup::Tint const& TotalSize,
                                           [[maybe_unused]] std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  using Tidx = typename Telt::Tidx;
  Tidx n = eList.size();
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
    Tint SizeGen = vf.size();
    if (SizeGen == TotalSize)
      break;
    size_t pos = get_undone();
    if (pos == miss_val)
      break;
    status[pos] = 1;
    Face f = vf[pos];
    for (auto &eGen : BigGens) {
      OnFace_inplace(eFaceImg, f, eGen);
      f_insert(eFaceImg);
    }
  }
  return OrbitSplittingSet(vf, SmaGRP);
}

template <typename Tgroup>
#ifdef CHECK_DOUBLE_COSET
vectface DoubleCosetDescription_SingleCoset(Tgroup const &SmaGRP,
                                            Face const &eList, typename Tgroup::Tint const& TotalSize,
                                            std::vector<typename Tgroup::Telt> const& ListCos,
                                            [[maybe_unused]] std::ostream &os) {
#else
vectface DoubleCosetDescription_SingleCoset(Tgroup const &SmaGRP,
                                            Face const &eList, std::vector<typename Tgroup::Telt> const& ListCos,
                                            [[maybe_unused]] std::ostream &os) {
#endif
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
#ifdef CHECK_DOUBLE_COSET
  using Tint = typename Tgroup::Tint;
  Tint GenSize = 0;
#endif
  Tidx n = eList.size();
  //
  vectface vf(n);
  Face eFaceImg(n);
  std::unordered_set<Face> SetFace;
  auto f_insert=[&](Face NewF) -> void {
    Face f_can = SmaGRP.OptCanonicalImage(NewF);
    if (SetFace.count(f_can) == 0) {
      vf.push_back(f_can);
      SetFace.insert(f_can);
#ifdef CHECK_DOUBLE_COSET
      Tgroup TheStab = SmaGRP.Stabilizer_OnSets(NewF);
      Tint OrbitSize = SmaGRP.size() / TheStab.size();
      GenSize += OrbitSize;
#endif
    }
  };
  for (auto & eCos : ListCos) {
    OnFace_inplace(eFaceImg, eList, eCos);
    f_insert(eFaceImg);
  }
#ifdef CHECK_DOUBLE_COSET
  if (GenSize != TotalSize) {
    std::cerr << "GenSize=" << GenSize << " TotalSize=" << TotalSize << "\n";
    throw TerminalException{1};
  }
#endif
  return vf;
}


template <typename Tgroup>
void PrintDoubleCosetCasesTestProblem(Tgroup const &BigGRP, Tgroup const &SmaGRP,
                                      const vectface &eListBig) {
  std::string strSizeSma = std::to_string(SmaGRP.size());
  std::string strSizeBig = std::to_string(BigGRP.size());
  std::string strLen = std::to_string(eListBig.size());
  std::string strN = std::to_string(int(BigGRP.n_act()));
  std::string Prefix = "DoubleCoset_n" + strN + "_big" + strSizeBig + "_sma" + strSizeSma + "_vf" + strLen + "_idx";
  std::string FileOut = FindAvailableFileFromPrefix(Prefix);
  std::ofstream os(FileOut);
  WriteGroup(os, BigGRP);
  WriteGroup(os, SmaGRP);
  WriteListFace(os, eListBig);
}

template <typename Tgroup>
vectface DoubleCosetDescription_Representation_Block(Tgroup const &BigGRP,
                                      Tgroup const &SmaGRP,
                                      const vectface &eListBig,
                                      std::ostream &os) {
  // Use generators to build more and more elements and check for equivalence.
  // In some rare cases it does not work and another scheme has to be used.
  using Tidx_value = uint16_t;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  WeightMatrix<true, int, Tidx_value> WMat = WeightMatrixFromPairOrbits<Tgroup, Tidx_value>(SmaGRP);
  Tidx n = BigGRP.n_act();
  vectface eListSma(n);
  std::vector<Telt> BigGens = BigGRP.SmallGeneratingSet();
  for (auto &eSet : eListBig) {
    Tgroup TheStab = BigGRP.Stabilizer_OnSets(eSet);
    Tint TotalSize = BigGRP.size() / TheStab.size();
    vectface ListListSet =
      DoubleCosetDescription_Representation<Tgroup, Tidx_value>(BigGens, SmaGRP, WMat, eSet, TotalSize, os);
    eListSma.append(ListListSet);
  }
  return eListSma;
}

template <typename Tgroup>
vectface DoubleCosetDescription_Canonic_Block(Tgroup const &BigGRP,
                                      Tgroup const &SmaGRP,
                                      const vectface &eListBig,
                                      std::ostream &os) {
  // A reliable technique, it has the same issues as the representative method
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  Tidx n = BigGRP.n_act();
  vectface eListSma(n);
  std::vector<Telt> BigGens = BigGRP.SmallGeneratingSet();
  for (auto &eSet : eListBig) {
    Tgroup TheStab = BigGRP.Stabilizer_OnSets(eSet);
    Tint TotalSize = BigGRP.size() / TheStab.size();
    vectface ListListSet =
      DoubleCosetDescription_Canonic<Tgroup>(BigGens, SmaGRP, eSet, TotalSize, os);
    eListSma.append(ListListSet);
  }
  return eListSma;
}

template <typename Tgroup>
vectface DoubleCosetDescription_Exhaustive_Block(Tgroup const &BigGRP,
                                      Tgroup const &SmaGRP,
                                      const vectface &eListBig,
                                      std::ostream &os) {
  // Generate the full orbit for the big group and then split it using the small group.
  // To be preferred when SmaGRP is really small.
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  Tidx n = BigGRP.n_act();
  vectface eListSma(n);
  std::vector<Telt> BigGens = BigGRP.SmallGeneratingSet();
  for (auto &eSet : eListBig) {
    Tgroup TheStab = BigGRP.Stabilizer_OnSets(eSet);
    Tint TotalSize = BigGRP.size() / TheStab.size();
    vectface ListListSet =
      DoubleCosetDescription_Exhaustive<Tgroup>(BigGens, SmaGRP, eSet, TotalSize, os);
    eListSma.append(ListListSet);
  }
  return eListSma;
}

template <typename Tgroup>
vectface DoubleCosetDescription_SingleCoset_Block(Tgroup const &BigGRP,
                                      Tgroup const &SmaGRP,
                                      const vectface &eListBig,
                                      std::ostream &os) {
  // Compute the left cosets, that is BigGRP = cup_i g_i SmaGRP
  // Then form the representatives x g_i and test for equivalence.
  // Should be great when the index is small.
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = BigGRP.n_act();
  vectface eListSma(n);
  std::vector<Telt> ListCos = BigGRP.LeftTransversal_Direct(SmaGRP);
  for (auto &eSet : eListBig) {
#ifdef CHECK_DOUBLE_COSET
    Tgroup TheStab = BigGRP.Stabilizer_OnSets(eSet);
    Tint TotalSize = BigGRP.size() / TheStab.size();
    vectface ListListSet =
      DoubleCosetDescription_SingleCoset(SmaGRP, eSet, TotalSize, ListCos, os);
#else
    vectface ListListSet =
      DoubleCosetDescription_SingleCoset(SmaGRP, eSet, ListCos, os);
#endif
    eListSma.append(ListListSet);
  }
  return eListSma;
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
  os << "|BigGRP|=" << BigGRP.size() << " |SmaGRP|=" << SmaGRP.size() << " |vf|=" << eListBig.size() << " method_split=" << method_split << std::endl;
#ifdef PRINT_DOUBLE_COSETS_TEST_PROBLEM
  PrintDoubleCosetCasesTestProblem(BigGRP, SmaGRP, eListBig);
#endif
  auto get_split=[&]() -> vectface {
    if (method_split == "repr") {
      return DoubleCosetDescription_Representation_Block(BigGRP, SmaGRP, eListBig, os);
    }
    if (method_split == "canonic") {
      return DoubleCosetDescription_Canonic_Block(BigGRP, SmaGRP, eListBig, os);
    }
    if (method_split == "exhaustive") {
      return DoubleCosetDescription_Exhaustive_Block(BigGRP, SmaGRP, eListBig, os);
    }
    if (method_split == "single_cosets") {
      return DoubleCosetDescription_SingleCoset_Block(BigGRP, SmaGRP, eListBig, os);
    }
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  };
  vectface eListSma = get_split();
#ifdef TIMINGS
    os << "OrbitSplitting elapsed_microseconds=" << time
       << " |eListBig|=" << eListBig.size() << " |eListSma|=" << eListSma.size()
       << std::endl;
#endif
  return eListSma;
}


template <typename Tgroup>
vectface OrbitSplittingListOrbit(Tgroup const &BigGRP, Tgroup const &SmaGRP,
                                 const vectface &eListBig, std::ostream &os) {
  std::string method_split = "canonic";
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
