// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_GRP_DOUBLECOSET_H_
#define SRC_POLY_GRP_DOUBLECOSET_H_

// clang-format off
#include "GRP_GroupFct.h"
#include "WeightMatrix.h"
#include "hopscotch_set.h"
#include "robin_set.h"
#include "sparse_set.h"
#include <limits>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef TIMINGS
# define TIMINGS_GROUP
#endif

#ifdef DEBUG
# define DEBUG_DOUBLE_COSET
#endif


template <typename Tgroup> struct FaceOrbitsizeGrpContainer {
private:
  Tgroup GRP;
  vectface vf;

public:
  FaceOrbitsizeGrpContainer(Tgroup _GRP, vectface &&_vf)
      : GRP(_GRP), vf(std::move(_vf)) {}
  std::pair<Face, typename Tgroup::Tint> GetPair(size_t const &idx_orb) const {
    using Tint = typename Tgroup::Tint;
    Face f = vf[idx_orb];
    Tint orbSize = GRP.OrbitSize_OnSets(f);
    return {f, orbSize};
  }
  size_t size() const { return vf.size(); }
};

//
// Double Coset Computation
//

template <typename Tgroup, typename Tidx_value>
vectface DoubleCosetDescription_Representation(
    std::vector<typename Tgroup::Telt> const &BigGens, Tgroup const &SmaGRP,
    WeightMatrix<true, int, Tidx_value> const &WMat, Face const &eList,
    typename Tgroup::Tint const &TotalSize, [[maybe_unused]] std::ostream &os) {
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
    SizeGen += SmaGRP.OrbitSize_OnSets(testList);
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
vectface DoubleCosetDescription_Canonic(
    std::vector<typename Tgroup::Telt> const &BigGens, Tgroup const &SmaGRP,
    Face const &eList, typename Tgroup::Tint const &TotalSize,
    std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  //
  Tint SizeGen = 0;
  std::unordered_set<Face> SetFace;
  Tidx n = eList.size();
  vectface CurrList(n);
  auto DoubleCosetInsertEntry_first = [&](Face const &testList) -> void {
    std::pair<Face, Tint> pairCan = SmaGRP.OptCanonicalImageOrbitSize(testList);
    if (SetFace.count(pairCan.first) > 0)
      return;
    CurrList.push_back(pairCan.first);
    SetFace.insert(pairCan.first);
    SizeGen += pairCan.second;
  };
  auto DoubleCosetInsertEntry_second = [&](Face const &testList) -> void {
    std::pair<Face, Tint> pairCan = SmaGRP.OptCanonicalImageOrbitSize(testList);
    if (SetFace.count(pairCan.first) > 0)
      return;
    SetFace.insert(pairCan.first);
    SizeGen += pairCan.second;
  };
  auto get_list_list_set = [&]() -> vectface {
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
  //  os << "DoubleCosetDescription_Canonic SizeGen=" << SizeGen << " TotalSize=" << TotalSize << "\n";
  //  os << "DoubleCosetDescription_Canonic |SetFace|=" << SetFace.size() << "\n";
  vectface ListListSet = get_list_list_set();
  //  os << "DoubleCosetDescription_Canonic |ListListSet|=" << ListListSet.size()
  //     << " n=" << ListListSet.get_n() << "\n";
  std::unordered_set<Face> PartialOrbit = SetFace;
  while (true) {
    //    os << "|ListListSet|=" << ListListSet.size() << "\n";
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
vectface DoubleCosetDescription_CanonicInitialTriv(
    std::vector<typename Tgroup::Telt> const &BigGens, Tgroup const &SmaGRP,
    Face const &eList, typename Tgroup::Tint const &TotalSize,
    std::ostream &os) {
  using Tidx = typename Tgroup::Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  //
  Tint SizeGen = 0;
  std::unordered_set<Face> SetFace;
  Tidx n = eList.size();
  vectface CurrList(n);
  auto DoubleCosetInsertEntry_first = [&](Face const &f) -> void {
    Face f_can = SmaGRP.CanonicalImageInitialTriv(f);
    if (SetFace.count(f_can) > 0)
      return;
    CurrList.push_back(f_can);
    SetFace.insert(f_can);
    SizeGen += SmaGRP.OrbitSize_OnSets(f);
  };
  auto DoubleCosetInsertEntry_second = [&](Face const &f) -> void {
    Face f_can = SmaGRP.CanonicalImageInitialTriv(f);
    if (SetFace.count(f_can) > 0)
      return;
    SetFace.insert(f_can);
    SizeGen += SmaGRP.OrbitSize_OnSets(f);
  };
  auto get_list_list_set = [&]() -> vectface {
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

template <typename Tgroup, typename T_hash_set>
vectface DoubleCosetDescription_Exhaustive_T(
    std::vector<typename Tgroup::Telt> const &BigGens, Tgroup const &SmaGRP,
    Face const &eList, typename Tgroup::Tint const &TotalSize,
    [[maybe_unused]] std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  using Tidx = typename Telt::Tidx;
  Tidx n = eList.size();
  //
  vectface vf(n);
  T_hash_set SetFace;
  size_t total_len = 0;
  auto f_insert = [&](const Face &f) -> void {
    if (SetFace.find(f) == SetFace.end()) {
      vf.push_back(f);
      SetFace.insert(f);
      total_len++;
    }
  };
  f_insert(eList);
  Face eFaceImg(n);
  size_t pos = 0;
  //  std::cerr << "TotalSize=" << TotalSize << "\n";
  while (true) {
    Tint SizeGen = total_len;
    //    std::cerr << "pos=" << pos << " total_len=" << total_len << "\n";
    if (SizeGen == TotalSize)
      break;
    size_t new_pos = total_len;
    for (size_t idx = pos; idx < total_len; idx++) {
      Face f = vf[idx];
      for (auto &eGen : BigGens) {
        OnFace_inplace(eFaceImg, f, eGen);
        f_insert(eFaceImg);
      }
    }
    pos = new_pos;
  }
  return OrbitSplittingSet_T(SetFace, SmaGRP);
}

template <typename Tgroup>
vectface DoubleCosetDescription_SingleCoset(
    Tgroup const &SmaGRP, Face const &eList,
    [[maybe_unused]] typename Tgroup::Tint const &TotalSize,
    std::vector<typename Tgroup::Telt> const &ListCos,
    [[maybe_unused]] std::ostream &os) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = eList.size();
  //
  vectface vf(n);
  Face eFaceImg(n);
  std::unordered_set<Face> SetFace;
  auto f_insert = [&](Face NewF) -> void {
    Face f_can = SmaGRP.OptCanonicalImage(NewF);
    if (SetFace.count(f_can) == 0) {
      vf.push_back(f_can);
      SetFace.insert(f_can);
    }
  };
  for (auto &eCos : ListCos) {
    OnFace_inplace(eFaceImg, eList, eCos);
    f_insert(eFaceImg);
  }
#ifdef DEBUG_DOUBLE_COSET
  using Tint = typename Tgroup::Tint;
  Tint GenSize = 0;
  for (auto & f : SetFace) {
    GenSize += SmaGRP.OrbitSize_OnSets(f);
  }
  if (GenSize != TotalSize) {
    std::cerr << "GenSize=" << GenSize << " TotalSize=" << TotalSize << "\n";
    throw TerminalException{1};
  }
#endif
  return vf;
}

template <typename Tgroup, typename Tface_orbitsize>
void PrintDoubleCosetCasesTestProblem(
    Tgroup const &BigGRP, Tgroup const &SmaGRP,
    const Tface_orbitsize &ListFaceOrbitsize) {
  using Tint = typename Tgroup::Tint;
  size_t nbOrbit = ListFaceOrbitsize.size();
  std::string strSizeSma = std::to_string(SmaGRP.size());
  std::string strSizeBig = std::to_string(BigGRP.size());
  std::string strLen = std::to_string(nbOrbit);
  std::string strN = std::to_string(int(BigGRP.n_act()));
  std::string Prefix = "DoubleCoset_n" + strN + "_big" + strSizeBig + "_sma" +
                       strSizeSma + "_vf" + strLen + "_idx";
  std::string FileOut = FindAvailableFileFromPrefix(Prefix);
  std::ofstream os(FileOut);
  WriteGroup(os, BigGRP);
  WriteGroup(os, SmaGRP);
  os << nbOrbit << "\n";
  for (size_t i_orbit = 0; i_orbit < nbOrbit; i_orbit++) {
    std::pair<Face, Tint> pair = ListFaceOrbitsize.GetPair(i_orbit);
    Face const &f = pair.first;
    size_t len = f.size();
    os << len;
    for (size_t i = 0; i < len; i++) {
      int val = f[i];
      os << " " << val;
    }
    os << "\n";
  }
}

template <typename Tgroup, typename Tface_orbitsize, typename Fterminal>
vectface DoubleCosetDescription_Representation_Block(
    Tgroup const &BigGRP, Tgroup const &SmaGRP,
    const Tface_orbitsize &ListFaceOrbitsize, Fterminal f_terminal,
    std::ostream &os) {
  // Use generators to build more and more elements and check for equivalence.
  // In some rare cases it does not work and another scheme has to be used.
  using Tidx_value = uint16_t;
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  WeightMatrix<true, int, Tidx_value> WMat =
      WeightMatrixFromPairOrbits<Tgroup, Tidx_value>(SmaGRP);
  Tidx n = BigGRP.n_act();
  vectface eListSma(n);
  std::vector<Telt> BigGens = BigGRP.SmallGeneratingSet();
  size_t nbOrbit = ListFaceOrbitsize.size();
  for (size_t i_orbit = 0; i_orbit < nbOrbit; i_orbit++) {
    std::pair<Face, Tint> pair = ListFaceOrbitsize.GetPair(i_orbit);
    Face const &eSet = pair.first;
    Tint const &TotalSize = pair.second;
    if (f_terminal())
      break;
    vectface ListListSet =
        DoubleCosetDescription_Representation<Tgroup, Tidx_value>(
            BigGens, SmaGRP, WMat, eSet, TotalSize, os);
    eListSma.append(ListListSet);
  }
  return eListSma;
}

template <typename Tgroup, typename Tface_orbitsize, typename Fterminal>
vectface
DoubleCosetDescription_Canonic_Block(Tgroup const &BigGRP, Tgroup const &SmaGRP,
                                     const Tface_orbitsize &ListFaceOrbitsize,
                                     Fterminal f_terminal, std::ostream &os) {
  // A reliable technique, it has the same issues as the representative method
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  Tidx n = BigGRP.n_act();
  vectface eListSma(n);
  std::vector<Telt> BigGens = BigGRP.SmallGeneratingSet();
  size_t nbOrbit = ListFaceOrbitsize.size();
  for (size_t i_orbit = 0; i_orbit < nbOrbit; i_orbit++) {
    std::pair<Face, Tint> pair = ListFaceOrbitsize.GetPair(i_orbit);
    Face const &eSet = pair.first;
    Tint const &TotalSize = pair.second;
    //    Tint orbSize = BigGRP.OrbitSize_OnSets(eSet);
    //    os << "osbSize=" << orbSize << " TotalSize=" << TotalSize << "\n";
    //    if (orbSize != TotalSize) {
    //      std::cerr << "orbSize=" << orbSize << " TotalSize=" << TotalSize << "\n";
    //      std::cerr << "Inconsistency in the Size computation\n";
    //      throw TerminalException{1};
    //    }
    if (f_terminal())
      break;
    vectface ListListSet = DoubleCosetDescription_Canonic<Tgroup>(
        BigGens, SmaGRP, eSet, TotalSize, os);
    eListSma.append(ListListSet);
  }
  return eListSma;
}

template <typename Tgroup, typename Tface_orbitsize, typename Fterminal>
vectface DoubleCosetDescription_CanonicInitialTriv_Block(
    Tgroup const &BigGRP, Tgroup const &SmaGRP,
    const Tface_orbitsize &ListFaceOrbitsize, Fterminal f_terminal,
    std::ostream &os) {
  // A reliable technique, it has the same issues as the representative method
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  Tidx n = BigGRP.n_act();
  vectface eListSma(n);
  std::vector<Telt> BigGens = BigGRP.SmallGeneratingSet();
  size_t nbOrbit = ListFaceOrbitsize.size();
  for (size_t i_orbit = 0; i_orbit < nbOrbit; i_orbit++) {
    std::pair<Face, Tint> pair = ListFaceOrbitsize.GetPair(i_orbit);
    Face const &eSet = pair.first;
    Tint const &TotalSize = pair.second;
    if (f_terminal())
      break;
    vectface ListListSet = DoubleCosetDescription_CanonicInitialTriv<Tgroup>(
        BigGens, SmaGRP, eSet, TotalSize, os);
    eListSma.append(ListListSet);
  }
  return eListSma;
}

template <typename Tgroup, typename T_hash_set, typename Tface_orbitsize,
          typename Fterminal>
vectface DoubleCosetDescription_Exhaustive_Block_T(
    Tgroup const &BigGRP, Tgroup const &SmaGRP,
    const Tface_orbitsize &ListFaceOrbitsize, Fterminal f_terminal,
    std::ostream &os) {
  // Generate the full orbit for the big group and then split it using the small
  // group. To be preferred when SmaGRP is really small.
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  Tidx n = BigGRP.n_act();
  vectface eListSma(n);
  std::vector<Telt> BigGens = BigGRP.SmallGeneratingSet();
  size_t nbOrbit = ListFaceOrbitsize.size();
  for (size_t i_orbit = 0; i_orbit < nbOrbit; i_orbit++) {
    std::pair<Face, Tint> pair = ListFaceOrbitsize.GetPair(i_orbit);
    Face const &eSet = pair.first;
    Tint const &TotalSize = pair.second;
    if (f_terminal())
      break;
    vectface ListListSet =
        DoubleCosetDescription_Exhaustive_T<Tgroup, T_hash_set>(
            BigGens, SmaGRP, eSet, TotalSize, os);
    eListSma.append(ListListSet);
  }
  return eListSma;
}

template <typename Tgroup, typename Tface_orbitsize, typename Fterminal>
vectface DoubleCosetDescription_CanonicInitialTriv_ExhaustiveLimit_Block(Tgroup const &BigGRP,
                                                                         Tgroup const &SmaGRP,
                                                                         const Tface_orbitsize &ListFaceOrbitsize,
                                                                         typename Tgroup::Tint const& n_limit,
                                                                         Fterminal f_terminal,
                                                                         std::ostream &os) {
  // Generate the full orbit for the big group and then split it using the small group.
  // To be preferred when SmaGRP is really small.
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  using Tint = typename Tgroup::Tint;
  using T_hash_set = std::unordered_set<Face>;
  Tidx n = BigGRP.n_act();
  vectface eListSma(n);
  std::vector<Telt> BigGens = BigGRP.SmallGeneratingSet();
  size_t nbOrbit = ListFaceOrbitsize.size();
  for (size_t i_orbit=0; i_orbit<nbOrbit; i_orbit++) {
    std::pair<Face,Tint> pair = ListFaceOrbitsize.GetPair(i_orbit);
    Face const& eSet = pair.first;
    Tint const& TotalSize = pair.second;
    if (f_terminal()) break;
    if (n_limit <= TotalSize) {
      vectface ListListSet =
        DoubleCosetDescription_Exhaustive_T<Tgroup,T_hash_set>(BigGens, SmaGRP, eSet, TotalSize, os);
      eListSma.append(ListListSet);
    } else {
      vectface ListListSet =
        DoubleCosetDescription_CanonicInitialTriv<Tgroup>(BigGens, SmaGRP, eSet, TotalSize, os);
      eListSma.append(ListListSet);
    }
  }
  return eListSma;
}

template <typename Tgroup, typename Tface_orbitsize, typename Fterminal>
vectface DoubleCosetDescription_SingleCoset_Block(
    Tgroup const &BigGRP, Tgroup const &SmaGRP,
    const Tface_orbitsize &ListFaceOrbitsize, Fterminal f_terminal,
    std::ostream &os) {
  // Compute the left cosets, that is BigGRP = cup_i g_i SmaGRP
  // Then form the representatives x g_i and test for equivalence.
  // Should be great when the index is small.
  using Telt = typename Tgroup::Telt;
  using Tint = typename Tgroup::Tint;
  using Tidx = typename Telt::Tidx;
  Tidx n = BigGRP.n_act();
  vectface eListSma(n);
  std::vector<Telt> ListCos = BigGRP.LeftTransversal_Direct(SmaGRP);
  size_t nbOrbit = ListFaceOrbitsize.size();
  for (size_t i_orbit = 0; i_orbit < nbOrbit; i_orbit++) {
    std::pair<Face, Tint> pair = ListFaceOrbitsize.GetPair(i_orbit);
    Face const &eSet = pair.first;
    if (f_terminal())
      break;
    Tint const &TotalSize = pair.second;
    vectface ListListSet = DoubleCosetDescription_SingleCoset(
        SmaGRP, eSet, TotalSize, ListCos, os);
    eListSma.append(ListListSet);
  }
  return eListSma;
}

template <typename Tgroup, typename Tface_orbitsize, typename Fterminal>
vectface
OrbitSplittingListOrbitKernel_spec(Tgroup const &BigGRP, Tgroup const &SmaGRP,
                                   const Tface_orbitsize &ListFaceOrbitsize,
                                   std::string const &method_split,
                                   Fterminal f_terminal, std::ostream &os) {
  using Tint = typename Tgroup::Tint;
#ifdef TIMINGS_GROUP
  MicrosecondTime time;
#endif
  os << "|BigGRP|=" << BigGRP.size() << " |SmaGRP|=" << SmaGRP.size()
     << " |vf|=" << ListFaceOrbitsize.size() << " method_split=" << method_split
     << std::endl;
#ifdef PRINT_DOUBLE_COSETS_TEST_PROBLEM
  PrintDoubleCosetCasesTestProblem(BigGRP, SmaGRP, ListFaceOrbitsize);
#endif
  auto get_split = [&]() -> vectface {
    if (method_split == "repr") {
      return DoubleCosetDescription_Representation_Block(
          BigGRP, SmaGRP, ListFaceOrbitsize, f_terminal, os);
    }
    if (method_split == "canonic") {
      return DoubleCosetDescription_Canonic_Block(
          BigGRP, SmaGRP, ListFaceOrbitsize, f_terminal, os);
    }
    if (method_split == "canonic_initial_triv") {
      return DoubleCosetDescription_CanonicInitialTriv_Block(
          BigGRP, SmaGRP, ListFaceOrbitsize, f_terminal, os);
    }
    if (method_split == "exhaustive_std") {
      return DoubleCosetDescription_Exhaustive_Block_T<
          Tgroup, std::unordered_set<Face>>(BigGRP, SmaGRP, ListFaceOrbitsize,
                                            f_terminal, os);
    }
    if (method_split == "exhaustive_sparse") {
      return DoubleCosetDescription_Exhaustive_Block_T<Tgroup,
                                                       tsl::sparse_set<Face>>(
          BigGRP, SmaGRP, ListFaceOrbitsize, f_terminal, os);
    }
    if (method_split == "exhaustive_robin" || method_split == "exhaustive") {
      // That variant seems a little bit faster
      return DoubleCosetDescription_Exhaustive_Block_T<Tgroup,
                                                       tsl::robin_set<Face>>(
          BigGRP, SmaGRP, ListFaceOrbitsize, f_terminal, os);
    }
    if (method_split == "exhaustive_hopscotch") {
      return DoubleCosetDescription_Exhaustive_Block_T<
          Tgroup, tsl::hopscotch_set<Face>>(BigGRP, SmaGRP, ListFaceOrbitsize,
                                            f_terminal, os);
    }
    if (method_split == "single_cosets") {
      return DoubleCosetDescription_SingleCoset_Block(
          BigGRP, SmaGRP, ListFaceOrbitsize, f_terminal, os);
    }
    std::optional<std::string> opt =
      get_postfix(method_split, "canonic_initial_triv_exhaustive_limit");
    if (opt) {
      // In this scheme, we use the exhaustive method when the orbit size is lower than n.
      // Otherwise, we use the canonic_initial_triv method.
      std::string const& n_str = *opt;
      Tint n = ParseScalar<Tint>(n_str);
      return DoubleCosetDescription_CanonicInitialTriv_ExhaustiveLimit_Block(BigGRP, SmaGRP, ListFaceOrbitsize, n, f_terminal, os);
    }
    std::cerr << "Failed to find a matching entry\n";
    throw TerminalException{1};
  };
  vectface eListSma = get_split();
#ifdef TIMINGS_GROUP
  os << "OrbitSplitting elapsed_microseconds=" << time
     << " |eListBig|=" << ListFaceOrbitsize.size()
     << " |eListSma|=" << eListSma.size() << std::endl;
#endif
  return eListSma;
}

template <typename Tgroup, typename Tface_orbitsize>
vectface OrbitSplittingListOrbit_spec(Tgroup const &BigGRP,
                                      Tgroup const &SmaGRP,
                                      const Tface_orbitsize &ListFaceOrbitsize,
                                      std::string const &method_split,
                                      std::ostream &os) {
  using Tint = typename Tgroup::Tint;
  auto f_direct = [&](std::string const &the_method) -> vectface {
    auto f_terminal = [&]() -> bool { return false; };
    return OrbitSplittingListOrbitKernel_spec(BigGRP, SmaGRP, ListFaceOrbitsize,
                                              the_method, f_terminal, os);
  };
  if (method_split != "guess") {
    return f_direct(method_split);
  } else {
    if (ListFaceOrbitsize.size() < 3000) {
      // Too small orbit size, sampling is too expensive.
      Tint index = BigGRP.size() / SmaGRP.size();
      if (index < 50) {
        return f_direct("single_cosets");
      }
      if (SmaGRP.size() < 200) {
        return f_direct("exhaustive");
      }
      return f_direct("canonic");
    }
    std::vector<std::string> Lmethod = {"canonic", "canonic_initial_triv",
                                        "exhaustive_std", "exhaustive_robin",
                                        "single_cosets"};
    int64_t max_val = std::numeric_limits<int64_t>::max();
    int64_t smallest_time = max_val;
    std::string chosen_method = "unset";
    auto evaluate_time = [&](std::string const &the_method) -> int64_t {
      int64_t the_time = max_val;
      int iter = 0;
      int n_iter = 100; // Enough for sampling the speed.
      NanosecondTime time;
      auto f_terminal = [&]() -> bool {
        iter++;
        int64_t duration = time.const_eval_int64();
        if (iter == n_iter) {
          the_time = duration;
          os << "Finished the n_iter=" << n_iter << "\n";
          return true;
        }
        if (duration > smallest_time) {
          os << "duration already larger, no need to continue for method="
             << the_method << "\n";
          return true;
        }
        return false;
      };
      (void)OrbitSplittingListOrbitKernel_spec(
          BigGRP, SmaGRP, ListFaceOrbitsize, the_method, f_terminal, os);
      return the_time;
    };
    for (auto &method : Lmethod) {
      int64_t time = evaluate_time(method);
      if (time < smallest_time) {
        chosen_method = method;
        smallest_time = time;
      }
    }
    return f_direct(chosen_method);
  }
}

template <typename Tgroup, typename Tface_orbitsize>
vectface OrbitSplittingListOrbit(Tgroup const &BigGRP, Tgroup const &SmaGRP,
                                 const Tface_orbitsize &eListBig,
                                 std::ostream &os) {
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
  using Tint = typename Tgroup::Tint;
  using Telt = typename Tgroup::Telt;
  std::cerr << "|BigGRP|=" << BigGRP.size() << " |SmaGRP|=" << SmaGRP.size()
            << "\n";
  size_t nb_orbit_big = eListBig.size();
  Tint nb_orbit_sma;
  size_t pos = 0;
  std::vector<Telt> BigGens = BigGRP.SmallGeneratingSet();
  for (auto &eSet : eListBig) {
    pos++;
    Tint TotalSize = BigGRP.OrbitSize_OnSets(eSet);
    vectface ListListSet = DoubleCosetDescription_Canonic<Tgroup>(
        BigGens, SmaGRP, eSet, TotalSize, std::cerr);
    Tint orb_siz = ListListSet.size();
    nb_orbit_sma += orb_siz;
    for (auto &eFace : ListListSet) {
      Tint res = getsetasint<Tint>(eFace);
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
