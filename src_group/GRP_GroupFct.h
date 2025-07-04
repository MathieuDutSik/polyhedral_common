// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_GRP_GROUPFCT_H_
#define SRC_GROUP_GRP_GROUPFCT_H_

// clang-format off
#include "Basic_file.h"
#include "Basic_string.h"
#include "Boost_bitset.h"
#include "COMB_Stor.h"
#include "Temp_common.h"
#include "Timings.h"
#include "hash_functions.h"
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_GROUP
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_GROUP_FCT
#endif

//
// permutation functions
//

template <typename Telt> std::vector<int> PermutationOrbit(Telt const &ePerm) {
  int siz = ePerm.size();
  std::vector<int> StatusOrbit(siz, -1);
  int idxOrbit = 0;
  auto GetUnsetPoint = [&](void) -> int {
    for (int i = 0; i < siz; i++)
      if (StatusOrbit[i] == -1)
        return i;
    return -1;
  };
  while (true) {
    int iPoint = GetUnsetPoint();
    if (iPoint == -1)
      break;
    int iPointWork = iPoint;
    while (true) {
      StatusOrbit[iPointWork] = idxOrbit;
      iPointWork = OnPoints(iPointWork, ePerm);
      if (iPointWork == iPoint)
        break;
    }
    idxOrbit++;
  }
  return StatusOrbit;
}

template <typename Telt>
void WritePermutationGAP(std::ostream &os, Telt const &ePerm) {
  if (ePerm.isIdentity()) {
    os << "()";
    return;
  }
  std::vector<int> eVectOrbit = PermutationOrbit(ePerm);
  int nbOrbit = VectorMax(eVectOrbit) + 1;
  int siz = eVectOrbit.size();
  for (int iOrbit = 0; iOrbit < nbOrbit; iOrbit++) {
    int nbMatch = 0;
    int eFirst = -1;
    for (int i = 0; i < siz; i++)
      if (eVectOrbit[i] == iOrbit) {
        if (nbMatch == 0) {
          eFirst = i;
        }
        nbMatch++;
      }
    std::vector<int> TheList(nbMatch);
    TheList[0] = eFirst;
    for (int i = 1; i < nbMatch; i++)
      TheList[i] = OnPoints(TheList[i - 1], ePerm);
    os << "(";
    for (int i = 0; i < nbMatch; i++) {
      if (i > 0)
        os << ",";
      int val = TheList[i] + 1;
      os << val;
    }
    os << ")";
  }
}

//
// Template general code for input output of groups
//

template <typename Tgroup> Tgroup ReadGroup(std::istream &is) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  if (!is.good()) {
    std::cerr << "ReadGroup operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  int nbGen;
  int n_i;
  is >> n_i;
  is >> nbGen;
  Tidx n = Tidx(n_i);
  std::vector<Telt> ListGen;
  for (int iGen = 0; iGen < nbGen; iGen++) {
    std::vector<Tidx> v(n);
    for (Tidx i = 0; i < n; i++) {
      int eVal_i;
      is >> eVal_i;
      Tidx eVal = Tidx(eVal_i);
      if (eVal >= n) {
        std::cerr << "n=" << n_i << " nbGen=" << nbGen << "\n";
        std::cerr << "Error in ReadGroup function at i=" << i << "/" << n_i
                  << "\n";
        std::cerr << "Number of elements acted on n=" << n << " iGen=" << iGen
                  << "/" << nbGen << "\n";
        std::cerr << "But eVal=" << eVal << " eVal_i=" << eVal_i << "\n";
        throw TerminalException{1};
      }
      v[i] = eVal;
    }
    ListGen.emplace_back(std::move(Telt(std::move(v))));
  }
  return Tgroup(ListGen, n);
}

template <typename Tgroup> Tgroup ReadGroupFile(std::string const &file_name) {
  if (!IsExistingFile(file_name)) {
    std::cerr << "Error in ReadGroupFile\n";
    std::cerr << "file_name=" << file_name << " does not appear to exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(file_name);
  return ReadGroup<Tgroup>(is);
}

template <typename Tgroup>
void WriteGroup(std::ostream &os, Tgroup const &TheGRP) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  int nbGen = ListGen.size();
  Tidx n = TheGRP.n_act();
  os << int(n) << " " << nbGen << "\n";
  for (auto &eGen : ListGen) {
    for (Tidx i = 0; i < n; i++) {
      Tidx eVal = OnPoints(i, eGen);
      os << " " << int(eVal);
    }
    os << "\n";
  }
}

template<typename Tgroup>
void PrintRepresentativeAction_OnSets_GRP_f1_f2(Tgroup const& GRP, Face const& f1, Face const& f2) {
  std::string prefix = "RepresentativeAction_OnSets_GRP_f1_f2_idx";
  std::string FileOut = FindAvailableFileFromPrefix(prefix);
  std::ofstream os(FileOut);
  WriteGroup(os, GRP);
  auto f_print=[&](Face const& f) -> void {
    for (size_t i=0; i<f.size(); i++) {
      os << " " << f[i];
    }
    os << "\n";
  };
  f_print(f1);
  f_print(f2);
}


template <typename Tgroup> std::string StringGroup(Tgroup const &TheGRP) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  int nbGen = ListGen.size();
  Tidx n = TheGRP.n_act();
  int n_i = static_cast<int>(n);
  std::string strOut = std::to_string(n_i) + " " + std::to_string(nbGen);
  for (auto &eGen : ListGen) {
    for (Tidx i = 0; i < n; i++) {
      Tidx eVal = OnPoints(i, eGen);
      int eVal_i = static_cast<int>(eVal);
      strOut += " " + std::to_string(eVal_i);
    }
  }
  return strOut;
}

template <typename Tgroup>
void WriteGroupFile(std::string const &eFile, Tgroup const &TheGRP) {
  std::ofstream os(eFile);
  WriteGroup(os, TheGRP);
}

template <typename Tgroup>
void WriteGroupMakeUp(std::ostream &os, Tgroup const &TheGRP) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  os << "nb acting element=" << TheGRP.n_act() << "\n";
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  int nbGen = 0;
  Tidx n = TheGRP.n_act();
  for (auto &eGen : ListGen) {
    for (Tidx i = 0; i < TheGRP.n_act(); i++) {
      Tidx eVal = OnPoints(i, eGen);
      os << " " << int(eVal);
    }
    os << "\n";
    nbGen++;
  }
  os << "nbGen=" << nbGen << "\n";
  os << "size=" << TheGRP.size() << "\n";
}

template <typename Tgroup>
void WriteGroupGAP(std::ostream &os, Tgroup const &TheGRP) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  os << "local eListList, ListGen, GRP;\n";
  os << "eListList:=[\n";
  bool IsFirst = true;
  Tidx n = TheGRP.n_act();
  for (auto &eGen : ListGen) {
    if (!IsFirst)
      os << ",\n";
    IsFirst = false;
    os << "[";
    for (Tidx i = 0; i < n; i++) {
      int eVal = 1 + OnPoints(i, eGen);
      if (i > 0)
        os << ",";
      os << eVal;
    }
    os << "]";
  }
  os << "];\n";
  os << "ListGen:=List(eListList, PermList);\n";
  os << "GRP:=Group(ListGen);\n";
  os << "SetSize(GRP, " << TheGRP.size() << ");\n";
  os << "return GRP;\n";
}

template <typename Tgroup>
void WriteGroupFileGAP(std::string const& eFile, Tgroup const &TheGRP) {
  std::ofstream osf(eFile);
  WriteGroupGAP(osf, TheGRP);
}

template <typename Tgroup>
void WriteGroupFormat(std::string const &FileGroup,
                      std::string const &OutFormat, Tgroup const &TheGRP) {
  auto f_print = [&](std::ostream &os) -> void {
    if (OutFormat == "CPP") {
      return WriteGroup(os, TheGRP);
    }
    if (OutFormat == "GAP") {
      return WriteGroupGAP(os, TheGRP);
    }
    std::cerr
        << "Failed to find a matching entry. Allowed types are GAP and CPP\n";
    throw TerminalException{1};
  };
  print_stderr_stdout_file(FileGroup, f_print);
}

//
// group combinatorial algorithms
//

void f_print(std::vector<int> const &V, std::string const &estr) {
  std::cerr << estr << " =";
  for (auto &val : V)
    std::cerr << " " << val;
  std::cerr << "\n";
}

template <typename Tgroup>
std::vector<int> OrbitIntersection(Tgroup const &TheGRP,
                                   std::vector<int> const &gList) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = TheGRP.n_act();
  std::vector<int> rList = gList;
  auto f_sum = [&]() -> Tidx {
    Tidx eSum = 0;
    for (Tidx i = 0; i < n; i++)
      eSum += rList[i];
    return eSum;
  };
  auto LGen = TheGRP.GeneratorsOfGroup();
  Tidx eSum = f_sum();
  f_print(rList, "input(rList)");
  while (true) {
    for (Tidx i = 0; i < n; i++) {
      if (rList[i] == 0) {
        for (auto &eGen : LGen) {
          int j = OnPoints(i, eGen);
          rList[j] = 0;
        }
      }
    }
    f_print(rList, "iter(rList)");
    Tidx eSumNew = f_sum();
    if (eSum == eSumNew)
      break;
    eSum = eSumNew;
  }
  f_print(rList, "returning(rList)");
  return rList;
}

template <typename Tgroup>
std::vector<int> OrbitUnion(Tgroup const &TheGRP,
                            std::vector<int> const &gList) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = TheGRP.n_act();
  std::vector<int> gListB(n);
  for (Tidx i = 0; i < n; i++)
    gListB[i] = 1 - gList[i];
  f_print(gList, "OrbitUnion(gList)");
  f_print(gListB, "OrbitUnion(gListB)");
  std::vector<int> rListB = OrbitIntersection(TheGRP, gListB);
  for (Tidx i = 0; i < n; i++)
    rListB[i] = 1 - rListB[i];
  return rListB;
}

template <typename Tgroup>
Face OrbitIntersection(Tgroup const &GRP, Face const &gList) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::vector<Telt> LGen = GRP.GeneratorsOfGroup();
  Tidx n = GRP.n_act();
  Face rList = gList;
  size_t eSum = rList.count();
  while (true) {
    for (Tidx i = 0; i < n; i++) {
      if (rList[i] == 0) {
        for (auto &eGen : LGen) {
          Tidx j = OnPoints(i, eGen);
          rList[j] = 0;
        }
      }
    }
    size_t eSumNew = rList.count();
    if (eSum == eSumNew)
      break;
    eSum = eSumNew;
  }
  return rList;
}

template <typename Tgroup>
Face OrbitUnion(Tgroup const &GRP, Face const &gList) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = GRP.n_act();
  Face gListB(n);
  for (Tidx i = 0; i < n; i++)
    gListB[i] = 1 - gList[i];
  Face rListB = OrbitIntersection(GRP, gListB);
  for (Tidx i = 0; i < n; i++)
    rListB[i] = 1 - rListB[i];
  return rListB;
}

//
// Several building of new groups.
//

template <typename Tidx> struct ReducingArray {
  size_t nb;
  std::vector<Tidx> ListPositionRev;
  std::vector<Tidx> ListPosition;
};

template <typename Tidx>
ReducingArray<Tidx> GetReducingArrayFace(Face const &eList) {
  Tidx size = eList.size();
  size_t nb = eList.count();
  Tidx nb_i = nb;
  if (nb_i == 0) {
    std::cerr << "Call of GetReducingArrayFace with 0 points\n";
    throw TerminalException{1};
  }
  std::vector<Tidx> ListPositionRev(size, -1);
  boost::dynamic_bitset<>::size_type aRow = eList.find_first();
  std::vector<Tidx> ListPosition(nb);
  for (Tidx iRow = 0; iRow < nb_i; iRow++) {
    ListPositionRev[aRow] = iRow;
    ListPosition[iRow] = Tidx(aRow);
    aRow = eList.find_next(aRow);
  }
  return {nb, std::move(ListPositionRev), std::move(ListPosition)};
}

template <typename Telt>
Telt SingleElementReduction(Telt const &eElt,
                            ReducingArray<typename Telt::Tidx> const &ra) {
  using Tidx = typename Telt::Tidx;
  std::vector<Tidx> v(ra.nb);
  for (size_t i = 0; i < ra.nb; i++) {
    Tidx eVal1 = ra.ListPosition[i];
    Tidx eVal2 = OnPoints(eVal1, eElt);
    Tidx eVal3 = ra.ListPositionRev[eVal2];
    v[i] = eVal3;
  }
  return Telt(std::move(v));
}

template <typename Telt>
Telt ReduceElementActionFace(Telt const &eElt, Face const &eList) {
  using Tidx = typename Telt::Tidx;
  ReducingArray<Tidx> ra = GetReducingArrayFace<Tidx>(eList);
  return SingleElementReduction(eElt, ra);
}

template <typename Tgroup>
Tgroup ReducedGroupActionRa(Tgroup const &TheGRP, Face const &eList) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  ReducingArray<Tidx> ra = GetReducingArrayFace<Tidx>(eList);
  std::vector<Telt> ListGen;
  for (auto &eGen : TheGRP.GeneratorsOfGroup()) {
    Telt eGenRed = SingleElementReduction(eGen, ra);
    ListGen.emplace_back(std::move(eGenRed));
  }
  return Tgroup(ListGen, ra.nb);
}

template <typename Tidx>
ReducingArray<Tidx> GetReducingArrayVect(std::vector<Tidx> const &eList,
                                         Tidx const &size) {
  size_t nb = eList.size();
  Tidx nb_i = nb;
  if (nb_i == 0) {
    std::cerr << "Call of GetReducingArrayVect with 0 points\n";
    throw TerminalException{1};
  }
  std::vector<Tidx> ListPositionRev(size, -1);
  for (Tidx pos = 0; pos < nb_i; pos++) {
    Tidx val = eList[pos];
    ListPositionRev[val] = pos;
  }
  return {nb, std::move(ListPositionRev), eList};
}

template <typename Tgroup>
Tgroup
ReducedGroupActionRa(Tgroup const &TheGRP,
                     ReducingArray<typename Tgroup::Telt::Tidx> const &ra) {
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> ListGen;
  for (auto &eGen : TheGRP.GeneratorsOfGroup()) {
    Telt eGenRed = SingleElementReduction(eGen, ra);
    ListGen.emplace_back(std::move(eGenRed));
  }
  return Tgroup(ListGen, ra.nb);
}

template <typename Tgroup>
Tgroup ReducedGroupActionFace(Tgroup const &TheGRP, Face const &eList) {
  using Tidx = typename Tgroup::Telt::Tidx;
  ReducingArray<Tidx> ra = GetReducingArrayFace<Tidx>(eList);
  return ReducedGroupActionRa<Tgroup>(TheGRP, ra);
}

template <typename Tgroup>
Tgroup
ReducedGroupActionVect(Tgroup const &TheGRP,
                       std::vector<typename Tgroup::Telt::Tidx> const &eList) {
  using Tidx = typename Tgroup::Telt::Tidx;
  Tidx size = TheGRP.n_act();
  ReducingArray<Tidx> ra = GetReducingArrayVect<Tidx>(eList, size);
  return ReducedGroupActionRa<Tgroup>(TheGRP, ra);
}

template <typename Tgroup>
Tgroup ConjugateGroup(Tgroup const &TheGRP,
                      typename Tgroup::Telt const &ePerm) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tidx n = TheGRP.n_act();
  Telt ePermR = ~ePerm;
  std::vector<Telt> ListGen;
  for (auto &eGen : TheGRP.GeneratorsOfGroup()) {
    std::vector<Tidx> v(n);
    for (Tidx i = 0; i < n; i++) {
      Tidx eVal1 = OnPoints(i, ePermR);
      Tidx eVal2 = OnPoints(eVal1, eGen);
      Tidx eVal3 = OnPoints(eVal2, ePerm);
      v[i] = eVal3;
    }
    ListGen.emplace_back(std::move(Telt(std::move(v))));
  }
  return Tgroup(ListGen, n);
}

//
// Some enumeration code
//

template <typename Tgroup>
void GROUP_FuncInsertInSet(Tgroup const &TheGRP, Face const &eList,
                           vectface &ListListSet) {
  using Telt = typename Tgroup::Telt;
  for (auto &fList : ListListSet) {
    std::optional<Telt> test = TheGRP.RepresentativeAction_OnSets(eList, fList);
    if (test)
      return;
  }
  ListListSet.push_back(eList);
}

template <typename Tgroup>
void GROUP_FuncInsertInSet_UseInv(Tgroup const &TheGRP, Face const &eList,
                                  std::vector<int> const &eInv,
                                  vectface &ListSet,
                                  std::vector<std::vector<int>> &ListInv) {
  using Telt = typename Tgroup::Telt;
  size_t nb = ListSet.size();
  for (size_t iList = 0; iList < nb; iList++)
    if (eInv == ListInv[iList]) {
      std::optional<Telt> test =
          TheGRP.RepresentativeAction_OnSets(eList, ListSet[iList]);
      if (test)
        return;
    }
  ListSet.push_back(eList);
  ListInv.push_back(eInv);
}

//
// Some combinatorial algorithms using only the generators of the group.
//

template <typename Telt>
std::vector<int> ComputeFullOrbitPoint(const size_t &n,
                                       const std::vector<Telt> &ListGen,
                                       int const &ePoint) {
  IntegerSubsetStorage Vorb = VSLT_InitializeStorage(n);
  IntegerSubsetStorage Vactive = VSLT_InitializeStorage(n);
  std::vector<int> eList;
  VSLT_StoreValue(Vactive, ePoint);
  while (true) {
    if (VSLT_IsEmpty(Vactive))
      break;
    int TheFirst = VSLT_TheFirstPosition(Vactive);
    VSLT_RemoveValue(Vactive, TheFirst);
    VSLT_StoreValue(Vorb, TheFirst);
    eList.push_back(TheFirst);
    for (auto &eGen : ListGen) {
      int NewPt = OnPoints(TheFirst, eGen);
      if (!VSLT_IsItInSubset(Vorb, NewPt) && !VSLT_IsItInSubset(Vactive, NewPt))
        VSLT_StoreValue(Vactive, NewPt);
    }
  }
  return eList;
}

template <typename Telt>
vectface DecomposeOrbitPoint_Kernel(const std::vector<Telt> &LGen,
                                    Face const &eList) {
  size_t nbPoint = eList.size();
  IntegerSubsetStorage Vlist = VSLT_InitializeStorage(nbPoint);
  vectface ListOrb(nbPoint);
  size_t len = eList.count();
  boost::dynamic_bitset<>::size_type aRow = eList.find_first();
  for (size_t i = 0; i < len; i++) {
    VSLT_StoreValue(Vlist, size_t(aRow));
    aRow = eList.find_next(aRow);
  }
  while (true) {
    if (VSLT_IsEmpty(Vlist))
      break;
    size_t TheFirst = VSLT_TheFirstPosition(Vlist);
    std::vector<int> eOrb = ComputeFullOrbitPoint(nbPoint, LGen, TheFirst);
    Face vectOrb(nbPoint);
    for (auto &ePt : eOrb) {
      vectOrb[ePt] = 1;
      VSLT_RemoveValue(Vlist, ePt);
    }
    ListOrb.push_back(vectOrb);
  }
  return ListOrb;
}

template <typename Tgroup>
vectface DecomposeOrbitPoint(Tgroup const &TheGRP, Face const &eList) {
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> LGen = TheGRP.GeneratorsOfGroup();
  return DecomposeOrbitPoint_Kernel(LGen, eList);
}

template <typename Tgroup>
vectface DecomposeOrbitPoint_Full(Tgroup const &TheGRP) {
  size_t n = TheGRP.n_act();
  Face eList(n);
  for (size_t i = 0; i < n; i++)
    eList[i] = 1;
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> LGen = TheGRP.GeneratorsOfGroup();
  return DecomposeOrbitPoint_Kernel(LGen, eList);
}

template <typename Tgroup>
std::vector<size_t> DecomposeOrbitPoint_FullRepr(Tgroup const &TheGRP) {
  vectface vf = DecomposeOrbitPoint_Full(TheGRP);
  std::vector<size_t> l_idx;
  for (auto &eFace : vf) {
    boost::dynamic_bitset<>::size_type aRow = eFace.find_first();
    l_idx.push_back(aRow);
  }
  return l_idx;
}

template <typename Telt>
vectface DecomposeOrbitPoint_KernelFull(const size_t &n,
                                        const std::vector<Telt> &LGen) {
  Face eList(n);
  for (size_t i = 0; i < n; i++)
    eList[i] = 1;
  return DecomposeOrbitPoint_Kernel(LGen, eList);
}

template <typename Tobj, typename Tgen>
std::vector<Tobj> OrbitSplittingGeneralized(
    std::vector<Tobj> const &PreListTotal, std::vector<Tgen> const &ListGen,
    std::function<Tobj(Tobj const &, Tgen const &)> const &TheAct) {
  std::vector<Tobj> TheReturn;
  std::unordered_set<Tobj> ListTotal;
  for (auto eObj : PreListTotal)
    ListTotal.insert(eObj);
  while (true) {
    auto iter = ListTotal.begin();
    if (iter == ListTotal.end())
      break;
    Tobj eObj = *iter;
    TheReturn.push_back(eObj);
    std::unordered_set<Tobj> Additional;
    Additional.insert(eObj);
    ListTotal.erase(eObj);
    std::unordered_set<Tobj> SingleOrbit;
    while (true) {
      std::unordered_set<Tobj> NewElts;
      for (auto const &gObj : Additional) {
        for (auto const &eGen : ListGen) {
          Tobj fObj = TheAct(gObj, eGen);
          if (SingleOrbit.count(fObj) == 0 && Additional.count(fObj) == 0) {
            if (NewElts.count(fObj) == 0) {
#ifdef DEBUG_GROUP
              if (ListTotal.count(fObj) == 0) {
                std::cerr << "Orbit do not match, PANIC!!!\n";
                throw TerminalException{1};
              }
#endif
              NewElts.insert(fObj);
              ListTotal.erase(fObj);
            }
          }
        }
      }
      for (auto &uSet : Additional)
        SingleOrbit.insert(uSet);
      if (NewElts.size() == 0)
        break;
      Additional = std::move(NewElts);
    }
  }
  return TheReturn;
}

template <typename Telt, typename Tobj, typename Tact>
std::vector<std::pair<Tobj, Telt>>
OrbitWithRepresentative(const Telt &id, std::vector<Telt> const &ListGen,
                        Tobj const &x, Tact act) {
  std::vector<std::pair<Tobj, Telt>> ListObj{{x, id}};
  std::unordered_set<Tobj> SetObj{x};
  size_t curr_pos = 0;
  while (true) {
    size_t len = ListObj.size();
    if (curr_pos == len)
      break;
    for (size_t u = curr_pos; u < len; u++) {
      for (auto &eElt : ListGen) {
        Tobj eImg = act(ListObj[u].first, eElt);
        if (SetObj.count(eImg) == 0) {
          Telt NewElt = ListObj[u].second * eElt;
          ListObj.push_back({eImg, NewElt});
          SetObj.insert(eImg);
        }
      }
    }
    curr_pos = len;
  }
  return ListObj;
}

template <typename Telt> Face OnFace(Face const &eSet, Telt const &eElt) {
  using Tidx = typename Telt::Tidx;
  size_t nbExt = eSet.size();
  Face fSet(nbExt);
  boost::dynamic_bitset<>::size_type pos = eSet.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    Tidx jExt = eElt.at(pos);
    fSet[jExt] = 1;
    pos = eSet.find_next(pos);
  }
  return fSet;
}

template <typename Telt>
void OnFace_inplace(Face &fSet, Face const &eSet, Telt const &eElt) {
  using Tidx = typename Telt::Tidx;
  fSet.reset();
  boost::dynamic_bitset<>::size_type pos = eSet.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    Tidx jExt = eElt.at(pos);
    fSet[jExt] = 1;
    pos = eSet.find_next(pos);
  }
}

template <typename Telt>
vectface OrbitFace(const Face &f, const std::vector<Telt> &LGen) {
  size_t len = f.size();
  std::unordered_set<Face> set_f;
  vectface vf(len);
  Face gw(len);
  std::vector<uint8_t> status;
  auto func_insert = [&](const Face &f) -> void {
    if (set_f.find(f) == set_f.end()) {
      vf.push_back(f);
      set_f.insert(f);
      status.push_back(0);
    }
  };
  func_insert(f);
  while (true) {
    bool is_finished = true;
    size_t len = vf.size();
    for (size_t i = 0; i < len; i++) {
      if (status[i] == 0) {
        is_finished = false;
        status[i] = 1;
        Face fw = vf[i];
        for (auto &eGen : LGen) {
          OnFace_inplace(gw, fw, eGen);
          func_insert(gw);
        }
      }
    }
    if (is_finished)
      break;
  }
  return vf;
}

template <typename Tgroup, typename T>
std::vector<std::pair<Face, T>>
OrbitSplittingMap(std::vector<std::pair<Face, T>> &PreListTotal,
                  Tgroup const &TheGRP) {
  using Telt = typename Tgroup::Telt;
  std::unordered_map<Face, T> ListTotal;
  for (auto &ePair : PreListTotal)
    ListTotal[std::move(ePair.first)] = std::move(ePair.second);
  std::vector<std::pair<Face, T>> ListReturn;
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  Face fSet(TheGRP.n_act());
  while (true) {
    typename std::unordered_map<Face, T>::iterator iter = ListTotal.begin();
    if (iter == ListTotal.end())
      break;
    Face eSet = iter->first;
    T val = std::move(iter->second);
    std::unordered_set<Face> Additional{eSet};
    ListTotal.erase(eSet);
    std::unordered_set<Face> SingleOrbit;
    while (true) {
      std::unordered_set<Face> NewElts;
      for (auto const &gSet : Additional) {
        for (auto const &eGen : ListGen) {
          OnFace_inplace(fSet, gSet, eGen);
          if (SingleOrbit.count(fSet) == 0 && Additional.count(fSet) == 0) {
            if (NewElts.count(fSet) == 0) {
#ifdef DEBUG_GROUP
              if (ListTotal.count(fSet) == 0) {
                std::cerr << "Orbit do not matched, PANIC!!!\n";
                throw TerminalException{1};
              }
#endif
              NewElts.insert(fSet);
              ListTotal.erase(fSet);
            }
          }
        }
      }
      for (auto &uSet : Additional)
        SingleOrbit.insert(uSet);
      if (NewElts.size() == 0)
        break;
      Additional = std::move(NewElts);
    }
    ListReturn.push_back({std::move(eSet), std::move(val)});
  }
  return ListReturn;
}

template <typename Tgroup, typename T_hash_set>
vectface OrbitSplittingSet_T(T_hash_set &ListTotal, Tgroup const &TheGRP) {
  using Telt = typename Tgroup::Telt;
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  size_t n = TheGRP.n_act();
  Face fSet(n);
  vectface vf_ret(n);
  while (true) {
    typename T_hash_set::iterator iter = ListTotal.begin();
    if (iter == ListTotal.end())
      break;
    Face eSet = *iter;
    T_hash_set SingleOrbit;
    vectface vf(n);
    size_t total_len = 0;
    auto f_insert = [&](Face const &f) -> void {
      if (SingleOrbit.count(f) == 0) {
        SingleOrbit.insert(f);
        ListTotal.erase(f);
        vf.push_back(f);
        total_len++;
      }
    };
    f_insert(eSet);
    size_t pos = 0;
    while (true) {
      if (pos == total_len) {
        break;
      }
      size_t curr_len = total_len;
      for (size_t idx = pos; idx < curr_len; idx++) {
        Face f = vf[idx];
        for (auto const &eGen : ListGen) {
          OnFace_inplace(fSet, f, eGen);
          f_insert(fSet);
        }
      }
      pos = curr_len;
    }
    vf_ret.push_back(eSet);
  }
  return vf_ret;
}

template <typename Tgroup, typename F>
void OrbitSplittingSet_Kernel(vectface const &PreListTotal,
                              Tgroup const &TheGRP, F f) {
  using Telt = typename Tgroup::Telt;
  std::unordered_set<Face> ListTotal;
  for (auto eFace : PreListTotal)
    ListTotal.insert(eFace);
  std::vector<Telt> ListGen = TheGRP.GeneratorsOfGroup();
  Face fSet(TheGRP.n_act());
  size_t n = PreListTotal.get_n();
  //  size_t total_size = PreListTotal.size();
  //  std::cerr << "|ListTotal|=" << ListTotal.size() << " |PreListTotal|=" <<
  //  PreListTotal.size() << "\n"; size_t tot_sum = 0;
  while (true) {
    std::unordered_set<Face>::iterator iter = ListTotal.begin();
    if (iter == ListTotal.end())
      break;
    Face eSet = *iter;
    std::unordered_set<Face> SingleOrbit;
    vectface vf(n);
    size_t total_len = 0;
    auto f_insert = [&](Face const &f) -> void {
      if (SingleOrbit.count(f) == 0) {
        SingleOrbit.insert(f);
        ListTotal.erase(f);
        vf.push_back(f);
        //        std::cerr << "erasing f=" << f << "\n";
        total_len++;
      }
    };
    f_insert(eSet);
    size_t pos = 0;
    while (true) {
      //      std::cerr << "pos=" << pos << " total_len=" << total_len << "\n";
      if (pos == total_len) {
        break;
      }
      size_t curr_len = total_len;
      for (size_t idx = pos; idx < curr_len; idx++) {
        Face f = vf[idx];
        for (auto const &eGen : ListGen) {
          OnFace_inplace(fSet, f, eGen);
          //          std::cerr << "f=" << f << "\n";
          //          std::cerr << "fSet=" << f << "\n";
          f_insert(fSet);
        }
      }
      pos = curr_len;
    }
    //    std::cerr << "   |SingleOrbit|=" << SingleOrbit.size() << "
    //    |PreListTotal|=" << PreListTotal.size() << " |ListTotal|=" <<
    //    ListTotal.size() << "\n"; tot_sum += SingleOrbit.size();
    f(eSet, SingleOrbit);
  }
  //  std::cerr << "tot_sum=" << tot_sum << " total_size=" << total_size <<
  //  "\n";
}

template <typename Tgroup>
vectface OrbitSplittingSet(vectface const &PreListTotal, Tgroup const &TheGRP) {
  vectface TheReturn(TheGRP.n_act());
  auto f = [&](Face const &eSet,
               [[maybe_unused]] std::unordered_set<Face> const &SingleOrbit)
      -> void { TheReturn.push_back(eSet); };
  OrbitSplittingSet_Kernel(PreListTotal, TheGRP, f);
  return TheReturn;
}

template <typename Tgroup>
vectface OrbitSplittingSet_GetMinimalOrbit(vectface const &PreListTotal,
                                           Tgroup const &TheGRP) {
  using Tidx = typename Tgroup::Telt::Tidx;
  Tidx len = TheGRP.n_act();
  vectface TheReturn(len);
  Face TheMin;
  bool HasMin = false;
  std::cerr << "OrbitSizes =";
  auto f = [&]([[maybe_unused]] Face const &eSet,
               std::unordered_set<Face> const &SingleOrbit) -> void {
    //    std::cerr << "f : begin\n";
    vectface orbit(len);
    Face minF;
    bool IsFirst = true;
    std::cerr << " " << SingleOrbit.size();
    for (auto &uSet : SingleOrbit) {
      orbit.push_back(uSet);
      if (IsFirst) {
        minF = uSet;
        IsFirst = false;
      } else {
        if (uSet < minF)
          minF = uSet;
      }
    }
    //    std::cerr << "f : orbit and minF built\n";
    // Now doing the comparison with existing data
    auto set_return = [&]() -> void {
      TheReturn = std::move(orbit);
      TheMin = minF;
    };
    if (!HasMin) {
      set_return();
      HasMin = true;
    } else {
      if (orbit.size() < TheReturn.size()) {
        set_return();
      } else {
        if (orbit.size() == TheReturn.size()) {
          if (minF < TheMin) {
            set_return();
          }
        }
      }
    }
    //    std::cerr << "f : end\n";
  };
  OrbitSplittingSet_Kernel(PreListTotal, TheGRP, f);
  std::cerr << " -- ";
  return TheReturn;
}

// Test if f1 is a subset of f2
bool is_subset(Face const &f1, Face const &f2) {
  boost::dynamic_bitset<>::size_type pos = f1.find_first();
  while (pos != boost::dynamic_bitset<>::npos) {
    if (f2[pos] == 0)
      return false;
    pos = f1.find_next(pos);
  }
  return true;
}

/*
  This is a combinatorial algorithm for finding all group elements g such that
  set2 \subset g . set1 We could probably do better with double coset
  decompositions.
 */
template <typename Tgroup>
std::vector<std::pair<Face, typename Tgroup::Telt>>
FindContainingOrbit(Tgroup const &GRP_ext, Face const &set1, Face const &set2) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  Tgroup stab1 = GRP_ext.Stabilizer_OnSets(set1);
  Tgroup stab2 = GRP_ext.Stabilizer_OnSets(set2);
#ifdef DEBUG_GROUP
  std::cerr << "|stab1=" << stab1.size() << " |stab2|=" << stab2.size() << "\n";
#endif
  std::vector<Telt> LGen = GRP_ext.GeneratorsOfGroup();
  Tidx n = GRP_ext.n_act();
  Telt id(n);
  auto f_act = [](const Face &x, const Telt &u) -> Face {
    return OnFace(x, u);
  };
  std::vector<std::pair<Face, Telt>> LPair =
      OrbitWithRepresentative(id, LGen, set1, f_act);
#ifdef DEBUG_GROUP
  std::cerr << "|LPair|=" << LPair.size() << "\n";
#endif
  std::unordered_map<Face, Telt> map_face_elt;
  vectface list_face(n);
  for (auto &ePair : LPair) {
    if (is_subset(set2, ePair.first)) {
      map_face_elt[ePair.first] = ePair.second;
      list_face.push_back(ePair.first);
    }
  }
#ifdef DEBUG_GROUP
  std::cerr << "|list_face|=" << list_face.size() << "\n";
#endif
  vectface vfRepr = OrbitSplittingSet(list_face, stab2);
  std::vector<std::pair<Face, typename Tgroup::Telt>> list_ret;
  for (auto &eRepr : vfRepr) {
    Telt eElt = map_face_elt[eRepr];
    list_ret.push_back({eRepr, eElt});
  }
  return list_ret;
}

// This is a remake of the "Partition<Tidx>" in partition.h of
// permutalib. But that code is super sensitive and we do not
// want to manipulate it.
template<typename Telt>
struct PartitionStorage {
private:
  using Tidx = typename Telt::Tidx;
public:
  std::ostream& os;
  // The list of points in the list
  std::vector<Tidx> points;
  // The indices of the first points.
  std::vector<Tidx> firsts;
  // The lengths of the part.
  std::vector<Tidx> lengths;
  // The cell to which it belongs
  std::vector<Tidx> cellno;
  // The position in the cell to which it belongs
  std::vector<Tidx> cellpos;
  // number of parts
  Tidx n_part;
  // number of elements
  Tidx n_elt;
  // scratch arrays
  std::vector<Tidx> scratch1;
  std::vector<Tidx> scratch2;
  std::vector<Tidx> scratch3;
  std::vector<Tidx> scratch4;
  PartitionStorage(Face const& f, std::ostream& _os) : os(_os) {
    size_t siz_in = f.count();
    size_t siz_out = f.size() - siz_in;
    firsts.push_back(0);
    firsts.push_back(siz_in);
    lengths.push_back(siz_in);
    lengths.push_back(siz_out);
    Tidx cell_in;
    Tidx cell_out;
    if (f.count() > 0 && f.count() < f.size()) {
      n_part = 2;
      cell_in = 0;
      cell_out = 1;
    } else {
      n_part = 1;
      // All goes to the same place. Put 0 all the way.
      cell_in = 0;
      cell_out = 0;
    }
    n_elt = f.size();
#ifdef DEBUG_GROUP
    os << "GRPFCT: PartitionStorage(Const), n_elt=" << n_elt << "\n";
#endif
    Tidx pos_in = 0;
    Tidx pos_out = 0;
    points = std::vector<Tidx>(n_elt);
    cellno = std::vector<Tidx>(n_elt);
    cellpos = std::vector<Tidx>(n_elt);
    for (Tidx i=0; i<n_elt; i++) {
      if (f[i] == 1) {
        points[pos_in] = i;
        cellno[i] = cell_in;
        cellpos[i] = pos_in;
        pos_in += 1;
      } else {
        points[pos_out + siz_in] = i;
        cellno[i] = cell_out;
        cellpos[i] = pos_out;
        pos_out += 1;
      }
    }
    scratch1 = std::vector<Tidx>(n_elt);
    scratch2 = std::vector<Tidx>(n_elt);
    scratch3 = std::vector<Tidx>(n_elt, 0);
    scratch4 = std::vector<Tidx>(n_elt);
  }
  bool RefinePartitionByElement(Telt const& g) {
    for (Tidx ipart=0; ipart<n_part; ipart++) {
      for (Tidx j=0; j<n_part; j++) {
        scratch1[j] = 0;
      }
      Tidx pos_first = firsts[ipart];
      Tidx len = lengths[ipart];
      for (Tidx u=0; u<len; u++) {
        Tidx x = points[pos_first + u];
        Tidx y = g.at(x);
        Tidx i_cell = cellno[y];
        Tidx first = firsts[i_cell];
        Tidx n_belong = scratch1[i_cell];
        scratch2[first + n_belong] = y;
        n_belong += 1;
        scratch1[i_cell] = n_belong;
      }
      Tidx n_part_more = 0;
      for (Tidx i_cell=0; i_cell<n_part; i_cell++) {
        Tidx siz_part = scratch1[i_cell];
        if (siz_part != 0 && siz_part != lengths[i_cell]) {
          Tidx len = lengths[i_cell];
          Tidx first = firsts[i_cell];
          Tidx siz_out_part = len - siz_part;
          for (Tidx i=0; i<siz_part; i++) {
            Tidx x = scratch2[first + i];
            Tidx pos = cellpos[x];
            scratch3[pos] = 1;
          }
          for (Tidx i=0; i<len; i++) {
            scratch4[i] = points[first + i];
          }
          Tidx pos_in = 0;
          Tidx pos_out = 0;
          for (Tidx u=0; u<len; u++) {
            Tidx x = scratch4[u];
            if (scratch3[u] == 1) {
              points[first + pos_in] = x;
              cellno[x] = i_cell;
              cellpos[x] = pos_in;
              pos_in += 1;
            } else {
              points[first + siz_part + pos_out] = x;
              cellno[x] = n_part + n_part_more;
              cellpos[x] = pos_out;
              pos_out += 1;
            }
          }
          lengths[i_cell] = siz_part;
          lengths.push_back(siz_out_part);
          firsts.push_back(firsts[i_cell] + siz_part);
          for (Tidx u=0; u<len; u++) {
            scratch3[u] = 0;
          }
          n_part_more += 1;
        }
      }
      n_part += n_part_more;
      if (n_part_more > 0) {
        return false;
      }
    }
    return true;
  }
  void RefinePartitionByListElt(std::vector<Telt> const& list_gen) {
    while(true) {
      bool IsFinished = true;
      for (auto & e_gen: list_gen) {
        bool test = RefinePartitionByElement(e_gen);
        if (!test) {
          IsFinished = false;
        }
      }
      if (IsFinished) {
        break;
      }
    }
  }
  Face map_face(Face const& f) const {
#ifdef SANITY_CHECK_GROUP_FCT
    std::vector<Tidx> n_occur(n_part, 0);
#endif
    Face fret(n_part);
    for (Tidx i=0; i<n_elt; i++) {
      if (f[i] == 1) {
        Tidx i_cell = cellno[i];
#ifdef SANITY_CHECK_GROUP_FCT
        n_occur[i_cell] += 1;
#endif
        fret[i_cell] = 1;
      }
    }
#ifdef SANITY_CHECK_GROUP_FCT
    for (Tidx i_part=0; i_part<n_part; i_part++) {
      if (n_occur[i_part] != 0 && n_occur[i_part] != lengths[i_part]) {
        std::cerr << "The size is not what we expect\n";
        throw TerminalException{1};
      }
    }
#endif
    return fret;
  }
  std::optional<Face> map_face_opt(Face const& f) const {
    std::vector<Tidx> n_occur(n_part, 0);
    Face fret(n_part);
    for (Tidx i=0; i<n_elt; i++) {
      if (f[i] == 1) {
        Tidx i_cell = cellno[i];
        n_occur[i_cell] += 1;
        fret[i_cell] = 1;
      }
    }
    for (Tidx i_part=0; i_part<n_part; i_part++) {
      if (n_occur[i_part] != 0 && n_occur[i_part] != lengths[i_part]) {
        return {};
      }
    }
    return fret;
  }
  Telt map_permutation(Telt const& g) const {
    std::vector<Tidx> eList(n_part);
    for (Tidx i_cell=0; i_cell<n_part; i_cell++) {
      Tidx x = points[firsts[i_cell]];
      Tidx y = g.at(x);
      Tidx j_cell = cellno[y];
      eList[i_cell] = j_cell;
    }
    return Telt(eList);
  }
};



// clang-format off
#endif  // SRC_GROUP_GRP_GROUPFCT_H_
// clang-format on
