#ifndef SRC_POLY_POLY_KSKELETTON_H_
#define SRC_POLY_POLY_KSKELETTON_H_

#include "GRP_GroupFct.h"
#include "MAT_Matrix.h"
#include "Namelist.h"
#include "POLY_LinearProgramming.h"
#include "POLY_PolytopeFct.h"
#include "Temp_PolytopeEquiStab.h"
#include <map>
#include <utility>
#include <vector>
#include <string>
#include <limits>

// We follow here the conventions of SPAN_face_LinearProgramming
// in Kskeleton.g for the computation.
template <typename T, typename Tgroup>
vectface
SPAN_face_LinearProgramming(Face const &face_fac, Tgroup const &StabFace_fac,
                            MyMatrix<T> const &FAC, Tgroup const &GRP_fac) {
  int nbRow = FAC.rows();
  size_t n_act_stabface = StabFace_fac.n_act();
  size_t n_act_fullgrp = GRP_fac.n_act();
  if (size_t(nbRow) != n_act_stabface || size_t(nbRow) != n_act_fullgrp) {
    std::cerr << "Inconsistency in nbRow, n_act_stabface, n_act_fullgrp\n";
    throw TerminalException{1};
  }
  int nbCol = FAC.cols();
  vectface TheReturn(nbRow);
  MyMatrix<T> eMatId = IdentityMat<T>(nbCol);
  Face Treated(nbRow);
  size_t sizFace = face_fac.count();
  MyMatrix<T> TestMat(sizFace + 1, nbCol);
  boost::dynamic_bitset<>::size_type ePt = face_fac.find_first();
  for (size_t iRow = 0; iRow < sizFace; iRow++) {
    Treated[ePt] = 1;
    TestMat.row(iRow) = FAC.row(ePt);
    ePt = face_fac.find_next(ePt);
  }
  for (int iRow = 0; iRow < nbRow; iRow++)
    if (Treated[iRow] == 0) {
      TestMat.row(sizFace) = FAC.row(iRow);
      SelectionRowCol<T> eSelect = TMat_SelectRowCol(TestMat);
      MyMatrix<T> NSP = eSelect.NSP;
      int nbEqua = NSP.rows();
      Face eCand(nbRow);
      Face eCandCompl(nbRow);
      Face gList(nbRow);
      for (int jRow = 0; jRow < nbRow; jRow++) {
        auto get_test = [&]() -> bool {
          for (int iEqua = 0; iEqua < nbEqua; iEqua++) {
            T eSum = 0;
            for (int iCol = 0; iCol < nbCol; iCol++)
              eSum += FAC(jRow, iCol) * NSP(iEqua, iCol);
            if (eSum != 0)
              return false;
          }
          return true;
        };
        if (get_test()) {
          eCand[jRow] = 1;
          gList[jRow] = 1;
        } else {
          eCandCompl[jRow] = 1;
        }
      }
      Face rList = OrbitUnion(StabFace_fac, gList);
      for (int jRow = 0; jRow < nbRow; jRow++)
        if (rList[jRow] == 1)
          Treated[jRow] = 1;
      Tgroup TheStab = GRP_fac.Stabilizer_OnSets(eCand);
      vectface ListOrb = DecomposeOrbitPoint(TheStab, eCand);
      size_t nbOrb = ListOrb.size();
      MyMatrix<T> ListVectSpann(nbOrb, nbCol);
      for (size_t iOrb = 0; iOrb < nbOrb; iOrb++) {
        Face eOrb = ListOrb[iOrb];
        MyVector<T> eVec = SumMatrixLineSubset(FAC, eOrb);
        AssignMatrixRow(ListVectSpann, iOrb, eVec);
      }
      MyMatrix<T> BasisSpann = RowReduction(ListVectSpann);
      int nbRowSpann = BasisSpann.rows();
      int LPdim = nbCol - nbRowSpann;
      MyMatrix<T> PreTheTot = Concatenate(BasisSpann, eMatId);
      MyMatrix<T> TheTot = RowReduction(PreTheTot);
      // the complement
      vectface ListOrbCompl = DecomposeOrbitPoint(TheStab, eCandCompl);
      size_t nbOrbCompl = ListOrbCompl.size();
      MyMatrix<T> PreListVectors(nbOrbCompl, nbCol);
      for (size_t iOrb = 0; iOrb < nbOrbCompl; iOrb++) {
        Face eOrb = ListOrbCompl[iOrb];
        MyVector<T> eVec = SumMatrixLineSubset(FAC, eOrb);
        AssignMatrixRow(PreListVectors, iOrb, eVec);
      }
      MyMatrix<T> TheTotInv = Inverse(TheTot);
      MyMatrix<T> PreListVectorsB = PreListVectors * TheTotInv;
      MyMatrix<T> ListVectors(nbOrbCompl, LPdim);
      for (size_t jRow = 0; jRow < nbOrbCompl; jRow++)
        for (int iCol = 0; iCol < LPdim; iCol++) {
          int jCol = iCol + nbRowSpann;
          ListVectors(jRow, iCol) = PreListVectorsB(jRow, jCol);
        }
      bool eTestExist = TestPositiveRelationSimple(ListVectors);
      if (eTestExist == 0)
        TheReturn.push_back(eCand);
    }
  return TheReturn;
}

/*
  FAC is the list of facets.
  EXT is the list of vertices (only used for the testing of faces)
  - face is the face as encoded in FAC.
  - StabFace and FullGRP are also encoded in FAC.
  RankFace must be the rank of the face in the comple.
  So, for example if face has just 1 element then RankFace=1.
  ---
  The elements being returned in vectface all have more entries.
  It could be just 1 more or more.
 */
template <typename T, typename Tgroup, typename Ffilt>
vectface SPAN_face_ExtremeRays_F(Face const &face_fac,
                                 Tgroup const &StabFace_fac,
                                 int const &RankFace, const Face &extfac_incd,
                                 MyMatrix<T> const &FAC, MyMatrix<T> const &EXT,
                                 Ffilt f_filt) {
  int nbFac = FAC.rows();
  size_t n_act_stabface = StabFace_fac.n_act();
  if (size_t(nbFac) != n_act_stabface) {
    std::cerr << "Inconsistency in nbFac, n_act_stabface\n";
    std::cerr << "nbFac=" << nbFac << " n_act_stabface=" << n_act_stabface
              << "\n";
    throw TerminalException{1};
  }
  int nbCol = FAC.cols();
  vectface TheReturn(nbFac);
  Face Treated(nbFac);
  size_t sizFace = face_fac.count();
  boost::dynamic_bitset<>::size_type iFac = face_fac.find_first();
  std::vector<size_t> face_vect;
  for (size_t iRow = 0; iRow < sizFace; iRow++) {
    Treated[iFac] = 1;
    face_vect.push_back(iFac);
    iFac = face_fac.find_next(iFac);
  }
  int nbExt = EXT.rows();
  Face EXTincd(nbExt);
  auto get_stat = [&](int const &iExt) -> bool {
    for (auto &iFac : face_vect)
      if (extfac_incd[iFac * nbExt + iExt] == 0)
        return false;
    return true;
  };
  for (int iExt = 0; iExt < nbExt; iExt++)
    EXTincd[iExt] = get_stat(iExt);
  int RankFace_fac = RankFace; // expressed from the facets
  int RankFace_ext = nbCol - RankFace_fac;
  int RankFaceTarget_ext = RankFace_ext - 1;
  std::vector<size_t> EXTincd_face_vect;
  EXTincd_face_vect.reserve(nbExt);
  for (int iFac = 0; iFac < nbFac; iFac++)
    if (Treated[iFac] == 0) {
      Face EXTincd_face = EXTincd;
      for (int iExt = 0; iExt < nbExt; iExt++)
        if (extfac_incd[iFac * nbExt + iExt] == 0)
          EXTincd_face[iExt] = 0;
      //
      size_t sizEXT = EXTincd_face.count();
      MyMatrix<T> EXTface(sizEXT, nbCol);
      boost::dynamic_bitset<>::size_type iExt = EXTincd_face.find_first();
      std::vector<size_t> EXTincd_face_vect;
      for (size_t iRow = 0; iRow < sizEXT; iRow++) {
        EXTface.row(iRow) = EXT.row(iExt);
        EXTincd_face_vect.push_back(iExt);
        iExt = EXTincd_face.find_next(iExt);
      }
      Face gList(nbFac);
      if (f_filt(EXTface, RankFaceTarget_ext)) {
        auto test_corr = [&](int const &jFac) -> bool {
          for (auto &iExt : EXTincd_face_vect)
            if (extfac_incd[jFac * nbExt + iExt] == 0)
              return false;
          return true;
        };
        for (int jFac = 0; jFac < nbFac; jFac++)
          if (test_corr(jFac))
            gList[jFac] = 1;
        TheReturn.push_back(gList);
      } else {
        gList[iFac] = 1;
      }
      Face rList = OrbitUnion(StabFace_fac, gList);
      for (int jFac = 0; jFac < nbFac; jFac++)
        if (rList[jFac] == 1)
          Treated[jFac] = 1;
    }
  return TheReturn;
}

template <typename T, typename Tgroup>
vectface SPAN_face_ExtremeRays(Face const &face, Tgroup const &StabFace,
                               int const &RankFace, const Face &extfac_incd,
                               MyMatrix<T> const &FAC, MyMatrix<T> const &EXT) {
  auto f_filt = [&](MyMatrix<T> const &M, int const &RankTarget) -> bool {
    if (M.rows() < RankTarget)
      return false; // The number of rows cannot match the rank, so reject
    return RankMat(M) == RankTarget;
  };
  return SPAN_face_ExtremeRays_F(face, StabFace, RankFace, extfac_incd, FAC,
                                 EXT, f_filt);
}

template <typename T, typename Tgroup>
vectface SPAN_face_ExtremeRaysNonSimplicial(
    Face const &face, Tgroup const &StabFace, int const &RankFace,
    const Face &extfac_incd, MyMatrix<T> const &FAC, MyMatrix<T> const &EXT) {
  auto f_filt = [&](MyMatrix<T> const &M, int const &RankTarget) -> bool {
    //    std::cerr << "|M|=" << M.rows() << " / " << M.cols() << " RankTarget="
    //    << RankTarget << "\n";
    if (M.rows() == RankTarget)
      return false; // If it were to be a face, it would be a simplicial one. So
                    // we remove it from consideration
    if (M.rows() < RankTarget)
      return false; // The number of rows cannot match the rank, so reject
    return RankMat(M) == RankTarget;
  };
  return SPAN_face_ExtremeRays_F(face, StabFace, RankFace, extfac_incd, FAC,
                                 EXT, f_filt);
}

template <typename T, typename Tgroup, typename Fspann, typename Ffinal>
std::vector<vectface>
EnumerationFaces_Fspann_Ffinal(Tgroup const &TheGRP, MyMatrix<T> const &FAC,
                               int LevSearch, Fspann f_spann, Ffinal f_final) {
  std::vector<vectface> RetList;
  int n = TheGRP.n_act();
  vectface ListOrb(n);
  Face eList(n);
  for (int i = 0; i < n; i++)
    eList[i] = 1;
  vectface vvO = DecomposeOrbitPoint(TheGRP, eList);
  for (auto &eOrb : vvO) {
    boost::dynamic_bitset<>::size_type MinVal = eOrb.find_first();
    Face nList(n);
    nList[MinVal] = 1;
    ListOrb.push_back(nList);
  }
  std::cerr << "iLevel=0 |NListOrb|=" << ListOrb.size() << "\n";
  bool test = f_final(0, ListOrb);
  RetList.emplace_back(std::move(ListOrb));
  if (test)
    return RetList;
  for (int iLevel = 1; iLevel <= LevSearch; iLevel++) {
    vectface NListOrb(n);
    for (auto &eOrb : RetList[iLevel - 1]) {
      Tgroup StabFace = TheGRP.Stabilizer_OnSets(eOrb);
      vectface TheSpann = f_spann(eOrb, StabFace, iLevel, FAC, TheGRP);
      for (Face fOrb : TheSpann) {
        Face fOrbCan = TheGRP.CanonicalImage(fOrb);
        NListOrb.push_back(fOrbCan);
      }
    }
    std::cerr << "iLevel=" << iLevel << " |NListOrb|=" << NListOrb.size()
              << "\n";
    bool test = f_final(iLevel, NListOrb);
    RetList.emplace_back(std::move(NListOrb));
    if (test)
      break;
  }
  return RetList;
}

template <typename T>
Face Compute_extfac_incd(const MyMatrix<T> &FAC, const MyMatrix<T> &EXT) {
  int nbFac = FAC.rows();
  int nbExt = EXT.rows();
  int nbCol = EXT.cols();
  Face extfac_incd(nbFac * nbExt);
  for (int iFac = 0; iFac < nbFac; iFac++)
    for (int iExt = 0; iExt < nbExt; iExt++) {
      T sum = 0;
      for (int i = 0; i < nbCol; i++)
        sum += FAC(iFac, i) * EXT(iExt, i);
      if (sum == 0)
        extfac_incd[iFac * nbExt + iExt] = 1;
    }
  return extfac_incd;
}

Face GetFacet_extfac(const Face &extfac_incd, [[maybe_unused]] size_t nbFac,
                     size_t nbExt, size_t iFac) {
  Face f(nbExt);
  for (size_t iExt = 0; iExt < nbExt; iExt++)
    f[iExt] = extfac_incd[iFac * nbExt + iExt];
  return f;
}

Face Compute_facext(const Face &extfac_incd, int nbFac, int nbExt) {
  Face facext_incd(nbFac * nbExt);
  for (int iFac = 0; iFac < nbFac; iFac++)
    for (int iExt = 0; iExt < nbExt; iExt++)
      facext_incd[iExt * nbFac + iFac] = extfac_incd[iFac * nbExt + iExt];
  return facext_incd;
}

Face Compute_faceEXT_from_faceFAC(const Face &extfac_incd, int nbFac, int nbExt,
                                  const Face &face_fac) {
  std::vector<int> eV;
  for (int iFac = 0; iFac < nbFac; iFac++)
    if (face_fac[iFac] == 1)
      eV.push_back(iFac);
  auto is_in = [&](int iExt) -> bool {
    for (auto &iFac : eV)
      if (extfac_incd[iFac * nbExt + iExt] == 0)
        return false;
    return true;
  };
  Face f_ext(nbExt);
  for (int iExt = 0; iExt < nbExt; iExt++)
    f_ext[iExt] = is_in(iExt);
  return f_ext;
}

Face Compute_faceFAC_from_faceEXT(const Face &extfac_incd, int nbFac, int nbExt,
                                  const Face &face_ext) {
  std::vector<int> eV;
  for (int iExt = 0; iExt < nbExt; iExt++)
    if (face_ext[iExt] == 1)
      eV.push_back(iExt);
  auto is_in = [&](int iFac) -> bool {
    for (auto &iExt : eV)
      if (extfac_incd[iFac * nbExt + iExt] == 0)
        return false;
    return true;
  };
  Face f_fac(nbFac);
  for (int iFac = 0; iFac < nbFac; iFac++)
    f_fac[iFac] = is_in(iFac);
  return f_fac;
}

template <typename T, typename Tgroup, typename Final>
std::vector<vectface>
EnumerationFaces_Ffinal(Tgroup const &TheGRP, MyMatrix<T> const &FAC,
                        MyMatrix<T> const &EXT, int LevSearch,
                        std::string const &method_spann, Final f_final) {
  if (method_spann == "LinearProgramming") {
    auto f_spann = [&](Face const &face, Tgroup const &StabFace,
                       [[maybe_unused]] int const &RankFace,
                       MyMatrix<T> const &FAC,
                       Tgroup const &FullGRP) -> vectface {
      return SPAN_face_LinearProgramming(face, StabFace, FAC, FullGRP);
    };
    return EnumerationFaces_Fspann_Ffinal<T, Tgroup, decltype(f_spann),
                                          decltype(f_final)>(
        TheGRP, FAC, LevSearch, f_spann, f_final);
  }
  if (method_spann == "ExtremeRays") {
    Face extfac_incd = Compute_extfac_incd(FAC, EXT);
    auto f_spann = [&](Face const &face, Tgroup const &StabFace,
                       int const &RankFace, MyMatrix<T> const &FAC,
                       [[maybe_unused]] Tgroup const &FullGRP) -> vectface {
      return SPAN_face_ExtremeRays(face, StabFace, RankFace, extfac_incd, FAC,
                                   EXT);
    };
    return EnumerationFaces_Fspann_Ffinal<T, Tgroup, decltype(f_spann),
                                          decltype(f_final)>(
        TheGRP, FAC, LevSearch, f_spann, f_final);
  }
  if (method_spann == "ExtremeRaysNonSimplicial") {
    int nbFac = FAC.rows();
    int nbExt = EXT.rows();
    int nbCol = EXT.cols();
    Face extfac_incd(nbFac * nbExt);
    for (int iFac = 0; iFac < nbFac; iFac++)
      for (int iExt = 0; iExt < nbExt; iExt++) {
        T sum = 0;
        for (int i = 0; i < nbCol; i++)
          sum += FAC(iFac, i) * EXT(iExt, i);
        if (sum == 0)
          extfac_incd[iFac * nbExt + iExt] = 1;
      }
    auto f_spann = [&](Face const &face, Tgroup const &StabFace,
                       int const &RankFace, MyMatrix<T> const &FAC,
                       [[maybe_unused]] Tgroup const &FullGRP) -> vectface {
      return SPAN_face_ExtremeRaysNonSimplicial(face, StabFace, RankFace,
                                                extfac_incd, FAC, EXT);
    };
    return EnumerationFaces_Fspann_Ffinal<T, Tgroup, decltype(f_spann),
                                          decltype(f_final)>(
        TheGRP, FAC, LevSearch, f_spann, f_final);
  }
  std::cerr << "We failed to find a matching method_spann\n";
  throw TerminalException{1};
}

template <typename T, typename Tgroup>
std::vector<vectface> EnumerationFaces(Tgroup const &TheGRP,
                                       MyMatrix<T> const &FAC,
                                       MyMatrix<T> const &EXT, int LevSearch,
                                       std::string const &method_spann,
                                       std::string const &method_final) {
  if (method_final == "all") {
    auto f_final = [&]([[maybe_unused]] int const &level,
                       [[maybe_unused]] vectface const &RetList) -> bool {
      return false;
    };
    return EnumerationFaces_Ffinal(TheGRP, FAC, EXT, LevSearch, method_spann,
                                   f_final);
  }
  if (method_final == "stop_nonsimplicial") {
    auto f_final = [&](int const &level, vectface const &RetList) -> bool {
      size_t auth_len = level + 1;
      for (auto &eFace : RetList) {
        if (eFace.count() != auth_len)
          return true;
      }
      return false;
    };
    return EnumerationFaces_Ffinal(TheGRP, FAC, EXT, LevSearch, method_spann,
                                   f_final);
  }
  std::cerr << "We failed to find a matching method_final\n";
  throw TerminalException{1};
}

// We test if eSet is included in a proper face of the polytope
template <typename T>
bool TestInclusionProperFace(std::vector<int> const &eSet,
                             MyMatrix<T> const &FAC) {
  int nbCol = FAC.cols();
  int nbRow = FAC.rows();
  std::vector<int> eVectCand(nbRow, 0);
  for (auto eVal : eSet)
    eVectCand[eVal] = 1;
  MyMatrix<T> eMatId = IdentityMat<T>(nbCol);
  while (true) {
    int len = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      if (eVectCand[iRow] == 1)
        len++;
    MyMatrix<T> TestMat(len, nbCol);
    int jRow = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      if (eVectCand[iRow] == 1) {
        for (int iCol = 0; iCol < nbCol; iCol++) {
          T eVal = FAC(iRow, iCol);
          TestMat(jRow, iCol) = eVal;
        }
        jRow++;
      }
    SelectionRowCol<T> eSelect = TMat_SelectRowCol(TestMat);
    MyMatrix<T> NSP = eSelect.NSP;
    int nbEqua = NSP.rows();
    std::vector<int> eCand, eCandCompl;
    for (int kRow = 0; kRow < nbRow; kRow++) {
      int test = 1;
      for (int iEqua = 0; iEqua < nbEqua; iEqua++)
        if (test == 1) {
          T eSum = 0;
          for (int iCol = 0; iCol < nbCol; iCol++)
            eSum += FAC(kRow, iCol) * NSP(iEqua, iCol);
          if (eSum != 0)
            test = 0;
        }
      if (test == 1)
        eCand.push_back(kRow);
      else
        eCandCompl.push_back(kRow);
    }
    int nbEltCompl = eCandCompl.size();
    if (nbEltCompl == 0)
      return false;
    int nbElt = eCand.size();
    MyMatrix<T> ListVectSpann(nbElt, nbCol);
    for (int iElt = 0; iElt < nbElt; iElt++) {
      int iRow = eCand[iElt];
      for (int iCol = 0; iCol < nbCol; iCol++) {
        T eVal = FAC(iRow, iCol);
        ListVectSpann(iElt, iCol) = eVal;
      }
    }
    MyMatrix<T> BasisSpann = RowReduction(ListVectSpann);
    int nbRowSpann = BasisSpann->nbRow;
    int LPdim = nbCol - nbRowSpann;
    MyMatrix<T> PreTheTot = Concatenate(BasisSpann, eMatId);
    MyMatrix<T> TheTot = RowReduction(PreTheTot);
    // the complement
    MyMatrix<T> PreListVectors(nbEltCompl, nbCol);
    for (int iElt = 0; iElt < nbEltCompl; iElt++) {
      int iRow = eCandCompl[iElt];
      for (int iCol = 0; iCol < nbCol; iCol++) {
        T eVal = FAC(iRow, iCol);
        PreListVectors(iElt, iCol) = eVal;
      }
    }
    MyMatrix<T> TheTotInv = Inverse(TheTot);
    MyMatrix<T> PreListVectorsB = PreListVectors * TheTotInv;
    MyMatrix<T> ListVectors(nbEltCompl, LPdim);
    for (int iElt = 0; iElt < nbEltCompl; iElt++)
      for (int iCol = 0; iCol < LPdim; iCol++) {
        int jCol = iCol + nbRowSpann;
        T eVal = PreListVectorsB(iElt, jCol);
        ListVectors(iElt, iCol) = eVal;
      }
    PosRelRes<T> eResult = SearchPositiveRelationSimple(ListVectors);
    if (eResult.eTestExist) {
      for (int iElt = 0; iElt < nbEltCompl; iElt++) {
        T eVal = eResult.TheRelat(iElt);
        if (eVal > 0)
          eVectCand[eCandCompl[iElt]] = 1;
      }
    } else {
      return true;
    }
  }
}

template <typename T>
bool TestPositiveRelationSimple(MyMatrix<T> const &ListVect) {
  PosRelRes<T> eResult = SearchPositiveRelationSimple(ListVect);
  return eResult.eTestExist;
}

template <typename Tgroup>
std::vector<std::vector<int>> GetMinimalReprVertices(Tgroup const &TheGRP) {
  Face eList;
  int n = TheGRP.n_act();
  for (int i = 0; i < n; i++)
    eList[i] = 1;
  vectface vvO = DecomposeOrbitPoint(TheGRP, eList);
  std::vector<std::vector<int>> RetList;
  for (auto eOrb : vvO) {
    boost::dynamic_bitset<>::size_type MinVal = eOrb.find_first();
    std::vector<int> nList{int(MinVal)};
    RetList.push_back(nList);
  }
  return RetList;
}

void PrintListOrb_GAP(std::ostream &os,
                      std::vector<std::vector<int>> const &ListOrb) {
  int nbOrb, iOrb, len, i;
  std::vector<int> eOrb;
  os << "return ";
  os << "rec(ListRepresentent:=[";
  nbOrb = ListOrb.size();
  for (iOrb = 0; iOrb < nbOrb; iOrb++) {
    if (iOrb > 0)
      os << ",";
    os << "[";
    eOrb = ListOrb[iOrb];
    len = eOrb.size();
    for (i = 0; i < len; i++) {
      if (i > 0)
        os << ",";
      os << eOrb[i] + 1;
    }
    os << "]";
  }
  os << "]);\n";
}

void PrintListListOrb_IntGAP(std::ostream &os,
                             std::vector<vectface> const &ListListOrb) {
  int nbLev = ListListOrb.size();
  os << "return [";
  for (int iLev = 0; iLev < nbLev; iLev++) {
    if (iLev > 0)
      os << ",\n";
    os << "rec(ListRepresentent:=[";
    int nbOrb = ListListOrb[iLev].size();
    for (int iOrb = 0; iOrb < nbOrb; iOrb++) {
      if (iOrb > 0)
        os << ",";
      os << "[";
      Face eOrb = ListListOrb[iLev][iOrb];
      size_t len = eOrb.count();
      boost::dynamic_bitset<>::size_type aRow = eOrb.find_first();
      for (size_t i = 0; i < len; i++) {
        if (i > 0)
          os << ",";
        os << aRow + 1;
        aRow = eOrb.find_next(aRow);
      }
      os << "]";
    }
    os << "])";
  }
  os << "];\n";
}

void OutputFaces(const std::vector<vectface> &TheOutput,
                 const std::string &OUTfile, const std::string &OutFormat) {
  if (OutFormat == "GAP") {
    std::ofstream os(OUTfile);
    os << "return ";
    os << "[";
    int len = TheOutput.size();
    for (int i = 0; i < len; i++) {
      if (i > 0)
        os << ",\n";
      VectVectInt_Gap_Print(os, TheOutput[i]);
    }
    os << "]";
    os << ";\n";
    return;
  }
  std::cerr << "No option has been chosen\n";
  throw TerminalException{1};
}

FullNamelist NAMELIST_GetStandard_FaceLattice() {
  std::map<std::string, SingleBlock> ListBlock;
  // DATA
  std::map<std::string, std::string> ListStringValues1;
  std::map<std::string, bool> ListBoolValues1;
  std::map<std::string, int> ListIntValues1;
  ListStringValues1["EXTfile"] = "unset.ext";
  ListStringValues1["FACfile"] = "unset.ext";
  ListStringValues1["GRPfile"] = "unset.grp";
  ListStringValues1["OUTfile"] = "unset.out";
  ListStringValues1["OutFormat"] = "GAP";
  ListStringValues1["method_spann"] = "unset.out";
  ListStringValues1["method_final"] = "unset.out";
  ListIntValues1["LevSearch"] = -1;
  SingleBlock BlockDATA;
  BlockDATA.ListStringValues = ListStringValues1;
  BlockDATA.ListBoolValues = ListBoolValues1;
  BlockDATA.ListIntValues = ListIntValues1;
  ListBlock["PROC"] = BlockDATA;
  // Merging all data
  return {std::move(ListBlock), "undefined"};
}

template <typename T, typename Tgroup>
void MainFunctionFaceLattice(FullNamelist const &eFull) {
  using Telt = typename Tgroup::Telt;
  using Tidx = typename Telt::Tidx;
  std::cerr << "Reading PROC\n";
  SingleBlock BlockPROC = eFull.ListBlock.at("PROC");
  std::string FACfile = BlockPROC.ListStringValues.at("FACfile");
  IsExistingFileDie(FACfile);
  std::cerr << "FACfile=" << FACfile << "\n";
  std::ifstream FACfs(FACfile);
  MyMatrix<T> FAC = ReadMatrix<T>(FACfs);
  if (size_t(FAC.rows()) > size_t(std::numeric_limits<Tidx>::max())) {
    std::cerr << "We have |FAC|=" << FAC.rows() << "\n";
    std::cerr << "But <Tidx>::max()="
              << size_t(std::numeric_limits<Tidx>::max()) << "\n";
    throw TerminalException{1};
  }
  if (RankMat(FAC) != FAC.cols()) {
    std::cerr << "The matrix EXT should be of full rank\n";
    throw TerminalException{1};
  }
  //
  std::string method_spann = BlockPROC.ListStringValues.at("method_spann");
  std::cerr << "method_spann=" << method_spann << "\n";
  std::string method_final = BlockPROC.ListStringValues.at("method_final");
  std::cerr << "method_final=" << method_final << "\n";
  //
  MyMatrix<T> EXT;
  if (method_spann == "ExtremeRays" ||
      method_spann == "ExtremeRaysNonSimplicial") {
    std::string EXTfile = BlockPROC.ListStringValues.at("EXTfile");
    IsExistingFileDie(EXTfile);
    std::cerr << "EXTfile=" << EXTfile << "\n";
    std::ifstream EXTfs(EXTfile);
    EXT = ReadMatrix<T>(EXTfs);
    if (FAC.cols() != EXT.cols()) {
      std::cerr << "The dimension of EXT and FAC should be the same\n";
      throw TerminalException{1};
    }
  }
  //
  std::string GRPfile = BlockPROC.ListStringValues.at("GRPfile");
  IsExistingFileDie(GRPfile);
  std::cerr << "GRPfile=" << GRPfile << "\n";
  std::ifstream GRPfs(GRPfile);
  Tgroup GRP = ReadGroup<Tgroup>(GRPfs);
  //
  int LevSearch = BlockPROC.ListIntValues.at("LevSearch");
  std::cerr << "LevSearch=" << LevSearch << "\n";
  if (LevSearch == -1) {
    int nbCol = FAC.cols();
    LevSearch = nbCol - 2;
  }
  //
  std::string OUTfile = BlockPROC.ListStringValues.at("OUTfile");
  std::string OutFormat = BlockPROC.ListStringValues.at("OutFormat");
  std::cerr << "OUTfile=" << OUTfile << " OutFormat=" << OutFormat << "\n";
  //
  std::vector<vectface> TheOutput =
      EnumerationFaces(GRP, FAC, EXT, LevSearch, method_spann, method_final);
  //
  OutputFaces(TheOutput, OUTfile, OutFormat);
}

#endif //  SRC_POLY_POLY_KSKELETTON_H_
