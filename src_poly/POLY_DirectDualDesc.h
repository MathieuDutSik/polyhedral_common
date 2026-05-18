// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DIRECTDUALDESC_H_
#define SRC_POLY_POLY_DIRECTDUALDESC_H_

// clang-format off
#include "Basic_string.h"
#include "POLY_c_cddlib.h"
#include "POLY_cddlib.h"
#include "POLY_lrslib.h"
#include "MAT_MatrixInt.h"
#ifndef WASM_PLATFORM
#include "POLY_DirectDualDesc_External.h"
#endif
#include "POLY_DualDescription_PrimalDual.h"
#include "SmallPolytopes.h"
#include <algorithm>
#include <utility>
#include <string>
#include <vector>
// clang-format on

#ifdef TIMINGS
#define TIMINGS_DUAL_DESC
#endif

#ifdef DEBUG
#define DEBUG_DUAL_DESC
#endif

#ifdef DISABLE_DEBUG_DUAL_DESC
#undef DEBUG_DUAL_DESC
#endif

#ifdef KEY_VALUE
#define KEY_VALUE_DUAL_DESC
#endif

template <typename T> std::vector<size_t> Convert_T_To_Set(T const &val) {
  size_t pos = 0;
  std::vector<size_t> V;
  T valWork = val;
  T two = 2;
  while (true) {
    if (valWork == 0)
      break;
    T res = ResInt(valWork, two);
    if (res == 1) {
      V.push_back(pos);
    }
    valWork = (valWork - res) / two;
    pos++;
  }
  return V;
}

template <typename T> T Convert_Set_To_T(std::vector<size_t> const &V) {
  T ThePow = 1;
  size_t pos = 0;
  T retval = 0;
  for (auto &eV : V) {
    while (true) {
      if (pos == eV)
        break;
      ThePow *= 2;
      pos++;
    }
    retval += ThePow;
  }
  return retval;
}

template <typename T> bool is_method_supported(std::string const &prog) {
  if (prog == "cdd_cbased") {
#ifdef USE_CDDLIB
    return true;
#else
    return false;
#endif
  }
  //
  if constexpr (is_ring_field<T>::value) {
    if (prog == "cdd")
      return true;
    // If it is a field, then it makes sense to look at the internal ring
    if (prog == "lrs_ring")
      return true;
  }
  //
  if (prog == "small_polytopes")
    return true;
  // It applies to the field case or ring
  if (prog == "lrs")
    return true;
  // It applies to the field case or ring
  if (prog == "pd_lrs")
    return true;
  //
#ifndef WASM_PLATFORM
  if constexpr (is_implementation_of_Q<T>::value) {
    if (prog == "glrs")
      return true;
    //
    if (prog == "ppl_ext")
      return true;
    //
    if (prog == "cdd_ext")
      return true;
    //
    if (prog == "normaliz")
      return true;
  }
#endif
  return false;
}

[[noreturn]] void
terminate_direct_dual_desc(std::string const &ansProg,
                           std::vector<std::string> const &ListProg) {
  std::cerr << "DDD: ERROR: No right program found with ansProg=" << ansProg
            << "\n";
  std::cerr << "DDD: List of authorized programs :";
  bool IsFirst = true;
  for (auto &eP : ListProg) {
    if (!IsFirst)
      std::cerr << " ,";
    IsFirst = false;
    std::cerr << " " << eP;
  }
  std::cerr << "\n";
  throw TerminalException{1};
}

template <typename T>
vectface DirectFacetComputationIncidence(MyMatrix<T> const &EXT,
                                         std::string const &ansProg,
                                         std::ostream &os) {
#ifdef DEBUG_DUAL_DESC
  os << "DDD: DirectFacetComputationIncidence, ansProg=" << ansProg << "\n";
#endif
  std::string eProg;
  std::vector<std::string> ListProg;
  //
  eProg = "cdd_cbased";
  ListProg.push_back(eProg);
  if (ansProg == eProg) {
#ifdef USE_CDDLIB
    return cbased_cdd::DualDescription_incd(EXT);
#else
    std::cerr << "DDD: The code has been compiled without the CDDLIB library\n";
    throw TerminalException{1};
#endif
  }
  //
  if constexpr (is_ring_field<T>::value) {
    // CDD certainly requires the ring to be a field
    eProg = "cdd";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return cdd::DualDescription_incd(EXT, os);
    // If it is a field, then it makes sense to look at the internal ring
    eProg = "lrs_ring";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return lrs::DualDescription_incd_reduction(EXT);
  }
  // Small polytopes have special solutions
  eProg = "small_polytopes";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return SmallPolytope_Incidence(EXT, os);
  // It applies to the field case or ring
  eProg = "lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescription_incd(EXT);
  // It applies to the field case or ring
  eProg = "pd_lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return POLY_DualDescription_PrimalDualIncidence(EXT, os);
  //
  //
  // The external programs are available only for rationl types
  //
#ifndef WASM_PLATFORM
  if constexpr (is_implementation_of_Q<T>::value) {
    eProg = "glrs";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIncidence(EXT, "glrs", os);
    //
    eProg = "ppl_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIncidence(EXT, "ppl_lcdd", os);
    //
    eProg = "cdd_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIncidence(EXT, "lcdd_gmp", os);
    //
    eProg = "normaliz";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIncidence(EXT, "normaliz", os);
  }
#endif
  //
  std::cerr << "DDD: ERROR in DirectFacetComputationIncidence\n";
  terminate_direct_dual_desc(ansProg, ListProg);
}

template <typename T>
MyMatrix<T> DirectFacetComputationInequalities(MyMatrix<T> const &EXT,
                                               std::string const &ansProg,
                                               std::ostream &os) {
#ifdef DEBUG_DUAL_DESC
  os << "DDD: DirectFacetComputationInequalities, ansProg=" << ansProg << "\n";
#endif
  std::string eProg;
  std::vector<std::string> ListProg;
  //
  if constexpr (is_ring_field<T>::value) {
    // CDD certainly requires a field for working out.
    eProg = "cdd";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return cdd::DualDescription(EXT, os);
    // For lrs_ring, we certainly need to have a field for T
    eProg = "lrs_ring";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return lrs::DualDescription_reduction(EXT);
  }
  // Small polytopes have special solutions
  eProg = "small_polytopes";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return SmallPolytope_Ineq(EXT, os);
  // lrs does not use divisions, so work even if not field.
  eProg = "lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescription(EXT);
  // It applies to the field case or ring
  eProg = "pd_lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return POLY_DualDescription_PrimalDualInequalities(EXT, os);
  //
  // The external programs are available only for rationl types
  //
#ifndef WASM_PLATFORM
  if constexpr (is_implementation_of_Q<T>::value) {
    eProg = "glrs";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIneq(EXT, "glrs", os);
    //
    eProg = "ppl_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIneq(EXT, "ppl_lcdd", os);
    //
    eProg = "cdd_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIneq(EXT, "lcdd_gmp", os);
    //
    eProg = "normaliz";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramIneq(EXT, "normaliz", os);
  }
#endif
  //
  std::cerr << "DDD: ERROR in DirectFacetComputationInequalities\n";
  terminate_direct_dual_desc(ansProg, ListProg);
}

template <typename T, typename Fprocess>
void DirectFacetComputationFaceIneq(MyMatrix<T> const &EXT,
                                    std::string const &ansProg,
                                    Fprocess f_process, std::ostream &os) {
#ifdef DEBUG_DUAL_DESC
  os << "DDD: DirectFacetComputationFaceIneq, ansProg=" << ansProg << "\n";
#endif
  std::string eProg;
  std::vector<std::string> ListProg;
  //
  if constexpr (is_ring_field<T>::value) {
    // CDD requires for T to be a field
    eProg = "cdd";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return cdd::DualDescriptionFaceIneq(EXT, f_process, os);
    // For lrs_ring that computes in a subring, we need T to be a field.
    eProg = "lrs_ring";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return lrs::DualDescriptionFaceIneq_reduction(EXT, f_process);
  }
  // We need to make it work also for ase without reduction if that makes sense
  // which is not sure at all.
  // Small polytopes can have special solutions
  eProg = "small_polytopes";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return SmallPolytope_FaceIneq(EXT, f_process, os);
  // T can be a field or a ring here
  eProg = "lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return lrs::DualDescriptionFaceIneq(EXT, f_process);
  // It applies to the field case or ring
  eProg = "pd_lrs";
  ListProg.push_back(eProg);
  if (ansProg == eProg)
    return POLY_DualDescription_PrimalDualFaceIneq(EXT, f_process, os);
  //
  // The external programs are available only for rationl types
  //
#ifndef WASM_PLATFORM
  if constexpr (is_implementation_of_Q<T>::value) {
    eProg = "glrs";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramFaceIneq(EXT, "glrs", f_process, os);
    //
    eProg = "ppl_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramFaceIneq(EXT, "ppl_lcdd", f_process, os);
    //
    eProg = "cdd_ext";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramFaceIneq(EXT, "lcdd_gmp", f_process, os);
    //
    eProg = "normaliz";
    ListProg.push_back(eProg);
    if (ansProg == eProg)
      return DualDescExternalProgramFaceIneq(EXT, "normaliz", f_process, os);
  }
#endif
  //
  std::cerr << "DDD: ERROR in DirectFacetComputationFaceIneq\n";
  terminate_direct_dual_desc(ansProg, ListProg);
}

template <typename T, typename Tgroup>
vectface DirectFacetOrbitComputation(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                     std::string const &ansProg,
                                     std::ostream &os) {
#ifdef TIMINGS_DUAL_DESC
  MicrosecondTime time;
#endif
#ifdef KEY_VALUE_DUAL_DESC
  MicrosecondTime time_total;
#endif
  vectface ListIncd = DirectFacetComputationIncidence(EXT, ansProg, os);
#ifdef DEBUG_DUAL_DESC
  os << "DDD: |ListIncd|=" << ListIncd.size() << "\n";
#endif
#ifdef TIMINGS_DUAL_DESC
  os << "|DDD: DualDescription|=" << time << "\n";
#endif
  if (ListIncd.empty()) {
    std::cerr << "We found ListIncd to be empty. A clear error\n";
    throw TerminalException{1};
  }
  if (GRP.size() == 1) {
    return ListIncd;
  }
  vectface TheOutput = OrbitSplittingSet(ListIncd, GRP);
#ifdef KEY_VALUE_DUAL_DESC
  os << "DDD: KEY=(DirectFacetOrbitComputation_" << EXT.rows() << "_"
     << EXT.cols() << "_" << GRP.size() << "_" << ansProg << "_"
     << ListIncd.size() << "_" << TheOutput.size() << ") VALUE=(" << time_total
     << ")\n";
#endif
  return TheOutput;
}

template <typename T, typename Tgroup>
std::vector<std::pair<Face, MyVector<T>>>
DirectFacetIneqOrbitComputation(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                std::string const &ansProg, std::ostream &os) {
#ifdef TIMINGS_DUAL_DESC
  MicrosecondTime time;
#endif
#ifdef KEY_VALUE_DUAL_DESC
  MicrosecondTime time_total;
#endif
  std::vector<std::pair<Face, MyVector<T>>> ListReturn;
  auto f_process = [&](std::pair<Face, MyVector<T>> const &pair_face) -> void {
    ListReturn.push_back(pair_face);
  };
  DirectFacetComputationFaceIneq(EXT, ansProg, f_process, os);
#ifdef TIMINGS_DUAL_DESC
  os << "|DDD: DualDescription|=" << time << "\n";
#endif
#ifdef DEBUG_DUAL_DESC
  os << "DDD: Found  |ListIncd|=" << ListReturn.size() << "\n";
#endif
  if (ListReturn.empty()) {
    std::cerr << "We found ListIncd to be empty. A clear error\n";
    throw TerminalException{1};
  }
  if (GRP.size() == 1) {
    return ListReturn;
  }
  std::vector<std::pair<Face, MyVector<T>>> TheOutput =
      OrbitSplittingMap(ListReturn, GRP);
#ifdef TIMINGS_DUAL_DESC
  os << "|DDD: OrbitSplittingMap|=" << time << "\n";
#endif
#ifdef KEY_VALUE_DUAL_DESC
  os << "DDD: KEY=(DirectFacetIneqOrbitComputation_" << EXT.rows() << "_"
     << EXT.cols() << "_" << GRP.size() << "_" << ansProg << "_"
     << ListReturn.size() << "_" << TheOutput.size() << ") VALUE=("
     << time_total << ")\n";
#endif
  return TheOutput;
}

/*
  Returns the method that we can use for the dual description.
  We need somewhat better method for choosing the heuristics.
 */
template <typename T>
std::string get_dual_desc_method([[maybe_unused]] MyMatrix<T> const &EXT, [[maybe_unused]] std::ostream &os) {
  std::string ansProg = "lrs";
  return ansProg;
}


template <typename T>
MyMatrix<T> DirectDualDescription_mat(MyMatrix<T> const &EXT, std::ostream &os) {
  std::string ansProg = get_dual_desc_method(EXT, os);
  return DirectFacetComputationInequalities(EXT, ansProg, os);
}

template <typename T>
vectface DirectDualDescription_vf(MyMatrix<T> const &EXT, std::ostream &os) {
  std::string ansProg = get_dual_desc_method(EXT, os);
  return DirectFacetComputationIncidence(EXT, ansProg, os);
}

template <typename T, typename Fprocess>
void DirectDualDescription_f(MyMatrix<T> const &EXT,
                             Fprocess f_process, std::ostream &os) {
  std::string ansProg = get_dual_desc_method(EXT, os);
  return DirectFacetComputationFaceIneq(EXT, ansProg, f_process, os);
}


// clang-format off
#endif  // SRC_POLY_POLY_DIRECTDUALDESC_H_
// clang-format on
