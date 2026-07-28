// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_DIRECTDUALDESC_H_
#define SRC_POLY_POLY_DIRECTDUALDESC_H_

// clang-format off
#include "Basic_string.h"
#include "POLY_double_description.h"
#include "POLY_lrslib.h"
#include "MAT_MatrixInt.h"
#ifndef WASM_PLATFORM
#include "POLY_DirectDualDesc_External.h"
#endif
#include "POLY_DualDescription_PrimalDual.h"
#include "POLY_DualDescription_BeneathBeyond.h"
#include "SmallPolytopes.h"
#include <algorithm>
#include <optional>
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

// The dual description program to use for a direct (non-recursive)
// computation. This is the type-safe replacement for the "ansProg" strings.
// The string form is still available (parsed from the heuristics / command
// line) via dual_desc_program_from_string.
enum class DualDescProgram {
  cdd,
  small_polytopes,
  lrs,
  pd_lrs,
  beneath_beyond,
  glrs,
  ppl_ext,
  cdd_ext,
  normaliz,
};

// The string encoding of the program, used for logging and error messages and
// for writing the value back into the string-based heuristic machinery.
namespace std {
inline std::string to_string(DualDescProgram prog) {
  switch (prog) {
  case DualDescProgram::cdd:
    return "cdd";
  case DualDescProgram::small_polytopes:
    return "small_polytopes";
  case DualDescProgram::lrs:
    return "lrs";
  case DualDescProgram::pd_lrs:
    return "pd_lrs";
  case DualDescProgram::beneath_beyond:
    return "beneath_beyond";
  case DualDescProgram::glrs:
    return "glrs";
  case DualDescProgram::ppl_ext:
    return "ppl_ext";
  case DualDescProgram::cdd_ext:
    return "cdd_ext";
  case DualDescProgram::normaliz:
    return "normaliz";
  }
  return "unknown";
}
}  // namespace std

// Parses the string form of the program. Returns nullopt if the string is not
// a known dual description program (e.g. "fullrankfacetset" handled elsewhere).
inline std::optional<DualDescProgram>
dual_desc_program_from_string_opt(std::string const &prog) {
  if (prog == "cdd")
    return DualDescProgram::cdd;
  if (prog == "small_polytopes")
    return DualDescProgram::small_polytopes;
  if (prog == "lrs")
    return DualDescProgram::lrs;
  if (prog == "pd_lrs")
    return DualDescProgram::pd_lrs;
  if (prog == "beneath_beyond" || prog == "bb")
    return DualDescProgram::beneath_beyond;
  if (prog == "glrs")
    return DualDescProgram::glrs;
  if (prog == "ppl_ext")
    return DualDescProgram::ppl_ext;
  if (prog == "cdd_ext")
    return DualDescProgram::cdd_ext;
  if (prog == "normaliz")
    return DualDescProgram::normaliz;
  return {};
}

// Same as above but errors out on an unknown program.
inline DualDescProgram dual_desc_program_from_string(std::string const &prog) {
  std::optional<DualDescProgram> opt = dual_desc_program_from_string_opt(prog);
  if (!opt) {
    std::cerr << "DDD: ERROR: unknown dual description program prog=" << prog
              << "\n";
    throw TerminalException{1};
  }
  return *opt;
}

template <typename T> bool is_method_supported(DualDescProgram prog) {
  switch (prog) {
  case DualDescProgram::cdd:
    // Served by the double description method, which runs division-free.
  case DualDescProgram::small_polytopes:
  case DualDescProgram::lrs:
  case DualDescProgram::pd_lrs:
    // Applies to the field or ring case.
    return true;
  case DualDescProgram::beneath_beyond:
    // Beneath-and-beyond computes facet normals through a nullspace, so it
    // requires T to be a field.
    return is_ring_field<T>::value;
  case DualDescProgram::glrs:
  case DualDescProgram::ppl_ext:
  case DualDescProgram::cdd_ext:
  case DualDescProgram::normaliz:
    // The external programs are available only for rational types.
#ifndef WASM_PLATFORM
    return is_implementation_of_Q<T>::value;
#else
    return false;
#endif
  }
  return false;
}

template <typename T> bool is_method_supported(std::string const &prog) {
  std::optional<DualDescProgram> opt = dual_desc_program_from_string_opt(prog);
  return opt.has_value() && is_method_supported<T>(*opt);
}

[[noreturn]] inline void
terminate_direct_dual_desc(std::string const &context, DualDescProgram prog) {
  std::cerr << "DDD: ERROR in " << context
            << ": no available handler for program "
            << std::to_string(prog) << "\n";
  throw TerminalException{1};
}

template <typename T>
vectface DirectFacetComputationIncidence(MyMatrix<T> const &EXT,
                                         DualDescProgram prog,
                                         std::ostream &os) {
#ifdef DEBUG_DUAL_DESC
  os << "DDD: DirectFacetComputationIncidence, prog="
     << std::to_string(prog) << "\n";
#endif
  switch (prog) {
  case DualDescProgram::cdd:
    // The double description method, applicable to the field or ring case
    return double_desc::DualDescription_incd(EXT, os);
  case DualDescProgram::small_polytopes:
    // Small polytopes have special solutions
    return SmallPolytope_Incidence(EXT, os);
  case DualDescProgram::lrs:
    return lrs::DualDescription_incd(EXT);
  case DualDescProgram::pd_lrs:
    // It applies to the field case or ring
    return POLY_DualDescription_PrimalDualIncidence(EXT, os);
  case DualDescProgram::beneath_beyond:
    // Native beneath-and-beyond, full-dimensional pointed cone (field only)
    if constexpr (is_ring_field<T>::value)
      return POLY_DualDescription_BeneathBeyondIncidence(EXT, os);
    break;
  case DualDescProgram::glrs:
  case DualDescProgram::ppl_ext:
  case DualDescProgram::cdd_ext:
  case DualDescProgram::normaliz:
    // The external programs are available only for rational types
#ifndef WASM_PLATFORM
    if constexpr (is_implementation_of_Q<T>::value) {
      if (prog == DualDescProgram::glrs)
        return DualDescExternalProgramIncidence(EXT, "glrs", os);
      if (prog == DualDescProgram::ppl_ext)
        return DualDescExternalProgramIncidence(EXT, "ppl_lcdd", os);
      if (prog == DualDescProgram::cdd_ext)
        return DualDescExternalProgramIncidence(EXT, "lcdd_gmp", os);
      if (prog == DualDescProgram::normaliz)
        return DualDescExternalProgramIncidence(EXT, "normaliz", os);
    }
#endif
    break;
  }
  terminate_direct_dual_desc("DirectFacetComputationIncidence", prog);
}

template <typename T>
vectface DirectFacetComputationIncidence(MyMatrix<T> const &EXT,
                                         std::string const &ansProg,
                                         std::ostream &os) {
  return DirectFacetComputationIncidence(
      EXT, dual_desc_program_from_string(ansProg), os);
}

template <typename T>
MyMatrix<T> DirectFacetComputationInequalities(MyMatrix<T> const &EXT,
                                               DualDescProgram prog,
                                               std::ostream &os) {
#ifdef DEBUG_DUAL_DESC
  os << "DDD: DirectFacetComputationInequalities, prog="
     << std::to_string(prog) << "\n";
#endif
  switch (prog) {
  case DualDescProgram::cdd:
    // The double description method, applicable to the field or ring case
    return double_desc::DualDescription(EXT, os);
  case DualDescProgram::small_polytopes:
    // Small polytopes have special solutions
    return SmallPolytope_Ineq(EXT, os);
  case DualDescProgram::lrs:
    return lrs::DualDescription(EXT);
  case DualDescProgram::pd_lrs:
    // It applies to the field case or ring
    return POLY_DualDescription_PrimalDualInequalities(EXT, os);
  case DualDescProgram::beneath_beyond:
    // Native beneath-and-beyond, full-dimensional pointed cone (field only)
    if constexpr (is_ring_field<T>::value)
      return POLY_DualDescription_BeneathBeyondInequalities(EXT, os);
    break;
  case DualDescProgram::glrs:
  case DualDescProgram::ppl_ext:
  case DualDescProgram::cdd_ext:
  case DualDescProgram::normaliz:
    // The external programs are available only for rational types
#ifndef WASM_PLATFORM
    if constexpr (is_implementation_of_Q<T>::value) {
      if (prog == DualDescProgram::glrs)
        return DualDescExternalProgramIneq(EXT, "glrs", os);
      if (prog == DualDescProgram::ppl_ext)
        return DualDescExternalProgramIneq(EXT, "ppl_lcdd", os);
      if (prog == DualDescProgram::cdd_ext)
        return DualDescExternalProgramIneq(EXT, "lcdd_gmp", os);
      if (prog == DualDescProgram::normaliz)
        return DualDescExternalProgramIneq(EXT, "normaliz", os);
    }
#endif
    break;
  }
  terminate_direct_dual_desc("DirectFacetComputationInequalities", prog);
}

template <typename T>
MyMatrix<T> DirectFacetComputationInequalities(MyMatrix<T> const &EXT,
                                               std::string const &ansProg,
                                               std::ostream &os) {
  return DirectFacetComputationInequalities(
      EXT, dual_desc_program_from_string(ansProg), os);
}

template <typename T, typename Fprocess>
void DirectFacetComputationFaceIneq(MyMatrix<T> const &EXT,
                                    DualDescProgram prog, Fprocess f_process,
                                    std::ostream &os) {
#ifdef DEBUG_DUAL_DESC
  os << "DDD: DirectFacetComputationFaceIneq, prog="
     << std::to_string(prog) << "\n";
#endif
  switch (prog) {
  case DualDescProgram::cdd:
    // The double description method, applicable to the field or ring case
    return double_desc::DualDescriptionFaceIneq(EXT, f_process, os);
  case DualDescProgram::small_polytopes:
    // Small polytopes can have special solutions
    return SmallPolytope_FaceIneq(EXT, f_process, os);
  case DualDescProgram::lrs:
    return lrs::DualDescriptionFaceIneq(EXT, f_process);
  case DualDescProgram::pd_lrs:
    // It applies to the field case or ring
    return POLY_DualDescription_PrimalDualFaceIneq(EXT, f_process, os);
  case DualDescProgram::beneath_beyond:
    // Native beneath-and-beyond, full-dimensional pointed cone (field only)
    if constexpr (is_ring_field<T>::value)
      return POLY_DualDescription_BeneathBeyondFaceIneq(EXT, f_process, os);
    break;
  case DualDescProgram::glrs:
  case DualDescProgram::ppl_ext:
  case DualDescProgram::cdd_ext:
  case DualDescProgram::normaliz:
    // The external programs are available only for rational types
#ifndef WASM_PLATFORM
    if constexpr (is_implementation_of_Q<T>::value) {
      if (prog == DualDescProgram::glrs)
        return DualDescExternalProgramFaceIneq(EXT, "glrs", f_process, os);
      if (prog == DualDescProgram::ppl_ext)
        return DualDescExternalProgramFaceIneq(EXT, "ppl_lcdd", f_process, os);
      if (prog == DualDescProgram::cdd_ext)
        return DualDescExternalProgramFaceIneq(EXT, "lcdd_gmp", f_process, os);
      if (prog == DualDescProgram::normaliz)
        return DualDescExternalProgramFaceIneq(EXT, "normaliz", f_process, os);
    }
#endif
    break;
  }
  terminate_direct_dual_desc("DirectFacetComputationFaceIneq", prog);
}

template <typename T, typename Fprocess>
void DirectFacetComputationFaceIneq(MyMatrix<T> const &EXT,
                                    std::string const &ansProg,
                                    Fprocess f_process, std::ostream &os) {
  return DirectFacetComputationFaceIneq(
      EXT, dual_desc_program_from_string(ansProg), f_process, os);
}

template <typename T, typename Tgroup>
vectface DirectFacetOrbitComputation(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                     DualDescProgram prog, std::ostream &os) {
#ifdef TIMINGS_DUAL_DESC
  MicrosecondTime time;
#endif
#ifdef KEY_VALUE_DUAL_DESC
  MicrosecondTime time_total;
#endif
  vectface ListIncd = DirectFacetComputationIncidence(EXT, prog, os);
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
     << EXT.cols() << "_" << GRP.size() << "_"
     << std::to_string(prog) << "_" << ListIncd.size() << "_"
     << TheOutput.size() << ") VALUE=(" << time_total << ")\n";
#endif
  return TheOutput;
}

template <typename T, typename Tgroup>
vectface DirectFacetOrbitComputation(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                     std::string const &ansProg,
                                     std::ostream &os) {
  return DirectFacetOrbitComputation(
      EXT, GRP, dual_desc_program_from_string(ansProg), os);
}

template <typename T, typename Tgroup>
std::vector<std::pair<Face, MyVector<T>>>
DirectFacetIneqOrbitComputation(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                DualDescProgram prog, std::ostream &os) {
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
  DirectFacetComputationFaceIneq(EXT, prog, f_process, os);
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
     << EXT.cols() << "_" << GRP.size() << "_"
     << std::to_string(prog) << "_" << ListReturn.size() << "_"
     << TheOutput.size() << ") VALUE=(" << time_total << ")\n";
#endif
  return TheOutput;
}

template <typename T, typename Tgroup>
std::vector<std::pair<Face, MyVector<T>>>
DirectFacetIneqOrbitComputation(MyMatrix<T> const &EXT, Tgroup const &GRP,
                                std::string const &ansProg, std::ostream &os) {
  return DirectFacetIneqOrbitComputation(
      EXT, GRP, dual_desc_program_from_string(ansProg), os);
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
