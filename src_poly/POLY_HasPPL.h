// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_HASPPL_H_
#define SRC_POLY_POLY_HASPPL_H_

// Compile-time + run-time gate for the external `ppl_lcdd` binary.  Returns
// true only when (a) the arithmetic T is a rational implementation and (b)
// `ppl_lcdd` is on PATH.  Under POLYHEDRAL_WASM the answer is always false
// because the WASM build cannot spawn external programs.

// clang-format off
#ifndef POLYHEDRAL_WASM
#include "Basic_file.h"
#endif
// clang-format on

#ifdef POLYHEDRAL_WASM

template <typename T> inline bool IsPPLpossible() { return false; }

#else

template <typename T>
  requires(is_implementation_of_Q<T>::value)
inline bool IsPPLpossible() {
  return IsProgramInPath("ppl_lcdd");
}

template <typename T>
  requires(!is_implementation_of_Q<T>::value)
inline bool IsPPLpossible() {
  return false;
}

#endif

// clang-format off
#endif  // SRC_POLY_POLY_HASPPL_H_
// clang-format on
