// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_RANKIN_ENUMERATION_K_SPACE_H_
#define SRC_RANKIN_ENUMERATION_K_SPACE_H_

// Returns the Hermite constant at the n-th power.
// That is the minimum of min(A)^n / det(A)
template<typename T>
T GetUpperBoundHermitePower(int n) {
  if (n == 1) {
    return 1;
  }
  if (n == 2) {
    return 4/3;
  }
  if (n == 3) {
    return 2;
  }
  if (n == 4) {
    return 4;
  }
  if (n == 5) {
    return 8;
  }
  if (n == 6) {
    return 64/3;
  }
  if (n == 7) {
    return 64;
  }
  if (n == 8) {
    return 256;
  }
  int h = n*(n-1)/2;
  return (4/3)^h;
}

// clang-format off
#endif  // SRC_RANKIN_ENUMERATION_K_SPACE_H_
// clang-format on
