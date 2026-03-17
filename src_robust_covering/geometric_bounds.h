// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_ROBUST_COVERING_GEOMETRIC_BOUNDS_H_
#define SRC_ROBUST_COVERING_GEOMETRIC_BOUNDS_H_


/*
  combined
  upper_bound <= (sqrt(bound1) + sqrt(bound2) )^2
    <= bound1 + bound2 + 2 sqrt(bound1 * bound2)
  If bound1 <= bound2 We bound it by
    bound1 + 3 bound2
 */
template <typename T> T combine_two_bounds(T const &bound1, T const &bound2) {
  auto test_correctness = [&](T const &val) -> bool {
    // Correct if sqrt(bound1) + sqrt(bound2) <= sqrt(val)
    // Equivalent to bound1 + bound2 + 2 sqrt(bound1 bound2) <= val
    // Equivalent to 2 sqrt(bound1 bound2) <= val - bound1 - bound2
    // Equivalent to 4 bound1 bound2 <= (val - bound1 - bound2)^2
    if (val < 0) {
      // Should not occur really.
      return false;
    }
    T delta = val - bound1 - bound2;
    if (delta < 0) {
      return false;
    }
    T val1 = 4 * bound1 * bound2;
    T val2 = delta * delta;
    if (val1 <= val2) {
      return true;
    } else {
      return false;
    }
  };
  T low(0);
  T upp(1);
  while (true) {
    bool test_low = test_correctness(low);
    bool test_upp = test_correctness(upp);
    if (!test_low && test_upp) {
      break;
    }
    low += 1;
    upp += 1;
  }
  for (int i = 0; i < 6; i++) {
    T mid = (low + upp) / 2;
    bool test_mid = test_correctness(mid);
    if (test_mid) {
      upp = mid;
    } else {
      low = mid;
    }
  }
  return upp;
}



// clang-format off
#endif  // SRC_ROBUST_COVERING_GEOMETRIC_BOUNDS_H_
// clang-format on
