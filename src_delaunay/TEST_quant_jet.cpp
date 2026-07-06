// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Compile-feasibility probe: force instantiation of the quantization integral
// over jet<T, N>, to discover which subroutines are (not) jet-able.
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "QuantizationIntegral.h"
#include "Permutation.h"
#include "Group.h"
// clang-format on

using T = mpq_class;
using Tint = mpz_class;
using Tidx = uint32_t;
using Telt = permutalib::SingleSidedPerm<Tidx>;
using Tint_grp = mpz_class;
using Tgroup = permutalib::Group<Telt, Tint_grp>;

// Not meant to run; only to force template instantiation of the jet path.
[[maybe_unused]] MyMatrix<jet<T, 2>> probe(DelaunayTesselation<T, Tgroup> const &DT,
                                           MyMatrix<T> const &SHV,
                                           MyMatrix<jet<T, 2>> const &Gram) {
  return QuantizationSecMomentMatJet<T, 2, Tint, Tgroup>(DT, SHV, Gram, 3,
                                                         T(1) / T(2), std::cerr);
}

int main() { return 0; }
