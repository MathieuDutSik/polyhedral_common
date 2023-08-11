// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_POLY_POLY_C_CDDLIB_MPQ_H_
#define SRC_POLY_POLY_C_CDDLIB_MPQ_H_

#define GMPRATIONAL

// clang-format off
#include "Boost_bitset_kernel.h"
#include "MatrixTypes.h"
#ifndef DISABLE_MPQ_CLASS
#include "gmpxx.h"
#endif
#include <boost/multiprecision/gmp.hpp>
// clang-format on

namespace cbased_cdd {
#ifndef DISABLE_MPQ_CLASS
vectface DualDescription_incd_mpq_class(MyMatrix<mpq_class> const &TheEXT);
#endif
vectface DualDescription_incd_boost_mpq_rational(
    MyMatrix<boost::multiprecision::mpq_rational> const &TheEXT);
// clang-format off
}  // namespace cbased_cdd
// clang-format on

// clang-format off
#endif  // SRC_POLY_POLY_C_CDDLIB_MPQ_H_
// clang-format on
