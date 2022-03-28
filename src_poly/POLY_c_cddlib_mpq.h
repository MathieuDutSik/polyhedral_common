#ifndef SRC_POLY_POLY_C_CDDLIB_MPQ_H_
#define SRC_POLY_POLY_C_CDDLIB_MPQ_H_

#define GMPRATIONAL

#include "Boost_bitset_kernel.h"
#include "MatrixTypes.h"
#include "gmpxx.h"
#include <boost/multiprecision/gmp.hpp>

namespace cbased_cdd {
vectface DualDescription_incd_mpq_class(MyMatrix<mpq_class> const &TheEXT);
vectface DualDescription_incd_boost_mpq_rational(
    MyMatrix<boost::multiprecision::mpq_rational> const &TheEXT);
} // namespace cbased_cdd

#endif //  SRC_POLY_POLY_C_CDDLIB_MPQ_H_
