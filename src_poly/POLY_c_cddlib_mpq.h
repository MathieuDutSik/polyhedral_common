#define GMPRATIONAL

#include "gmpxx.h"
#include <boost/multiprecision/gmp.hpp>
#include "Boost_bitset_kernel.h"
#include "MatrixTypes.h"

namespace cbased_cdd {

  vectface DualDescription_incd_mpq_class(MyMatrix<mpq_class> const& TheEXT);

#ifdef INCLUDE_NUMBER_THEORY_BOOST_GMP_INT
  vectface DualDescription_incd_boost_mpq_rational(MyMatrix<boost::multiprecision::mpq_rational> const& TheEXT);
#endif

}
