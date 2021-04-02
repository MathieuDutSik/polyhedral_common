#define GMPRATIONAL

#include "gmpxx.h"

#include <vector>
#include "Boost_bitset.h"
#include "MAT_Matrix.h"

namespace cbased_cdd {

std::vector<Face> DualDescription_incd_mpq(MyMatrix<mpq_class> const& TheEXT);

}
