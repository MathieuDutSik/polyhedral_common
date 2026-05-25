// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for IsoDelaunay-domain enumeration.
// Models src_delaunay/LATT_SerialLattice_IsoDelaunayDomain.cpp at the
// type-instantiation level: builds a canonical T-space (sym n x n matrices),
// constructs DataIsoDelaunayDomain, and asks for the initial generic
// tesselation.  The full FullEnumeration sweep is intentionally skipped --
// goal is to verify IsoDelaunayDomains.h templates compile under WASM_PLATFORM.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "IsoDelaunayDomains.h"
#include "Tspace_Generation.h"
#include "Permutation.h"
#include "Group.h"
#include <iostream>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

int main() {
  using T = cpp_rational;
  using Tint = cpp_int;
  using Tidx = uint32_t;
  using Telt = permutalib::SingleSidedPerm<Tidx>;
  using TintGroup = cpp_int;
  using Tgroup = permutalib::Group<Telt, TintGroup>;

  // Canonical T-space of symmetric n x n matrices, n = 2.
  LinSpaceMatrix<T> LinSpa = ComputeCanonicalSpace<T>(2);
  int dimEXT = LinSpa.n + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, std::cerr);
  RecordDualDescOperation<T, Tgroup> rddo(AllArr, std::cerr);
  std::optional<MyMatrix<T>> CommonGramMat;

  DataIsoDelaunayDomains<T, Tint, Tgroup> data{
    LinSpa, std::move(rddo), CommonGramMat};
  DelaunayTesselation<T, Tgroup> initial =
      GetInitialGenericDelaunayTesselation<T, Tint, Tgroup>(data);
  std::cerr << "Test_wasm_isodelaunay: |initial cells|="
            << initial.l_dels.size() << " OK\n";
  return 0;
}
