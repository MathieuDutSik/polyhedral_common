// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test for Delaunay-polytope enumeration on a lattice.
// Models src_delaunay/LATT_SerialComputeDelaunay.cpp at the type-instantiation
// level only: builds DataLattice + PolyHeuristicSerial and computes ONE
// initial Delaunay polytope via FindDelaunayPolytopeExtended.  Skipping the
// full EnumerateAndStore_Serial sweep keeps the smoke test fast; the goal is
// just to confirm the LatticeDelaunay templates compile under WASM_PLATFORM.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "LatticeDelaunay.h"
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

  // Z^2 with the standard inner product.
  MyMatrix<T> GramMat = IdentityMat<T>(2);
  int dimEXT = GramMat.rows() + 1;
  PolyHeuristicSerial<TintGroup> AllArr =
      AllStandardHeuristicSerial<T, TintGroup>(dimEXT, std::cerr);
  DataLattice<T, Tint, Tgroup> data_lattice =
      GetDataLattice<T, Tint, Tgroup>(GramMat, AllArr, std::cerr);
  MyMatrix<Tint> EXT = FindDelaunayPolytopeExtended<T, Tint, Tgroup>(data_lattice);
  std::cerr << "Test_wasm_delaunay: |EXT|=" << EXT.rows() << "x" << EXT.cols()
            << " OK\n";
  return 0;
}
