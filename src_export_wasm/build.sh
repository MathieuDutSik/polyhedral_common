#!/bin/bash
set -euo pipefail

# Smoke-build every Test_wasm_*.cpp in this directory under WebAssembly via
# emscripten. Each test #includes the public headers of one src_* area to
# verify that the area's templates compile cleanly under POLYHEDRAL_WASM
# (which disables external-program invocation, signal handling, MPI, etc.).
#
# Requires emscripten (emcc) on PATH. On macOS:  brew install emscripten
#
# Companion harness in ../basic_common_cpp/src_wasm/ builds the lower-level
# matrix/number-theory smoke tests; this directory exercises the higher-level
# polyhedral areas: src_poly, src_copos, src_isotropy, src_group, src_latt,
# src_polygen, src_short, src_sparse_solver.

cd "$(dirname "$0")"

BOOST_INC="${BOOST_INC:-/opt/homebrew/opt/boost/include}"
EIGEN_INC="${EIGEN_INC:-/opt/homebrew/include/eigen3}"

if [ ! -d "$BOOST_INC" ]; then
  echo "Boost include dir not found at $BOOST_INC (set BOOST_INC env var)" >&2
  exit 1
fi
if [ ! -d "$EIGEN_INC" ]; then
  echo "Eigen include dir not found at $EIGEN_INC (set EIGEN_INC env var)" >&2
  exit 1
fi

ROOT="$(cd .. && pwd)"
BCPP="$ROOT/basic_common_cpp"

INCLUDES=(
  # Basic common (lower-level)
  -I"$BCPP/src_basic"
  -I"$BCPP/src_number"
  -I"$BCPP/src_matrix"
  -I"$BCPP/src_comb"
  -I"$BCPP/src_graph"
  -I"$BCPP/sparse-map/include/tsl"
  -I"$BCPP/robin-map/include/tsl"
  -I"$BCPP/hopscotch-map/include/tsl"
  # Polyhedral common (target areas + their internal deps)
  -I"$ROOT/src_poly"
  -I"$ROOT/src_copos"
  -I"$ROOT/src_isotropy"
  -I"$ROOT/src_group"
  -I"$ROOT/src_latt"
  -I"$ROOT/src_polygen"
  -I"$ROOT/src_short"
  -I"$ROOT/src_sparse_solver"
  -I"$ROOT/src_lorentzian"
  -I"$ROOT/permutalib/src"
  # External
  -I"$BOOST_INC"
  -I"$EIGEN_INC"
)

CXXFLAGS=(
  -std=c++20
  -O2
  --minify
  0
  -Wall
  -Wextra
  -Wno-deprecated-declarations
  -D_LIBCPP_ENABLE_CXX20_REMOVED_TYPE_TRAITS
  -DINCLUDE_NUMBER_THEORY_BOOST_CPP_INT
  -DPOLYHEDRAL_WASM
)

mkdir -p build
shopt -s nullglob
sources=(Test_wasm_*.cpp)
shopt -u nullglob

if [ "${#sources[@]}" -eq 0 ]; then
  echo "No Test_wasm_*.cpp sources found." >&2
  exit 1
fi

n_fail=0
for src in "${sources[@]}"; do
  name="${src%.cpp}"
  out="build/${name}.js"
  echo "==> Building $src -> $out"
  if ! emcc "${CXXFLAGS[@]}" "${INCLUDES[@]}" "$src" -o "$out"; then
    echo "    BUILD FAILED for $src" >&2
    n_fail=$((n_fail + 1))
    continue
  fi
  echo "==> Running $out under node"
  if ! node "$out"; then
    echo "    RUN FAILED for $src" >&2
    n_fail=$((n_fail + 1))
  fi
done

if [ "$n_fail" -ne 0 ]; then
  echo "$n_fail test(s) failed." >&2
  exit 1
fi
echo "All Wasm tests passed."
