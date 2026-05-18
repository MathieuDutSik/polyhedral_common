#!/bin/bash
set -euo pipefail

# Smoke-build every Test_wasm_*.cpp in this directory under WebAssembly via
# emscripten. Each test #includes the public headers of one src_* area to
# verify that the area's templates compile cleanly under WASM_PLATFORM
# (which disables external-program invocation, signal handling, MPI, etc.).
#
# Requires emscripten (emcc, emconfigure, emmake) on PATH.
# On macOS:  brew install emscripten
#
# nauty is built once into build/nauty-install/ by running its upstream
# autoconf-generated configure under emconfigure; the resulting libnauty.a
# is linked into the smoke tests that pull in GRAPH_traces.h (currently
# only Test_wasm_short via SHORT_Realizability).

cd "$(dirname "$0")"

BOOST_INC="${BOOST_INC:-/opt/homebrew/opt/boost/include}"
EIGEN_INC="${EIGEN_INC:-/opt/homebrew/include/eigen3}"
NAUTY_TAG="${NAUTY_TAG:-2.9.3}"
NAUTY_REPO="${NAUTY_REPO:-https://github.com/MathieuDutSik/nauty}"

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
WORKDIR="$(pwd)/build"
NAUTY_SRC="$WORKDIR/nauty-src"
NAUTY_INSTALL="$WORKDIR/nauty-install"
NAUTY_LIB="$NAUTY_INSTALL/lib/libnauty.a"
# nauty's `make install` drops the headers directly in $prefix/include
# (traces.h, nauty.h, ...), not in a `nauty/` subdir.
NAUTY_INC="$NAUTY_INSTALL/include"
mkdir -p "$WORKDIR"

ensure_wasm_nauty() {
  if [ -f "$NAUTY_LIB" ]; then
    echo "==> nauty (wasm) already built: $NAUTY_LIB"
    return
  fi
  if [ ! -d "$NAUTY_SRC" ]; then
    echo "==> Cloning nauty $NAUTY_TAG from $NAUTY_REPO"
    git clone --depth 1 --branch "$NAUTY_TAG" "$NAUTY_REPO" "$NAUTY_SRC"
  fi
  echo "==> Building nauty under emscripten in $NAUTY_SRC"
  (
    cd "$NAUTY_SRC"
    # --host=wasm32-unknown-emscripten forces autoconf into cross-compile
    # mode; without it configure picks up the local triple (e.g.
    # aarch64-apple-darwin) and runs Apple-specific linker tests
    # (-single_module, -force_load) that wasm-ld rejects.
    # --disable-popcnt/clz/tls suppress CPU/TLS feature autodetection,
    # which would otherwise try to run a probe binary on the target.
    emconfigure ./configure \
      --host=wasm32-unknown-emscripten \
      --prefix="$NAUTY_INSTALL" \
      --disable-popcnt \
      --disable-clz \
      --disable-tls
    emmake make
    emmake make install
  )
  if [ ! -f "$NAUTY_LIB" ]; then
    echo "ERROR: $NAUTY_LIB not produced by emmake install" >&2
    exit 1
  fi
}

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
  -DWASM_PLATFORM
)

# Tests that pull in GRAPH_traces.h (and therefore need libnauty linked).
# Anything going through the equivalence / canonicalization path
# (WeightMatrix → GRAPH_Bindings → GRAPH_traces) ends up here:
#   - Test_wasm_short      (SHORT_Realizability)
#   - Test_wasm_polytope_aut (LinPolytope_Automorphism)
#   - Test_wasm_gram_aut   (ArithmeticAutomorphismGroupMultiple)
needs_nauty() {
  case "$1" in
    Test_wasm_short.cpp|Test_wasm_polytope_aut.cpp|Test_wasm_gram_aut.cpp)
      return 0 ;;
    *) return 1 ;;
  esac
}

# Build nauty up-front if any source needs it, so the failure (if any)
# happens before the per-test loop.
shopt -s nullglob
sources=(Test_wasm_*.cpp)
shopt -u nullglob
if [ "${#sources[@]}" -eq 0 ]; then
  echo "No Test_wasm_*.cpp sources found." >&2
  exit 1
fi
need_nauty=0
for src in "${sources[@]}"; do
  if needs_nauty "$src"; then
    need_nauty=1
    break
  fi
done
if [ "$need_nauty" -eq 1 ]; then
  ensure_wasm_nauty
fi

n_fail=0
for src in "${sources[@]}"; do
  name="${src%.cpp}"
  out="build/${name}.js"
  extra=()
  if needs_nauty "$src"; then
    extra=(-I"$NAUTY_INC" "$NAUTY_LIB")
  fi
  echo "==> Building $src -> $out"
  # ${extra[@]+"${extra[@]}"} expands to nothing when `extra` is empty;
  # plain "${extra[@]}" tripped `set -u` (treats empty-array as unset).
  if ! emcc "${CXXFLAGS[@]}" "${INCLUDES[@]}" ${extra[@]+"${extra[@]}"} "$src" -o "$out"; then
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
