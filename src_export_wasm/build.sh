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
# boost_serialization is the only non-header-only piece of Boost that the
# polytope/lattice code needs (extended_type_info registry + archive base
# classes referenced by every `serialize()` member). Build it from the
# modular repo with emcc directly -- b2 cross-compile is overkill for the
# ~35 .cpp files we need.
BOOST_SERIAL_TAG="${BOOST_SERIAL_TAG:-boost-1.85.0}"
BOOST_SERIAL_REPO="${BOOST_SERIAL_REPO:-https://github.com/boostorg/serialization}"
BOOST_SERIAL_SRC="$WORKDIR/boost-serialization-src"
BOOST_SERIAL_LIB="$WORKDIR/libboost_serialization.a"
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

ensure_wasm_boost_serialization() {
  if [ -f "$BOOST_SERIAL_LIB" ]; then
    echo "==> boost_serialization (wasm) already built: $BOOST_SERIAL_LIB"
    return
  fi
  if [ ! -d "$BOOST_SERIAL_SRC" ]; then
    echo "==> Cloning boost_serialization $BOOST_SERIAL_TAG from $BOOST_SERIAL_REPO"
    git clone --depth 1 --branch "$BOOST_SERIAL_TAG" \
      "$BOOST_SERIAL_REPO" "$BOOST_SERIAL_SRC"
  fi
  echo "==> Building boost_serialization under emscripten"
  (
    cd "$BOOST_SERIAL_SRC/src"
    # BOOST_*_SOURCE defines tell the headers not to emit MSVC dllimport
    # attributes; BOOST_ALL_NO_LIB suppresses the auto-link pragmas that
    # would otherwise try to pull a Windows-only .lib symbol.
    # -I points at the modular repo's own headers first, then at the
    # full Boost include tree for transitive deps (throw_exception, io,
    # core, etc.).
    # -idirafter for $BOOST_INC: when BOOST_INC points at /usr/include (as in
    # the Linux CI), plain -I would inject the host glibc tree ahead of the
    # wasi-libc sysroot, so libc++'s `#include_next <stdlib.h>` would land on
    # /usr/include/stdlib.h and fail on glibc-only bits/libc-header-start.h.
    # -idirafter appends BOOST_INC strictly after the sysroot search path.
    emcc -std=c++20 -O2 \
      -I"$BOOST_SERIAL_SRC/include" \
      -idirafter "$BOOST_INC" \
      -DBOOST_ARCHIVE_SOURCE \
      -DBOOST_WARCHIVE_SOURCE \
      -DBOOST_SERIALIZATION_SOURCE \
      -DBOOST_ALL_NO_LIB \
      -c ./*.cpp
    emar rcs "$BOOST_SERIAL_LIB" ./*.o
  )
  if [ ! -f "$BOOST_SERIAL_LIB" ]; then
    echo "ERROR: $BOOST_SERIAL_LIB not produced" >&2
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
  -I"$ROOT/src_dualdesc"
  -I"$ROOT/src_enum_schemes"
  -I"$ROOT/src_polydecomp"
  -I"$ROOT/src_perfect"
  -I"$ROOT/src_delaunay"
  -I"$ROOT/src_sparse_solver"
  -I"$ROOT/src_lorentzian"
  -I"$ROOT/permutalib/src"
  # External: -idirafter (not -I) so the host include tree never shadows the
  # wasm sysroot. See ensure_wasm_boost_serialization above for details.
  -idirafter "$BOOST_INC"
  -idirafter "$EIGEN_INC"
)

# --- Optimization flags ---
#
# emscripten honors all the standard clang levels plus a couple of WASM-
# specific extras. Pick based on what matters:
#
#   -O0                  no optimization; fastest build, biggest+slowest .wasm
#   -O1 / -O2 / -O3      usual clang levels (-O3 == max speed)
#   -Os                  optimize for size; modest speed loss vs -O2
#   -Oz                  emscripten/clang extension; more aggressive than -Os.
#                        Smallest .wasm achievable without LTO.
#   -flto                link-time optimization. Biggest single-shot size win
#                        for template-heavy C++ (cross-TU inlining, dead-code
#                        elimination); typically 30-60% smaller .wasm.
#   -g                   keep DWARF debug info, function names in stack traces.
#
# emscripten runs wasm-opt internally after compile; the level it picks
# tracks the -O flag here, so -O3/-Os/-Oz also enable progressively more
# aggressive post-link wasm-opt passes.
#
# Typical recipes:
#   web deploy (size matters): -Oz -flto       (+ --closure 1 for JS glue)
#   node/server (speed):       -O3 -flto       (current setting)
#   debugging:                 -O0 -g
#
# Note: --minify 0 below applies only to the JS glue file emcc emits, not
# to the .wasm itself; it keeps the glue readable for debugging.
#
# Caveat for -flto: every input object has to be built with it. If you want
# libnauty.a included in LTO, pass CFLAGS=-flto to emconfigure inside
# ensure_wasm_nauty as well; otherwise libnauty links in but is skipped by
# the LTO pass (still works, just less compact).
CXXFLAGS=(
  -std=c++20
  -O3
  -flto
  --minify
  0
  -Wall
  -Wextra
  -Wno-deprecated-declarations
  -D_LIBCPP_ENABLE_CXX20_REMOVED_TYPE_TRAITS
  -DWASM_PLATFORM
)


# Tests that pull in the equivalence/canonicalization path and so need BOTH
# libnauty linked AND a wasm build of libboost_serialization (every
# Tspace/Delaunay/Perfect header defines `serialize()` members that
# reference the boost-archive registry). The path is
# WeightMatrix → GRAPH_Bindings → GRAPH_traces (for nauty) and
# any serialize() member function (for boost_serialization). The two
# always travel together in this code base, so a single flag covers both:
#   - Test_wasm_short        (SHORT_Realizability)
#   - Test_wasm_polytope_aut (LinPolytope_Automorphism)
#   - Test_wasm_gram_aut     (ArithmeticAutomorphismGroupMultiple)
#   - Test_wasm_delaunay     (Delaunay_Stabilizer → polytope equivalence)
#   - Test_wasm_isodelaunay  (initial-tesselation builder → Delaunay equivalence)
#   - Test_wasm_lorentz_aut  (LORENTZ_GetGeneratorsAutom → polytope automorphism)
#   - Test_wasm_perfect      (DataPerfectTspaceFunc → polytope equivalence)
needs_nauty() {
  case "$1" in
    Test_wasm_short.cpp|Test_wasm_polytope_aut.cpp|Test_wasm_gram_aut.cpp|\
Test_wasm_delaunay.cpp|Test_wasm_isodelaunay.cpp|Test_wasm_lorentz_aut.cpp|\
Test_wasm_perfect.cpp)
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
  ensure_wasm_boost_serialization
fi

n_fail=0
for src in "${sources[@]}"; do
  name="${src%.cpp}"
  out="build/${name}.js"
  extra=()
  if needs_nauty "$src"; then
    extra=(-I"$NAUTY_INC" "$NAUTY_LIB" "$BOOST_SERIAL_LIB")
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
