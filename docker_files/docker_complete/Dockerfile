# polyhedral_common — C++ polyhedral computation toolkit.
# This image contains the full build environment AND the compiled binaries,
# ready to run.
#
# Build (from the repository root):
#   docker build -t polyhedralcpp:build-env --file docker_files/docker_complete/Dockerfile .
#
# Run — drops you into a shell in the repo, with every binary already on PATH:
#   docker run -it polyhedralcpp:build-env
#   POLY_dual_description --help        # e.g. any of the POLY_/LATT_/GRP_ tools
#
# Prebuilt images are published to: ghcr.io/mathieudutsik/polyhedral_common
#   docker pull ghcr.io/mathieudutsik/polyhedral_common:latest
#
# Maintainer: Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

FROM ubuntu:24.04

LABEL maintainer="mathieu.dutour@gmail.com"
LABEL org.opencontainers.image.source="https://github.com/MathieuDutSik/polyhedral_common"
LABEL org.opencontainers.image.description="C++ polyhedral computation toolkit (build environment + compiled binaries)"

# noninteractive so tzdata and similar never block the build waiting on a prompt.
ENV DEBIAN_FRONTEND=noninteractive

# ---------------------------------------------------------------------------
# System packages. Installed in a single layer and the apt lists are removed
# afterwards so they do not bloat the image.
#   toolchain  : g++/gcc/gfortran/clang, make, cmake, pkg-config, autotools, git, wget
#   linker     : mold (see the dpkg-divert step below)
#   libraries  : openmpi, gmp, boost(+mpi,+serialization), flint(+ntl),
#                netcdf-c++4 (required by src_ctype)
#   convenience: sudo + a small editor + less for interactive use
# ---------------------------------------------------------------------------
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      ca-certificates wget git file \
      make cmake pkg-config autoconf automake libtool \
      g++ gcc gfortran clang \
      mold \
      openmpi-bin libopenmpi-dev \
      libgmp-dev \
      libboost-dev libboost-mpi-dev libboost-serialization-dev \
      libflint-dev libntl-dev \
      libnetcdf-c++4-dev \
      libgdbm-dev libsnappy-dev zlib1g-dev \
      sudo less nano rlwrap \
 && rm -rf /var/lib/apt/lists/*

# Replace GNU ld with mold. Linking the heavy MPI executables (e.g.
# LATT_MPI_Lattice_IsoDelaunayDomain) with -g -O3 + static Boost.MPI + bliss +
# nauty peaks at several GB per link — enough that GNU ld gets OOM-killed inside
# Docker Desktop's default VM. mold is a drop-in replacement using 2-5x less
# memory (and it is faster). dpkg-divert + symlink overrides both the `ld` and
# `ld.bfd` codepaths in gcc's collect2.
RUN dpkg-divert --divert /usr/bin/ld.orig --rename /usr/bin/ld \
 && ln -sf /usr/bin/mold /usr/bin/ld \
 && dpkg-divert --divert /usr/bin/ld.bfd.orig --rename /usr/bin/ld.bfd \
 && ln -sf /usr/bin/mold /usr/bin/ld.bfd

# ---------------------------------------------------------------------------
# Build-flag environment consumed by the polyhedral_common Makefiles.
# ---------------------------------------------------------------------------
ENV GMP_INCDIR="/usr/include" \
    GMP_CXX_LINK="-lgmp -lgmpxx"

ENV BOOST_INCDIR="/usr/include" \
    BOOST_LINK="-lboost_serialization" \
    MPI_LINK_CPP="-lboost_mpi"

ENV FLINT_INCDIR="/usr/include" \
    FLINT_LINK="-lflint"

# netcdf-c++4 from the Ubuntu package: headers under /usr/include, the library
# is libnetcdf_c++4 on the default linker search path.
ENV NETCDF_CXX_ALLINC="-I/usr/include" \
    NETCDF_CXX_LINK="-lnetcdf_c++4"

# Eigen is header-only. Pinned to 3.3.9 (the version the CI uses and is known
# good for the exact/algebraic scalar types; 3.4.0 regresses on them) and
# shallow-cloned to avoid pulling the full history. Change --branch to move it.
RUN git clone --depth 1 --branch 3.3.9 https://gitlab.com/libeigen/eigen.git /opt/eigen
ENV EIGEN_PATH="/opt/eigen"

# ---------------------------------------------------------------------------
# polyhedral_common source (with submodules) plus the two build-time deps that
# are not packaged: bliss and nauty.
# ---------------------------------------------------------------------------
RUN git clone --depth 1 --recursive --shallow-submodules \
      https://github.com/MathieuDutSik/polyhedral_common.git /GIT/polyhedral_common

# bliss (shallow). LIBBLISS_LINK rpaths the in-tree build directory, so it is
# kept in the image.
RUN git clone --depth 1 https://github.com/MathieuDutSik/bliss /GIT/bliss \
 && cmake -S /GIT/bliss -B /GIT/bliss/build \
 && cmake --build /GIT/bliss/build
ENV LIBBLISS_INCDIR="/GIT/bliss/src" \
    LIBBLISS_LINK="-L/GIT/bliss/build -lbliss -Wl,-rpath,/GIT/bliss/build"

# nauty/traces from the vendored copy, installed into /opt/nauty; the build
# directory is removed once installed.
ENV NAUTY_PATH="/opt/nauty" \
    NAUTY_INCLUDE="-I/opt/nauty/include" \
    NAUTY_LINK="-L/opt/nauty/lib -lnauty_shared -Wl,-rpath,/opt/nauty/lib"
RUN cmake -S /GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty \
          -B /GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty/build \
          -DCMAKE_INSTALL_PREFIX=/opt/nauty \
 && cmake --build /GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty/build \
 && cmake --install /GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty/build \
 && rm -rf /GIT/polyhedral_common/basic_common_cpp/ExternalLib/nauty/build

# ---------------------------------------------------------------------------
# Compile everything, then expose the binaries on PATH for convenience.
# ---------------------------------------------------------------------------
WORKDIR /GIT/polyhedral_common
RUN ./compile.sh

# Symlink every compiled ELF binary into /usr/local/bin so the tools can be
# invoked by name (POLY_*, LATT_*, GRP_*, …) from any directory.
RUN find /GIT/polyhedral_common/src_* -maxdepth 1 -type f -perm -u+x \
      -exec sh -c 'file -b "$1" | grep -q ELF && ln -sf "$1" /usr/local/bin/' _ {} \;

CMD ["/bin/bash"]
