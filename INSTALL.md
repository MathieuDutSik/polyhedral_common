Installation
============

Access to the source code
-------------------------

This repository uses submodules, so the recommended clone command is

```sh
git clone https://github.com/MathieuDutSik/polyhedral_common --recursive
```

If you already cloned the repository without `--recursive`, run

```sh
./init.sh
```

To update the submodules later, run

```sh
./update.sh
```

Overview
--------

There are three supported ways to compile the project:

1. Docker: use the prebuilt build environment from
   `docker_files/docker_complete/Dockerfile`.
2. CMake: build the executables from the top-level `CMakeLists.txt`.
3. Legacy Makefiles: compile each subsystem with the historical per-directory
   Makefiles. This is still the standard path used by `compile.sh` and by CI.

Dependencies
------------

The project uses the following external libraries:

* Eigen: dense linear algebra templates.
* GMP / GMPXX: exact integer and rational arithmetic.
* Boost.
  * `Boost.Serialization` is used broadly for serialization of matrices,
    combinatorial objects, and intermediate data structures.
  * `Boost.MPI` is used by the MPI-enabled programs in the dual description,
    Delaunay, perfect form, Lorentzian, C-type, and enumeration code.
* OpenMPI or another compatible MPI implementation: runtime and compiler
  support for the MPI-enabled binaries.
* GLPK: linear programming support.
* cddlib: double description computations, both rational and floating-point
  variants.
* FLINT: number-theoretic and arithmetic support in several subsystems.
* netCDF C / C++: used by the C-type tools.
* nauty / Traces: graph automorphism and canonical form computations.

Optional library
----------------

The codebase also contains support for the `bliss` graph automorphism library.
Historically, some graph routines could be run either with `bliss` or with
`Traces`. The default path in the current code is `Traces`, so `bliss` is not
required for the standard build, and the current CMake configuration does not
depend on it.

Method 1: Docker
----------------

The most complete way to get a working environment is to use
`docker_files/docker_complete/Dockerfile`.

Typical workflow:

```sh
cd docker_files/docker_complete
docker build -t polyhedralcpp:build-env .
docker run -it polyhedralcpp:build-env /bin/bash
```

Inside the image, the Dockerfile installs the external dependencies and clones
the repository with submodules. It then runs `./compile.sh` to build the full
legacy Makefile-based toolchain.

Use this method when:

* you want a reproducible Ubuntu-based environment,
* you do not want to install all dependencies on the host,
* or you want a reference setup close to the maintainer environment.

Method 2: CMake
---------------

The repository now provides a top-level `CMakeLists.txt` that can build the
project executables directly.

Typical workflow:

```sh
mkdir -p build
cd build
cmake ..
make -j
```

Notes:

* The CMake build creates the binaries under `build/bin/`.
* CMake discovers the system libraries with `find_package` and `pkg-config`.
* The CMake configuration also arranges for `nauty` to be built as an external
  dependency when needed.
* This is the most convenient native build method if you want a single build
  directory and standard CMake tooling.

Method 3: Legacy Makefiles
--------------------------

The historical build system is a collection of per-directory Makefiles. This is
still the build path used by the helper scripts and by CI.

To build the main collection of binaries, run

```sh
./compile.sh
```

This script enters the main source directories one by one and runs `make clean`
followed by `make`.

There is also a reduced variant

```sh
./compile_docker.sh
```

which builds a smaller subset of directories. It is mainly useful inside the
Docker environment.

Use the legacy Makefiles when:

* you want to match the existing CI path exactly,
* you are debugging a single subsystem Makefile,
* or you want to build directly inside a source subdirectory.

Remarks
-------

* The native Makefile-based builds place executables directly inside the
  corresponding `src_*` directories.
* The CMake build keeps generated files inside the build tree.
* If you switch between build methods, cleaning old binaries can avoid
  confusion.
