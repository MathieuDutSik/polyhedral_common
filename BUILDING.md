Building
========

Access to the source code
-------------------------

Since this repository uses submodules, the cloning command is

```sh
$ git clone https://github.com/MathieuDutSik/polyhedral_common --recursive
```

If you cloned but forgot to put the `--recursive` then run `./init.sh` so as
to get the subrepositories.

In order to update the submodule the command is `./update.sh`.

Compilation
-----------

The code depends on only a few libraries: Eigen, Boost and nauty.

Several ways to compile the program are made available:
* One is to compile via  **Dockerfile**. For the complete code, it is available at **docker_files/docker_complete/Dockerfile**.
* Another way is to use the **CMakeLists.txt** for building the Makefile and compiling.
* The standard way used in CI is to use the **Makefile**.


Dependencies
------------

Following dependencies are needed for compiling the code:

  * Eigen: http://eigen.tuxfamily.org/
  * Boost: http://www.boost.org/
  * GNU MultiPrecision Library (GMP): https://gmplib.org/
  * nauty : https://pallini.di.uniroma1.it/
