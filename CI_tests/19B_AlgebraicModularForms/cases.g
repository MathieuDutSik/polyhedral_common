# The cases of TestAlgebraicModularForms.g: the classes of the genus (all of
# them), the level, the Hecke primes, and the expected result
# expected_<name>.g. With mpi:=true, the case is also run on two MPI
# processes, which must give the same result.
#
# - E8, level 1, p = 2, 3: the numbers of neighbors, 135 and 1120.
# - E8, level 5, p = 2, 3: the two orbits of Aut(E8) on the isotropic lines
#   mod 5, of sizes 7560 and 12096.
# - II_16 (E8+E8 and D16+), level 1, p = 2: eigenvalues 32895 and 1800.
# - The 25 classes of the genus 4.4.945.945 of DefLatDB, level 1 at
#   p = 2, 11, and level 11 at p = 2, 13 (1365 states).
return [
  rec(name:="E8_L1", classes:="classes_E8.txt", level:=1, primes:=[2, 3], mpi:=false),
  rec(name:="E8_L5", classes:="classes_E8.txt", level:=5, primes:=[2, 3], mpi:=false),
  rec(name:="II16_L1", classes:="classes_II16.txt", level:=1, primes:=[2], mpi:=false),
  rec(name:="c04293_L1", classes:="classes_c04293.txt", level:=1, primes:=[2, 11], mpi:=false),
  rec(name:="c04293_L11", classes:="classes_c04293.txt", level:=11, primes:=[2, 13], mpi:=true)
];
