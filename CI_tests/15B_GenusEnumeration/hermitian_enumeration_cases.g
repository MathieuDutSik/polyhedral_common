# The cases of TestHermitianEnumeration.g: genera of positive definite
# Hermitian lattices over CM fields, by their genus symbols, with the class
# numbers and the orders of the unitary groups of Hecke.jl
# (genus_representatives, automorphism_group_order, mass). The field is
# abelian (conductor, generators of H; the mass computed by the program) or
# given by a polynomial (the mass given). They cover split, inert and
# ramified primes, the dyadic ramified primes (Q(i), Q(sqrt -2), Q(zeta_8)),
# fields of degree 2, 4 and 6, genera of up to 25 classes, two primes of K_0
# above 11 in Q(zeta_5), and Q(sqrt -7) given by two different polynomials.
return [
  rec(name:="Q(sqrt-3) 1,3,6", field:=rec(conductor:=3, gens:=[1]), symbol:="m=3; 2: (0,2,+)(1,1,-); 3: (0,1,+)(2,2,-)", mass:=1/72, auts:=[108, 216]),
  rec(name:="Q(sqrt-7) 1,14", field:=rec(conductor:=7, gens:=[2]), symbol:="m=2; 2: (0,1,+)(1,1,+); 7: (0,1,+)(2,1,+)", mass:=3/2, auts:=[2, 4, 4, 4, 4]),
  rec(name:="Q(sqrt-11) 2,5,1", field:=rec(conductor:=11, gens:=[3]), symbol:="m=3; 2: (0,2,+)(1,1,-); 5: (0,2,+)(1,1,+); 11: (0,3,-)", mass:=93/16, auts:=[2, 2, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 12, 16, 24]),
  rec(name:="Q(zeta5) 1,2", field:=rec(conductor:=5, gens:=[1]), symbol:="m=2; 2: (0,1,+)(1,1,-); 5: (0,2,-)", mass:=1/100, auts:=[100]),
  rec(name:="Q(zeta5) 1,11", field:=rec(conductor:=5, gens:=[1]), symbol:="m=2; 2: (0,2,+); 5: (0,2,+); (11,[-17,0,26,26]): (0,1,+)(1,1,+); (11,[-7,0,6,6]): (0,1,+)(1,1,+)", mass:=18/25, auts:=[10, 10, 10, 10, 10, 10, 20, 20, 100, 100]),
  rec(name:="Q(zeta7) 1,3", field:=rec(conductor:=7, gens:=[1]), symbol:="m=2; 2: (0,2,+); 3: (0,1,+)(1,1,-); 7: (0,2,-)", mass:=13/147, auts:=[14, 84, 196]),
  rec(name:="Q(i) 1,2,3", field:=rec(conductor:=4, gens:=[1]), symbol:="m=3; 2: (0,2,-,0)(2,1,+,1); 3: (0,2,+)(1,1,-)", mass:=7/128, auts:=[32, 64, 128]),
  rec(name:="Q(i) 1,1,1,2", field:=rec(conductor:=4, gens:=[1]), symbol:="m=4; 2: (0,3,+,0)(2,1,+,1)", mass:=1/768, auts:=[1536, 1536]),
  rec(name:="Q(sqrt-2) 1,3", field:=rec(conductor:=8, gens:=[3]), symbol:="m=2; 3: (0,1,+)(1,1,+); 2: (0,2,+,0)", mass:=1/2, auts:=[4, 4]),
  rec(name:="Q(zeta8) 1,3", field:=rec(conductor:=8, gens:=[1]), symbol:="m=2; 3: (0,1,+)(1,1,+); 2: (0,2,+,0)", mass:=5/64, auts:=[16, 64]),
  rec(name:="Q(sqrt-2) 1,1,1 by x^2+2", field:=rec(poly:=[2, 0, 1], conj:=[0, -1]), symbol:="m=3; 2: (0,3,+,0)", mass:=1/32, auts:=[48, 96]),
  rec(name:="Q(sqrt-7) 1,14 by x^2-x+2", field:=rec(poly:=[2, -1, 1], conj:=[1, -1]), symbol:="m=2; 2: (0,1,+)(1,1,+); 7: (0,1,+)(2,1,+)", mass:=3/2, auts:=[2, 4, 4, 4, 4])
];
