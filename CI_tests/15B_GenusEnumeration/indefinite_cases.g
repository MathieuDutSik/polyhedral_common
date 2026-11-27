# The cases of TestIndefiniteRepresentatives.g: indefinite lattices over Z
# of rank >= 3 and the number of classes of their genus (one per improper
# spinor genus, Eichler). The genera with 2 classes are those found by
# Hecke.jl (integer_genera, signatures (2,1) and (1,2), |det| <= 400; none
# of rank 4 or 5 there). The genera with 4 and 8 classes are built with
# primes p = 1 mod 8 that are squares modulo each other (17, 89, 409),
# so that each keeps a spinor bit. For all of them the number of classes is
# |A / K| of SageMath's improper spinor kernel, and for T, Y, the rank 4
# case with 4 classes and T doubled the representatives of
# GENUS_IndefiniteRepresentatives were checked with Hecke.jl: in the genus,
# pairwise non-isometric. Hecke.jl (v0.40.1): spinor_genera_in_genus
# returns 2^(C-1) lattices for C classes (one prime per element of the
# group is taken as a basis: 8 for T, 128 for Y, isometric in pairs by
# its own is_isometric). SageMath (10.9): Genus.spinor_generators returns
# a single prime (U = Q.subgroup((g,) + Q.gens()) is the whole group after
# the first one), so Genus.representatives() gives 2 lattices for T and for
# Y. Both are right only up to 2 classes.
return [
  rec(name:="Hecke.jl, signature (2,1), det -128", gram:=[[63, 64, -2], [64, 64, 0], [-2, 0, -2]], n_class:=2),
  rec(name:="Hecke.jl, signature (2,1), det -250", gram:=[[88, 95, -5], [95, 100, 0], [-5, 0, -10]], n_class:=2),
  rec(name:="Hecke.jl, signature (2,1), det -256", gram:=[[-930, 0, -976], [0, 1, 0], [-976, 0, -1024]], n_class:=2),
  rec(name:="Hecke.jl, signature (1,2), det 128", gram:=[[-63, -64, 0], [-64, -64, 0], [0, 0, -2]], n_class:=2),
  rec(name:="Hecke.jl, signature (1,2), det 250", gram:=[[-2, -15, 0], [-15, -100, 0], [0, 0, -10]], n_class:=2),
  rec(name:="Hecke.jl, signature (1,2), det 256", gram:=[[-225, -226, -240], [-226, -226, -240], [-240, -240, -256]], n_class:=2),
  rec(name:="T = diag(1, -P^2, P^4), P = 17*89", gram:=[[1, 0, 0], [0, -2289169, 0], [0, 0, 5240294710561]], n_class:=4),
  rec(name:="Y = diag(1, -P^2, P^4), P = 17*89*409", gram:=[[1, 0, 0], [0, -382934479489, 0], [0, 0, 146638815581511361701121]], n_class:=8),
  rec(name:="diag(1, -17^2, 17^4, -17^6)", gram:=[[1, 0, 0, 0], [0, -289, 0, 0], [0, 0, 83521, 0], [0, 0, 0, -24137569]], n_class:=2),
  rec(name:="diag(1, -P^2, P^4, -P^6), P = 17*89", gram:=[[1, 0, 0, 0], [0, -2289169, 0, 0], [0, 0, 5240294710561, 0], [0, 0, 0, -11995920202280213809]], n_class:=4),
  rec(name:="T doubled, even", gram:=[[2, 0, 0], [0, -4578338, 0], [0, 0, 10480589421122]], n_class:=4),
  rec(name:="I_{2,1}", gram:=[[1, 0, 0], [0, 1, 0], [0, 0, -1]], n_class:=1),
  rec(name:="I_{2,2}", gram:=[[1, 0, 0, 0], [0, -1, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]], n_class:=1),
  rec(name:="diag(1, -P^2, P^4, -1, 1)", gram:=[[1, 0, 0, 0, 0], [0, -2289169, 0, 0, 0], [0, 0, 5240294710561, 0, 0], [0, 0, 0, -1, 0], [0, 0, 0, 0, 1]], n_class:=1),
  rec(name:="II_{9,1} = E8 + U", gram:=[[2, -1, 0, 0, 0, 0, 0, 0, 0, 0], [-1, 2, -1, 0, 0, 0, 0, 0, 0, 0], [0, -1, 2, -1, 0, 0, 0, -1, 0, 0], [0, 0, -1, 2, -1, 0, 0, 0, 0, 0], [0, 0, 0, -1, 2, -1, 0, 0, 0, 0], [0, 0, 0, 0, -1, 2, -1, 0, 0, 0], [0, 0, 0, 0, 0, -1, 2, 0, 0, 0], [0, 0, -1, 0, 0, 0, 0, 2, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0]], n_class:=1)
];
