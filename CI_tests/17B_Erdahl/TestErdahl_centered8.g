Read("ErdahlTestFunctions.g");
Print("Beginning TestErdahl_centered8\n");

# The centrally symmetric perfect Delaunay polyhedra in dimension 8 of center
# e_1/2, i.e. those of the space of functions f with f(e_1 - x) = f(x):
# {0,1} x Z^7, 3_21 x Z and the 72-vertex polytope D_2^8 of Dutour, Erdahl and
# Rybnikov, the only centrally symmetric one of their 27 perfect 8-polytopes.
test1:=TestErdahl(rec(n:=8, space:="centered", FileCenter:="center8.txt",
                      method:="recursive",
                      expected:=[[2,7], [56,1], [72,0]]));
# The perfect Delaunay polytopes only, by flips from the symmetrization of
# ER_7, which is D_2^8.
test2:=TestErdahl(rec(n:=8, space:="centered", FileCenter:="center8.txt",
                      method:="polytopes", expected:=[[72,0]]));
ConcludeErdahl(test1 and test2);
