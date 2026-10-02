Read("ErdahlTestFunctions.g");
Print("Beginning TestErdahl_full7\n");

# The perfect Delaunay polyhedra in dimension 7 (full space of functions):
# {0,1} x Z^6, 2_21 x Z, the Erdahl-Rybnikov 35-tope ER_7 and the Gosset
# polytope 3_21.
test1:=TestErdahl(rec(n:=7, space:="full", method:="recursive",
                      expected:=[[2,6], [27,1], [35,0], [56,0]]));
# The perfect Delaunay polytopes only, by flips from ER_7.
test2:=TestErdahl(rec(n:=7, space:="full", method:="polytopes",
                      expected:=[[35,0], [56,0]]));
ConcludeErdahl(test1 and test2);
