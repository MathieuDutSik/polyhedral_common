Read("ErdahlTestFunctions.g");
Print("Beginning TestErdahl_full7\n");

# The perfect Delaunay polyhedra in dimension 7 (full space of functions):
# {0,1} x Z^6, 2_21 x Z, the Erdahl-Rybnikov 35-tope ER_7 and the Gosset
# polytope 3_21.
test:=TestErdahl(rec(n:=7, space_args:="full",
                     expected:=[[2,6], [27,1], [35,0], [56,0]]));
ConcludeErdahl(test);
