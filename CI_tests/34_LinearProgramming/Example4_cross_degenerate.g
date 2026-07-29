# The 3-dimensional cross polytope |x| + |y| + |z| <= 1 written with
# its 8 facets. Minimize x: the optimal vertex (-1,0,0) is degenerate,
# tight on 4 of the 8 rows, which exercises the tie-breaking of the
# simplex method. Value -1.
return rec(FAC:=[[1,  1,  1,  1],
                 [1,  1,  1, -1],
                 [1,  1, -1,  1],
                 [1,  1, -1, -1],
                 [1, -1,  1,  1],
                 [1, -1,  1, -1],
                 [1, -1, -1,  1],
                 [1, -1, -1, -1]],
           ineq:=[0, 1, 0, 0],
           scenario:="optimal",
           optimal_value:=-1);
