# The cube [0,1]^3: x_i >= 0 and 1 - x_i >= 0.
# Minimize 2 - x - 2y - 3z: unique optimal vertex (1,1,1), value -4,
# tight on the three rows 1 - x_i >= 0.
return rec(FAC:=[[0,  1,  0,  0],
                 [0,  0,  1,  0],
                 [0,  0,  0,  1],
                 [1, -1,  0,  0],
                 [1,  0, -1,  0],
                 [1,  0,  0, -1]],
           ineq:=[2, -1, -2, -3],
           scenario:="optimal",
           optimal_value:=-4);
