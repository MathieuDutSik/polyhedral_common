# Infeasible: x >= 1 and x <= 0 together with y >= 0.
# The dual solution is a Farkas certificate supported on the two
# contradicting rows.
return rec(FAC:=[[-1,  1,  0],
                 [ 0, -1,  0],
                 [ 0,  0,  1]],
           ineq:=[0, 1, 0],
           scenario:="infeasible");
