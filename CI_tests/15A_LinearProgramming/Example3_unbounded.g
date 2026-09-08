# Unbounded: minimize -x - y over the positive quadrant. The primal
# solution is a ray d of the quadrant with negative objective.
return rec(FAC:=[[0, 1, 0],
                 [0, 0, 1]],
           ineq:=[0, -1, -1],
           scenario:="unbounded");
