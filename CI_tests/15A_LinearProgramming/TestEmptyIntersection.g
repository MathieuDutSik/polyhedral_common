Read("../common.g");
Read("../access_points.g");
Print("Beginning TestEmptyIntersection\n");

# POLY_check_empty_intersection decides whether the convex hulls of two
# point sets (rows with first coordinate 1) are disjoint, by looking for
# a strictly separating hyperplane with a linear program. The expected
# answers are known geometrically.
Homog:=function(L)
    return List(L, x->Concatenation([1], x));
end;

TestCase:=function(name, L1, L2, expected)
    local result;
    result:=check_empty_intersection(Homog(L1), Homog(L2));
    if IsString(result) then
        Print(name, ": ", result, "\n");
        return false;
    fi;
    if result <> expected then
        Print(name, ": got ", result, ", expected ", expected, "\n");
        return false;
    fi;
    Print(name, ": ", result, " OK\n");
    return true;
end;

UnitSquare:=function(shift)
    return [[shift,0], [shift+1,0], [shift,1], [shift+1,1]];
end;
Cube:=Tuples([0,1], 3);
Simplex4:=function(shift)
    return [[shift,0,0,0], [shift+1,0,0,0], [shift,1,0,0], [shift,0,1,0],
            [shift,0,0,1]];
end;

ListCase:=[
  ["disjoint squares", UnitSquare(0), UnitSquare(3), true],
  ["overlapping squares", UnitSquare(0), UnitSquare(1/2), false],
  ["squares sharing an edge", UnitSquare(0), UnitSquare(1), false],
  ["crossing segments", [[0,0],[2,2]], [[0,2],[2,0]], false],
  ["parallel segments", [[0,0],[2,0]], [[0,1],[2,1]], true],
  ["cube and an outside point", Cube, [[2,2,2]], true],
  ["cube and its center", Cube, [[1/2,1/2,1/2]], false],
  ["cube and one of its vertices", Cube, [[1,1,1]], false],
  ["simplices apart in dimension 4", Simplex4(0), Simplex4(2), true],
  ["overlapping simplices in dimension 4", Simplex4(0), Simplex4(1/3), false]];

n_error:=0;
for eCase in ListCase
do
    if TestCase(eCase[1], eCase[2], eCase[3], eCase[4])=false then
        n_error:=n_error+1;
    fi;
od;

Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
