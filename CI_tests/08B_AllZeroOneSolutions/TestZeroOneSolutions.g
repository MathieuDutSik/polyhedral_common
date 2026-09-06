Read("../common.g");
Read("../access_points.g");

# The small cases are checked against a brute force enumeration of the
# 2^n vectors done here in GAP, so their expected answer does not come
# from the program being tested.
BruteForceZeroOne:=function(A, b)
    local n, m, ListSol, val, x, i, j, sum, is_sol;
    n:=Length(A[1]);
    m:=Length(A);
    ListSol:=[];
    for val in [0..2^n-1]
    do
        x:=List([1..n], j->QuoInt(val, 2^(j-1)) mod 2);
        is_sol:=true;
        for i in [1..m]
        do
            sum:=0;
            for j in [1..n]
            do
                sum:=sum + A[i][j] * x[j];
            od;
            if sum <> b[i] then
                is_sol:=false;
                break;
            fi;
        od;
        if is_sol then
            Add(ListSol, x);
        fi;
    od;
    return ListSol;
end;

CheckSolutions:=function(A, b, eResult, ListExpected)
    local n, i, j, sum, x;
    n:=Length(A[1]);
    if Length(eResult.solutions) <> eResult.n_solution then
        Print("The program reports ", eResult.n_solution, " solutions but ",
              Length(eResult.solutions), " were written\n");
        return false;
    fi;
    # Each returned vector has to be 0/1 and to satisfy the system
    for x in eResult.solutions
    do
        if Length(x) <> n then
            Print("A returned solution has length ", Length(x), " instead of ", n, "\n");
            return false;
        fi;
        if Length(Filtered(x, y->y<>0 and y<>1)) > 0 then
            Print("A returned solution is not 0/1\n");
            return false;
        fi;
        for i in [1..Length(A)]
        do
            sum:=0;
            for j in [1..n]
            do
                sum:=sum + A[i][j] * x[j];
            od;
            if sum <> b[i] then
                Print("A returned solution violates the row ", i, "\n");
                return false;
            fi;
        od;
    od;
    if Set(eResult.solutions) <> Set(ListExpected) then
        Print("The set of solutions differs from the expected one\n");
        Print("|returned|=", Length(Set(eResult.solutions)),
              " |expected|=", Length(Set(ListExpected)), "\n");
        return false;
    fi;
    return true;
end;

TestSmallCase:=function(eCase, arith)
    local options, eResult, ListExpected;
    options:=rec(arith:=arith, print_info:=true);
    if IsBound(eCase.lp_max_depth) then
        options.lp_max_depth:=eCase.lp_max_depth;
    fi;
    eResult:=get_zero_one_solutions(eCase.A, eCase.b, options);
    if is_error(eResult) then
        return false;
    fi;
    if eResult.resolved <> true then
        Print("The enumeration did not resolve\n");
        return false;
    fi;
    ListExpected:=BruteForceZeroOne(eCase.A, eCase.b);
    if Length(ListExpected) <> eCase.n_solution then
        Print("The brute force found ", Length(ListExpected),
              " solutions but the case declares ", eCase.n_solution, "\n");
        return false;
    fi;
    if CheckSolutions(eCase.A, eCase.b, eResult, ListExpected) = false then
        return false;
    fi;
    Print("  ", eResult.n_solution, " solutions, matching the brute force, n_node=",
          eResult.n_node, "\n");
    return true;
end;

# The instance Problem1 of the section by the branch and bound, whose
# solutions are the recorded ones. Problem2 is left out here: the
# branch and bound does not close its search tree, see README.md. It
# is covered below by the lattice enumeration.
TestProblem1:=function()
    local A, b, ListExpected, eResult;
    A:=ReadMatrixFile("Problem1.matrix");
    b:=ReadVectorFile("Problem1.rhs");
    ListExpected:=ReadMatrixFile("Problem1.solutions");
    Print("Problem1: A is ", Length(A), "x", Length(A[1]), " with ",
          Length(ListExpected), " recorded solutions\n");
    eResult:=get_zero_one_solutions(A, b, rec(print_info:=true));
    if is_error(eResult) then
        return false;
    fi;
    if eResult.resolved <> true then
        Print("The enumeration of Problem1 did not resolve\n");
        return false;
    fi;
    if CheckSolutions(A, b, eResult, ListExpected) = false then
        return false;
    fi;
    Print("  ", eResult.n_solution, " solutions, matching Problem1.solutions, n_node=",
          eResult.n_node, "\n");
    return true;
end;

# The node budget has to stop the enumeration of Problem1, which needs
# more than five million nodes.
TestNodeBudget:=function()
    local A, b, eResult;
    A:=ReadMatrixFile("Problem1.matrix");
    b:=ReadVectorFile("Problem1.rhs");
    eResult:=get_zero_one_solutions(A, b, rec(max_node:=1000, print_info:=true));
    if is_error(eResult) then
        return false;
    fi;
    if eResult.resolved <> false then
        Print("The enumeration should have been stopped by the node budget\n");
        return false;
    fi;
    Print("  UNRESOLVED as expected, n_node=", eResult.n_node, "\n");
    return true;
end;

# The same small cases through the pruned lattice enumeration of
# MILP_ZeroOneLattice, which is an independent implementation: it must
# agree with the brute force too.
TestSmallCaseLattice:=function(eCase)
    local eResult, ListExpected;
    eResult:=get_zero_one_lattice_solutions(eCase.A, eCase.b, rec(print_info:=true));
    if is_error(eResult) then
        return false;
    fi;
    if eResult.resolved <> true then
        Print("The lattice enumeration did not resolve\n");
        return false;
    fi;
    ListExpected:=BruteForceZeroOne(eCase.A, eCase.b);
    if CheckSolutions(eCase.A, eCase.b, eResult, ListExpected) = false then
        return false;
    fi;
    Print("  lattice: ", eResult.n_solution, " solutions, matching the brute force, n_node=",
          eResult.n_node, "\n");
    return true;
end;

# Problem1 and Problem2 of the section, by the lattice enumeration.
# Problem2 is out of reach of the branch and bound, its rows being sums
# of about ten small coefficients equal to 10, so the bound propagation
# has too much slack; the lattice does it in about five minutes.
TestProblemLattice:=function(idx)
    local A, b, ListExpected, eResult;
    A:=ReadMatrixFile(Concatenation("Problem", String(idx), ".matrix"));
    b:=ReadVectorFile(Concatenation("Problem", String(idx), ".rhs"));
    ListExpected:=ReadMatrixFile(Concatenation("Problem", String(idx), ".solutions"));
    Print("Problem", idx, ": A is ", Length(A), "x", Length(A[1]), " with ",
          Length(ListExpected), " recorded solutions\n");
    eResult:=get_zero_one_lattice_solutions(A, b, rec(print_info:=true));
    if is_error(eResult) then
        return false;
    fi;
    if eResult.resolved <> true then
        Print("The lattice enumeration of Problem", idx, " did not resolve\n");
        return false;
    fi;
    if CheckSolutions(A, b, eResult, ListExpected) = false then
        return false;
    fi;
    Print("  ", eResult.n_solution, " solutions, matching Problem", idx,
          ".solutions, n_node=", eResult.n_node, "\n");
    return true;
end;

# G_mat by brute force over Sym(n): the permutations of the columns
# for which the multiset of the rows, right hand side included, is
# unchanged. Only for the small cases.
BruteForceSystemSymmetry:=function(A, b)
    local n, m, ref, ListPerm, g, cur, i;
    n:=Length(A[1]);
    m:=Length(A);
    ref:=SortedList(List([1..m], i->[A[i], b[i]]));
    ListPerm:=[];
    for g in SymmetricGroup(n)
    do
        cur:=SortedList(List([1..m], i->[Permuted(A[i], g), b[i]]));
        if cur = ref then
            Add(ListPerm, g);
        fi;
    od;
    return Group(ListPerm);
end;

# G_aff by brute force: the permutations of the columns preserving the
# row space of [A | -b], which the row reduced echelon form decides.
BruteForceAffineSymmetry:=function(A, b)
    local n, m, ref, ListPerm, g, cur, i;
    n:=Length(A[1]);
    m:=Length(A);
    ref:=TriangulizedMat(List([1..m], i->Concatenation(A[i], [-b[i]])));
    ListPerm:=[];
    for g in SymmetricGroup(n)
    do
        cur:=TriangulizedMat(List([1..m], i->Concatenation(Permuted(A[i], g), [-b[i]])));
        if cur = ref then
            Add(ListPerm, g);
        fi;
    od;
    return Group(ListPerm);
end;

# The two groups of a case. G_mat is contained in G_aff always, G_mat
# depending on the rows that were written down and G_aff not. For a
# small enough case both are compared with the brute force above.
MAX_DEGREE_BRUTE_FORCE_SYMMETRY:=7;

TestSymmetry:=function(A, b, name)
    local n, resMat, resAff, Gmat, Gaff, gen, Gmat_bf, Gaff_bf;
    n:=Length(A[1]);
    resMat:=get_system_symmetry(A, b, rec(print_info:=true));
    if is_error(resMat) then
        return false;
    fi;
    resAff:=get_affine_symmetry(A, b, rec(print_info:=true));
    if is_error(resAff) then
        return false;
    fi;
    Gmat:=resMat.group;
    Gaff:=resAff.group;
    # The order that the program reports has to be the one GAP finds
    if Size(Gmat) <> resMat.size_reported then
        Print("For ", name, " MILP_SystemSymmetry reports the order ",
              resMat.size_reported, " but GAP finds ", Size(Gmat), "\n");
        return false;
    fi;
    if Size(Gaff) <> resAff.size_reported then
        Print("For ", name, " MILP_AffineSymmetry reports the order ",
              resAff.size_reported, " but GAP finds ", Size(Gaff), "\n");
        return false;
    fi;
    # G_mat is a subgroup of G_aff
    for gen in GeneratorsOfGroup(Gmat)
    do
        if not gen in Gaff then
            Print("For ", name, " a generator of G_mat is not in G_aff\n");
            return false;
        fi;
    od;
    if n <= MAX_DEGREE_BRUTE_FORCE_SYMMETRY then
        Gmat_bf:=BruteForceSystemSymmetry(A, b);
        Gaff_bf:=BruteForceAffineSymmetry(A, b);
        if Size(Gmat) <> Size(Gmat_bf) then
            Print("For ", name, " |G_mat|=", Size(Gmat),
                  " but the brute force finds ", Size(Gmat_bf), "\n");
            return false;
        fi;
        if Size(Gaff) <> Size(Gaff_bf) then
            Print("For ", name, " |G_aff|=", Size(Gaff),
                  " but the brute force finds ", Size(Gaff_bf), "\n");
            return false;
        fi;
        Print("  |G_mat|=", Size(Gmat), " |G_aff|=", Size(Gaff),
              ", both matching the brute force\n");
    else
        Print("  |G_mat|=", Size(Gmat), " |G_aff|=", Size(Gaff),
              ", G_mat contained in G_aff\n");
    fi;
    return true;
end;

# The groups of Problem1 and Problem2, and the check that the recorded
# solutions are a union of orbits: the group and the enumeration are
# independent pieces of code, so their agreement is worth something.
TestProblemSymmetry:=function(idx, size_expected)
    local A, b, sols, sets, resMat, resAff, Gmat, Gaff, s, g, img;
    A:=ReadMatrixFile(Concatenation("Problem", String(idx), ".matrix"));
    b:=ReadVectorFile(Concatenation("Problem", String(idx), ".rhs"));
    resMat:=get_system_symmetry(A, b, rec(print_info:=true));
    if is_error(resMat) then
        return false;
    fi;
    resAff:=get_affine_symmetry(A, b, rec(print_info:=true));
    if is_error(resAff) then
        return false;
    fi;
    Gmat:=resMat.group;
    Gaff:=resAff.group;
    if resMat.size_reported <> size_expected then
        Print("Problem", idx, ": |G_mat|=", resMat.size_reported,
              " but the expected order is ", size_expected, "\n");
        return false;
    fi;
    if resAff.size_reported <> size_expected then
        Print("Problem", idx, ": |G_aff|=", resAff.size_reported,
              " but the expected order is ", size_expected, "\n");
        return false;
    fi;
    sols:=ReadMatrixFile(Concatenation("Problem", String(idx), ".solutions"));
    sets:=Set(List(sols, r->Filtered([1..Length(r)], j->r[j]=1)));
    for s in sets
    do
        for g in Concatenation(GeneratorsOfGroup(Gmat), GeneratorsOfGroup(Gaff))
        do
            img:=OnSets(s, g);
            if not img in sets then
                Print("Problem", idx, ": the recorded solutions are not stable ",
                      "under the group\n");
                return false;
            fi;
        od;
    od;
    Print("  |G_mat|=|G_aff|=", size_expected,
          ", and the ", Length(sets), " recorded solutions are a union of orbits\n");
    return true;
end;

ListCases:=[];
# One row, all the subsets of a given size
Add(ListCases, rec(name:="choose_2_of_4", A:=[[1,1,1,1]], b:=[2], n_solution:=6));
# A single solution
Add(ListCases, rec(name:="single_solution", A:=[[1,2,3],[1,1,1]], b:=[4,2], n_solution:=1));
# No solution at all, by parity
Add(ListCases, rec(name:="no_solution_parity", A:=[[2,2,2]], b:=[3], n_solution:=0));
# Negative coefficients, so the lower bound of a row is not zero
Add(ListCases, rec(name:="negative_coefficients", A:=[[1,-1,2,-2],[1,1,1,1]], b:=[0,2], n_solution:=2));
# A denser case, still small enough for the brute force
Add(ListCases, rec(name:="dense_12", A:=[[3,1,4,1,5,9,2,6,5,3,5,8],[1,1,1,1,1,1,1,1,1,1,1,1],[2,0,2,0,2,0,2,0,2,0,2,0]], b:=[20,5,6], n_solution:=22));
# The same, with the linear programming pruning pushed deeper
Add(ListCases, rec(name:="dense_12_lp_deep", A:=[[3,1,4,1,5,9,2,6,5,3,5,8],[1,1,1,1,1,1,1,1,1,1,1,1],[2,0,2,0,2,0,2,0,2,0,2,0]], b:=[20,5,6], n_solution:=22, lp_max_depth:=20));
# The same, with the linear programming pruning disabled
Add(ListCases, rec(name:="dense_12_no_lp", A:=[[3,1,4,1,5,9,2,6,5,3,5,8],[1,1,1,1,1,1,1,1,1,1,1,1],[2,0,2,0,2,0,2,0,2,0,2,0]], b:=[20,5,6], n_solution:=22, lp_max_depth:=-1));

ListArithmetic:=["gmp", "gmp_boost", "multi_boost"];

FullTest:=function()
    local iCase, eCase, arith, test;
    iCase:=0;
    for eCase in ListCases
    do
        iCase:=iCase + 1;
        Print("iCase=", iCase, " / ", Length(ListCases), " name=", eCase.name, "\n");
        for arith in ListArithmetic
        do
            test:=TestSmallCase(eCase, arith);
            if test = false then
                Print("Failure for name=", eCase.name, " arith=", arith, "\n");
                return false;
            fi;
        od;
        if TestSmallCaseLattice(eCase) = false then
            Print("Failure of the lattice enumeration for name=", eCase.name, "\n");
            return false;
        fi;
    od;
    Print("Now the node budget\n");
    if TestNodeBudget() = false then
        return false;
    fi;
    Print("Now Problem1 by the branch and bound\n");
    if TestProblem1() = false then
        return false;
    fi;
    Print("Now Problem1 by the lattice\n");
    if TestProblemLattice(1) = false then
        return false;
    fi;
    Print("Now Problem2 by the lattice\n");
    if TestProblemLattice(2) = false then
        return false;
    fi;
    Print("Now the symmetry groups of the small cases\n");
    iCase:=0;
    for eCase in ListCases
    do
        iCase:=iCase + 1;
        Print("iCase=", iCase, " / ", Length(ListCases), " name=", eCase.name, "\n");
        if TestSymmetry(eCase.A, eCase.b, eCase.name) = false then
            return false;
        fi;
    od;
    Print("Now the symmetry groups of Problem1 and Problem2\n");
    if TestProblemSymmetry(1, 474989023199232) = false then
        return false;
    fi;
    if TestProblemSymmetry(2, 6) = false then
        return false;
    fi;
    return true;
end;

result:=FullTest();
Print("result=", result, "\n");

CI_Decision_Reset();
if result = false then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
