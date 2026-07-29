Read("../common.g");
Read("../access_points.g");
Print("Beginning TestLinearProgramming\n");

l_arith:=["mpq_class", "safe_rational", "cpp_rational", "mpq_rational"];

# The linear program is: minimize ineq[1] + Sum_j ineq[j+1] x_j over
# { x : FAC[i][1] + Sum_j FAC[i][j+1] x_j >= 0 }. The three scenarios
# and the certificates checked here, all in exact arithmetic:
# --- optimal: primal_solution and dual_solution are both set. The primal
#     point x is feasible and has objective value OptimalValue; the dual
#     multipliers lambda <= 0 satisfy ineq[j+1] + Sum_i lambda_i
#     FAC[i][j+1] = 0 and ineq[1] + Sum_i lambda_i FAC[i][1] =
#     OptimalValue, which proves the optimality of x; face is the set of
#     rows tight at x.
# --- infeasible: only dual_solution is set. It is a Farkas certificate:
#     lambda <= 0 with Sum_i lambda_i FAC[i][j+1] = 0 and
#     Sum_i lambda_i FAC[i][1] > 0, so Sum_i (-lambda_i) f_i(x) is a
#     negative constant and no feasible x exists.
# --- unbounded: only primal_solution is set. It is a ray d with
#     FAC[i]{[2..]} . d >= 0 for every row and ineq{[2..]} . d < 0, so
#     the objective decreases without bound along d.

eval_row:=function(row, x)
    local n, val, j;
    n:=Length(row);
    val:=row[1];
    for j in [2..n]
    do
        val:=val + row[j] * x[j-1];
    od;
    return val;
end;

TestOptimal:=function(eFile, eRec, TheLP, expected_value)
    local x, lambda, n, m, j, i, val, sum, face_check;
    if not IsBound(TheLP.primal_solution) or not IsBound(TheLP.dual_solution) then
        Print("eFile=", eFile, ": expected the optimal scenario with both primal and dual solutions\n");
        return false;
    fi;
    x:=TheLP.primal_solution;
    lambda:=TheLP.dual_solution;
    n:=Length(eRec.ineq);
    m:=Length(eRec.FAC);
    # The primal point: feasible and of the announced value.
    for i in [1..m]
    do
        if eval_row(eRec.FAC[i], x) < 0 then
            Print("eFile=", eFile, ": the primal solution violates row ", i, "\n");
            return false;
        fi;
    od;
    if eval_row(eRec.ineq, x) <> TheLP.OptimalValue then
        Print("eFile=", eFile, ": the primal objective value does not match OptimalValue\n");
        return false;
    fi;
    if TheLP.OptimalValue <> expected_value then
        Print("eFile=", eFile, ": OptimalValue=", TheLP.OptimalValue, " expected ", expected_value, "\n");
        return false;
    fi;
    # The dual multipliers: nonpositive, cancelling the objective gradient
    # and reproducing the optimal value, which proves optimality.
    for i in [1..m]
    do
        if lambda[i] > 0 then
            Print("eFile=", eFile, ": positive dual multiplier at row ", i, "\n");
            return false;
        fi;
    od;
    for j in [2..n]
    do
        sum:=eRec.ineq[j];
        for i in [1..m]
        do
            sum:=sum + lambda[i] * eRec.FAC[i][j];
        od;
        if sum <> 0 then
            Print("eFile=", eFile, ": the dual multipliers do not cancel the objective at column ", j, "\n");
            return false;
        fi;
    od;
    sum:=eRec.ineq[1];
    for i in [1..m]
    do
        sum:=sum + lambda[i] * eRec.FAC[i][1];
    od;
    if sum <> TheLP.OptimalValue then
        Print("eFile=", eFile, ": the dual value does not match OptimalValue\n");
        return false;
    fi;
    # The face: exactly the rows tight at the primal solution.
    face_check:=Filtered([1..m], i->eval_row(eRec.FAC[i], x)=0);
    if Set(TheLP.face) <> Set(face_check) then
        Print("eFile=", eFile, ": face=", TheLP.face, " but the tight rows are ", face_check, "\n");
        return false;
    fi;
    return true;
end;

TestInfeasible:=function(eFile, eRec, TheLP)
    local lambda, n, m, j, i, sum;
    if IsBound(TheLP.primal_solution) or not IsBound(TheLP.dual_solution) then
        Print("eFile=", eFile, ": expected the infeasible scenario with only a dual solution\n");
        return false;
    fi;
    lambda:=TheLP.dual_solution;
    n:=Length(eRec.ineq);
    m:=Length(eRec.FAC);
    for i in [1..m]
    do
        if lambda[i] > 0 then
            Print("eFile=", eFile, ": positive multiplier in the Farkas certificate at row ", i, "\n");
            return false;
        fi;
    od;
    for j in [2..n]
    do
        sum:=0;
        for i in [1..m]
        do
            sum:=sum + lambda[i] * eRec.FAC[i][j];
        od;
        if sum <> 0 then
            Print("eFile=", eFile, ": the Farkas certificate does not cancel column ", j, "\n");
            return false;
        fi;
    od;
    sum:=0;
    for i in [1..m]
    do
        sum:=sum + lambda[i] * eRec.FAC[i][1];
    od;
    if sum <= 0 then
        Print("eFile=", eFile, ": the Farkas certificate has nonpositive constant term\n");
        return false;
    fi;
    return true;
end;

TestUnbounded:=function(eFile, eRec, TheLP)
    local d, n, m, j, i, sum;
    if not IsBound(TheLP.primal_solution) or IsBound(TheLP.dual_solution) then
        Print("eFile=", eFile, ": expected the unbounded scenario with only a primal ray\n");
        return false;
    fi;
    d:=TheLP.primal_solution;
    n:=Length(eRec.ineq);
    m:=Length(eRec.FAC);
    for i in [1..m]
    do
        sum:=0;
        for j in [2..n]
        do
            sum:=sum + eRec.FAC[i][j] * d[j-1];
        od;
        if sum < 0 then
            Print("eFile=", eFile, ": the ray exits the feasible cone at row ", i, "\n");
            return false;
        fi;
    od;
    sum:=0;
    for j in [2..n]
    do
        sum:=sum + eRec.ineq[j] * d[j-1];
    od;
    if sum >= 0 then
        Print("eFile=", eFile, ": the objective does not decrease along the ray\n");
        return false;
    fi;
    return true;
end;

TestLinearProgram:=function(eFile)
    local eRec, arith, options, TheLP, test;
    eRec:=ReadAsFunction(eFile)();
    if eRec=fail then
        Print("Could not parse input file ", eFile, "\n");
        return false;
    fi;
    for arith in l_arith
    do
        options:=rec(print_info:=true, arith:=arith);
        TheLP:=get_linear_programming(eRec.FAC, eRec.ineq, options);
        if IsString(TheLP) then
            Print("eFile=", eFile, ": ", TheLP, "\n");
            return false;
        fi;
        if eRec.scenario="optimal" then
            test:=TestOptimal(eFile, eRec, TheLP, eRec.optimal_value);
        elif eRec.scenario="infeasible" then
            test:=TestInfeasible(eFile, eRec, TheLP);
        elif eRec.scenario="unbounded" then
            test:=TestUnbounded(eFile, eRec, TheLP);
        else
            Print("eFile=", eFile, ": unknown scenario ", eRec.scenario, "\n");
            return false;
        fi;
        if test=false then
            Print("eFile=", eFile, ": failure for arith=", arith, "\n");
            return false;
        fi;
    od;
    Print("eFile=", eFile, ": scenario=", eRec.scenario, " OK\n");
    return true;
end;

ListFiles:=["Example1_cube_optimal.g",
            "Example2_infeasible.g",
            "Example3_unbounded.g",
            "Example4_cross_degenerate.g"];

n_error:=0;
for eFile in ListFiles
do
    if TestLinearProgram(eFile)=false then
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
