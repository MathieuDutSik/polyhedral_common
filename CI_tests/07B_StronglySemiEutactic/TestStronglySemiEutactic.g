Read("../common.g");
Read("../access_points.g");

ListCases:=ReadAsFunction("ListCases.g")();

# Pairing the eutaxy relation A^{-1} = mu * sum_{v in S} v v^T with A
# itself gives n = <A^{-1}, A> = mu * sum_{v in S} v^T A v = mu * |S| * min.
# Any certificate returned by the program has to satisfy it.
CheckTraceIdentity:=function(n, eResult)
    if eResult.mu * eResult.siz * eResult.min <> n then
        Print("The certificate violates mu * |S| * min = n\n");
        Print("mu=", eResult.mu, " |S|=", eResult.siz, " min=", eResult.min, " n=", n, "\n");
        return false;
    fi;
    return true;
end;

TestOneCase:=function(eCase, arith)
    local n, options, eResult;
    n:=Length(eCase.gram);
    options:=rec(arith:=arith, print_info:=true);
    if IsBound(eCase.max_node) then
        options.max_node:=eCase.max_node;
    fi;
    eResult:=test_strongly_semi_eutactic(eCase.gram, options);
    if is_error(eResult) then
        return false;
    fi;
    if eResult.min <> eCase.min then
        Print("min=", eResult.min, " but the expected value is ", eCase.min, "\n");
        return false;
    fi;
    if eResult.n_shv <> eCase.n_shv then
        Print("|SHV|=", eResult.n_shv, " but the expected value is ", eCase.n_shv, "\n");
        return false;
    fi;
    if eResult.resolved <> eCase.resolved then
        Print("resolved=", eResult.resolved, " but the expected value is ", eCase.resolved, "\n");
        return false;
    fi;
    if eResult.resolved = false then
        # The node budget was exhausted, there is nothing more to check.
        Print("  UNRESOLVED as expected, n_node=", eResult.n_node, "\n");
        return true;
    fi;
    if eResult.is_sse <> eCase.is_sse then
        Print("is_sse=", eResult.is_sse, " but the expected value is ", eCase.is_sse, "\n");
        return false;
    fi;
    if eResult.is_sse = false then
        Print("  not strongly semi-eutactic as expected, n_node=", eResult.n_node, "\n");
        return true;
    fi;
    if eResult.mu <> eCase.mu then
        Print("mu=", eResult.mu, " but the expected value is ", eCase.mu, "\n");
        return false;
    fi;
    if eResult.siz <> eCase.siz then
        Print("|S|=", eResult.siz, " but the expected value is ", eCase.siz, "\n");
        return false;
    fi;
    if Length(eResult.subset) <> eCase.siz then
        Print("|S|=", eCase.siz, " but the printed subset has ", Length(eResult.subset), " entries\n");
        return false;
    fi;
    if Length(Set(eResult.subset)) <> eCase.siz then
        Print("The printed subset has repeated entries\n");
        return false;
    fi;
    if First(eResult.subset, x->x<0 or x>=eResult.n_shv) <> fail then
        Print("The printed subset has an entry outside of [0, ", eResult.n_shv, "[\n");
        return false;
    fi;
    if CheckTraceIdentity(n, eResult) = false then
        return false;
    fi;
    Print("  strongly semi-eutactic as expected, mu=", eResult.mu, " |S|=", eResult.siz, " n_node=", eResult.n_node, "\n");
    return true;
end;

# Every case is run with the gmp arithmetic. The cases marked "sweep"
# are run with the other arithmetics as well, which must not change
# the answer.
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
            if arith = "gmp" or (IsBound(eCase.sweep) and eCase.sweep) then
                test:=TestOneCase(eCase, arith);
                if test = false then
                    Print("Failure for name=", eCase.name, " arith=", arith, "\n");
                    return false;
                fi;
            fi;
        od;
    od;
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
