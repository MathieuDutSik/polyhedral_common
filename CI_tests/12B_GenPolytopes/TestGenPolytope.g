Read("../common.g");
Read("../access_points.g");
Print("Beginning Test copositivity\n");






TreatOneExample:=function(eCase)
    local GRP;
    GRP:=get_grp_automorphy(EXT);
    triang1:=get_triangulation_of_polytope(EXT);
    eElt:=Random(GRP);
    triang2:=List(triang1, x->OnSets(x, eElt));
    get_gen_polytope:=function(EXT, the_triang)
        local ListEXT;
        ListEXT:=[];
        for trig in the_triang
        do
            Add(ListEXT, EXT{trig});
        od;
        return ListEXT;
    end;
    gp1:=get_gen_polytope(triang1);
    gp2:=get_gen_polytope(triang2);
    result:=get_polygen_difference(gp1, gp2);
    if is_error(result) then
        return false;
    fi;
    if result.vol_d12<>0 then
        return false;
    fi;
    if result.vol_d21<>0 then
        return false;
    fi;
    return true;
end;

ListCase:=[];


ProcessAllCases:=function()
    local n_error, ListCaseError, i_case, eCase, test;
    n_error:=0;
    ListCaseError:=[];
    i_case:=0;
    for eCase in ListCase
    do
        i_case:=i_case + 1;
        Print("TestCopositivity, i_case=", i_case, " name=", eCase.name, "\n");
        test:=TreatOneExample(eCase);
        if test=false then
            n_error:=n_error + 1;
            Add(ListCaseError, eCase);
        fi;
    od;
    Print("n_case=", Length(ListCase), " n_error=", n_error, "\n");
    return rec(n_error:=n_error, ListCaseError:=ListCaseError);
end;

RecResult:=ProcessAllCases();
Print("RecResult.n_error=", RecResult.n_error, "\n");

if RecResult.n_error > 0 then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

