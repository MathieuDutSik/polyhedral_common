Read("../common.g");
Read("../access_points.g");
Print("Beginning Test generalizedpolytope difference\n");

Test_difference:=function(EXT)
    local GRP, triang1, eElt, triang2, get_gen_polytope, gp1, gp2, result;
    GRP:=get_grp_automorphy(EXT).GAPperm;
    triang1:=get_triangulation_of_polytope(EXT);
    eElt:=Random(GRP);
    triang2:=List(triang1, x->OnSets(x, eElt));
    get_gen_polytope:=function(the_triang)
        local ListEXT, trig;
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


Test_vertices:=function(EXT)
    local ListFacet, ListEXT, eIso, eFacet, NewEXT, LVert;
    ListFacet:=get_dual_desc_incidence(EXT, "lrs");
    ListEXT:=[];
    eIso:=Sum(EXT) / Length(EXT);
    for eFacet in ListFacet
    do
        NewEXT:=Concatenation([eIso], EXT{eFacet});
        Add(ListEXT, NewEXT);
    od;
    LVert:=get_polygen_vertices(ListEXT);
    if is_error(LVert) then
        return false;
    fi;
    if Set(LVert)<>Set(EXT) then
        return false;
    fi;
    return true;
end;





TreatOneExample:=function(EXT)
    local test1, test2;
    test1:=Test_difference(EXT);
    if test1=false then
        return false;
    fi;
    test2:=Test_vertices(EXT);
    if test2=false then
        return false;
    fi;
    return true;
end;



eCase1:=rec(EXT:=SpecialPolytopes("24cell"), the_volume:=8, name:="24cell");
eCase2:=rec(EXT:=Hypercube(3), the_volume:=1, name:="H3");
ListCase:=[eCase1, eCase2];


ProcessAllCases:=function()
    local n_error, ListCaseError, i_case, eCase, test;
    n_error:=0;
    ListCaseError:=[];
    i_case:=0;
    for eCase in ListCase
    do
        i_case:=i_case + 1;
        Print("TestDifference, i_case=", i_case, " name=", eCase.name, "\n");
        test:=TreatOneExample(eCase.EXT);
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

