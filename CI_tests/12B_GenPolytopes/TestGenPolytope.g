Read("../common.g");
Read("../access_points.g");
Print("Beginning Test generalizedpolytope difference\n");

Test_difference:=function(EXT)
    local GRP, triang1, eElt, triang2, get_gen_polytope, gp1, gp2, result;
    GRP:=get_grp_automorphy(EXT).GAPperm;
    triang1:=get_triangulation_of_polytope(EXT);
    eElt:=Random(GRP);
    # A simplex is a set of vertices, but the triangulation comes out in
    # the order lrs produced it, which OnSets does not accept.
    triang2:=List(triang1, x->OnSets(Set(x), eElt));
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

# Cases with known answers: overlapping and non-convex unions, including
# dimension 1 (whose boundary consists of points).
Box:=function(lo, hi)
    return List(Cartesian(List([1..Length(lo)], i->[lo[i], hi[i]])), x->Concatenation([1], x));
end;

Test_difference_exact:=function(name, gp1, gp2, vol_int, vol_d12, vol_d21)
    local result;
    result:=get_polygen_difference(gp1, gp2);
    if is_error(result) then
        Print("Test_difference_exact, name=", name, " error\n");
        return false;
    fi;
    if result.vol_int<>vol_int or result.vol_d12<>vol_d12 or result.vol_d21<>vol_d21 then
        Print("Test_difference_exact, name=", name, " vol_int=", result.vol_int, " vol_d12=", result.vol_d12, " vol_d21=", result.vol_d21, "\n");
        return false;
    fi;
    return true;
end;

Test_vertices_exact:=function(name, gp, LVert)
    local result;
    result:=get_polygen_vertices(gp);
    if is_error(result) then
        Print("Test_vertices_exact, name=", name, " error\n");
        return false;
    fi;
    if Set(result)<>Set(LVert) then
        Print("Test_vertices_exact, name=", name, " result=", result, "\n");
        return false;
    fi;
    return true;
end;

ListExact:=[
  rec(name:="cubes", test:=function()
        return Test_difference_exact("cubes", [Box([0,0,0],[2,2,2])], [Box([1,1,1],[3,3,3])], 1, 7, 7);
      end),
  rec(name:="L_square", test:=function()
        return Test_difference_exact("L_square", [Box([0,0],[2,1]), Box([0,1],[1,2])], [Box([0,0],[2,2])], 3, 0, 1);
      end),
  rec(name:="segments", test:=function()
        return Test_difference_exact("segments", [Box([0],[1]), Box([1],[2])], [Box([1/2],[3])], 3/2, 1/2, 1);
      end),
  rec(name:="L_vertices", test:=function()
        return Test_vertices_exact("L_vertices", [Box([0,0],[2,1]), Box([0,1],[1,2])], [[1,0,0],[1,2,0],[1,2,1],[1,1,1],[1,1,2],[1,0,2]]);
      end),
  rec(name:="segment_vertices", test:=function()
        return Test_vertices_exact("segment_vertices", [Box([0],[1]), Box([1],[2])], [[1,0],[1,2]]);
      end),
  rec(name:="rectangle_vertices", test:=function()
        return Test_vertices_exact("rectangle_vertices", [Box([0,0],[1,1]), Box([1,0],[2,1])], [[1,0,0],[1,2,0],[1,2,1],[1,0,1]]);
      end)];
n_error_exact:=0;
for eExact in ListExact
do
    Print("Exact case name=", eExact.name, "\n");
    if eExact.test()=false then
        n_error_exact:=n_error_exact + 1;
    fi;
od;
Print("n_error_exact=", n_error_exact, "\n");

if RecResult.n_error > 0 or n_error_exact > 0 then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

