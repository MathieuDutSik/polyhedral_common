Read("../common.g");
Read("../access_points.g");
Print("Beginning TestDelaunayEnumeration\n");


TestGeneration:=function(matrix, method)
    local U;
    U:=get_fullrank_invariant_family(matrix, method);
    if is_error(U) then
        return false;
    fi;
    Print("|U|=", Length(U), "\n");
    return true;
end;

ListRec:=ReadAsFunction("ListCases")();;
ListMethod:=["shortest", "relevant_voronoi", "filtered_relevant_voronoi", "fullrank", "spanning"];


FullTest:=function()
    local iRec, eRec, result, i_meth, n_meth, method;
    iRec:=0;
    for eRec in ListRec
    do
        Print("----------------------------------------------------------------------------\n");
        Print("iRec=", iRec, "/", Length(ListRec), " Treating lattice named ", eRec.name, "\n");
        n_meth:=Length(ListMethod);
        for i_meth in [1..n_meth]
        do
            method:=ListMethod[i_meth];
            Print("    i_meth=", i_meth, "/", n_meth, " method=", method, "\n");
            result:=TestGeneration(eRec.eG, method);
            if result=false then
                return false;
            fi;
        od;
        iRec:=iRec+1;
    od;
    return true;
end;

result:=FullTest();
Print("2: result=", result, "\n");
CI_Decision_Reset();
if result=false then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
CI_PrintExistConclusion();

