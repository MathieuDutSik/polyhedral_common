Read("../common.g");
Read("../access_points.g");
Print("Beginning TestDelaunayEnumeration\n");


# The methods that claim to return a characteristic vector set in the sense
# of Definition 1.2.1 of "A canonical form for positive definite matrices",
# that is a family generating Z^n. The other methods return families that
# are only required to be of full rank.
ListMethodSpanning:=["spanning", "cv"];

TestGeneration:=function(matrix, method)
    local U, n, divs;
    U:=get_fullrank_invariant_family(matrix, method);
    if is_error(U) then
        return false;
    fi;
    Print("|U|=", Length(U), "\n");
    n:=Length(matrix);
    if RankMat(U) <> n then
        Print("    The family is not of full rank\n");
        return false;
    fi;
    if Position(ListMethodSpanning, method) <> fail then
        divs:=ElementaryDivisorsMat(U);
        if Length(divs) <> n or Filtered(divs, x->x<>1) <> [] then
            Print("    The family does not generate Z^n, divisors=", divs, "\n");
            return false;
        fi;
    fi;
    return true;
end;

ListRec:=ReadAsFunction("ListCases")();;
ListMethod:=["shortest", "relevant_voronoi", "filtered_relevant_voronoi", "fullrank", "spanning", "cv", "cv_fullrank"];


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

