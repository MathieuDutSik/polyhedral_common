Read("../common.g");
Read("../access_points.g");
Print("Beginning TestDelaunayEnumeration\n");


# What each method promises, which is not the same for all of them.
#
# Every one of them is a canonical construction, so all are required to be
# equivariant: V(U A U^T) = V(A) U^(-1).
#
# Only some are required to be of full rank. "shortest" returns the shortest
# vectors and nothing more, and those span a proper subspace as soon as the
# lattice is not well rounded: on DualLambda9 there are two of them.
# "relevant_voronoi" returns the Voronoi relevant vectors, which do span, but
# the method is not built around that promise.
ListMethodFullRank:=["fullrank", "spanning", "cv", "cv_fullrank"];
#
# And only those built to generate Z^n, in the sense of Definition 1.2.1 of
# "A canonical form for positive definite matrices", are required to do so.
ListMethodSpanning:=["spanning", "cv"];

# A random conjugate U eMat U^T with U in GL_n(Z).
get_random_conjugate:=function(eMat)
    local n, GRP, LGen, len, eP, i;
    n:=Length(eMat);
    GRP:=GeneralLinearGroup(n, Integers);
    LGen:=GeneratorsOfGroup(GRP);
    len:=Random([1..2*n]);
    eP:=IdentityMat(n);
    for i in [1..len]
    do
        eP := eP * Random(LGen);
    od;
    return rec(eP:=eP, eMat:=eP * eMat * TransposedMat(eP));
end;

#
#  The second half of Definition 1.2.1 of "A canonical form for positive
#   definite matrices": a characteristic vector set function satisfies
#   V(U A U^T) = V(A) U^(-1). Without it the family is invariant only in the
#   sense that it has the same size, which is not what the canonical form
#   needs, and the failure is silent: E8 was enough to make a family that had
#   passed every other test give a canonical form depending on the labelling
#   of the input.
TestEquivariance:=function(matrix, method)
    local rec_conj, V, V_conj, iter;
    V:=get_fullrank_invariant_family(matrix, method);
    if is_error(V) then
        return false;
    fi;
    for iter in [1..3]
    do
        rec_conj:=get_random_conjugate(matrix);
        V_conj:=get_fullrank_invariant_family(rec_conj.eMat, method);
        if is_error(V_conj) then
            return false;
        fi;
        if Set(V_conj) <> Set(V * Inverse(rec_conj.eP)) then
            Print("    The family is not equivariant, method=", method, "\n");
            return false;
        fi;
    od;
    return true;
end;

TestGeneration:=function(matrix, method)
    local U, n, divs;
    U:=get_fullrank_invariant_family(matrix, method);
    if is_error(U) then
        return false;
    fi;
    Print("|U|=", Length(U), "\n");
    n:=Length(matrix);
    if Position(ListMethodFullRank, method) <> fail then
        if RankMat(U) <> n then
            Print("    The family is not of full rank\n");
            return false;
        fi;
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
            result:=TestEquivariance(eRec.eG, method);
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

