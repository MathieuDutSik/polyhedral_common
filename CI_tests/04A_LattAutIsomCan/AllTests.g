Read("../common.g");
Read("../access_points.g");
Print("Beginning TestCanonicalization\n");

ListMat:=[];

ClassicMatrices:=true;
ConwayExample:=true;
WellRoundedDim10:=true;

keep_err:=true;
#keep_err:=false;


if ClassicMatrices then
    Add(ListMat, ClassicalSporadicLattices("A4"));
    Add(ListMat, ClassicalSporadicLattices("A5"));
    Add(ListMat, ClassicalSporadicLattices("E6"));
    Add(ListMat, ClassicalSporadicLattices("D4"));
    Add(ListMat, ClassicalSporadicLattices("D5"));
fi;

if ConwayExample then
    # Conway/Sloane lattice in dimension 11 with no basis
    Add(ListMat, ClassicalSporadicLattices("ConwaySloane11"));
fi;

if WellRoundedDim10 then
    FileSave:="../21B_ShortRealizability/ListGram_n10_rnk10";
    if IsExistingFile(FileSave) then
        ListGram:=ReadAsFunction(FileSave)();
        Append(ListMat, ListGram);
    fi;
fi;


test_canonicalization_function:=function(eMat)
    local TheCan, n, GRP, iter, len, eP, eMat_B, TheCan_B, LGen, i;
    TheCan:=get_latt_canonical_form(eMat);
    if is_error(TheCan) then
        return fail;
    fi;
    n:=Length(eMat);
    GRP:=GeneralLinearGroup(n, Integers);
    LGen:=GeneratorsOfGroup(GRP);
    for iter in [1..10]
    do
        Print("  |eMat|=", Length(eMat), " iter=", iter, "\n");
        len:=Random([1..2*n]);
        eP:=IdentityMat(n);
        for i in [1..len]
        do
            eP := eP * Random(LGen);
        od;
        eMat_B:=eP * eMat * TransposedMat(eP);
        TheCan_B:=get_latt_canonical_form(eMat_B);
        if is_error(TheCan_B) then
            return fail;
        fi;
        if TheCan_B.eG<>TheCan.eG then
            return fail;
        fi;
    od;
    return true;
end;


test_all_cans:=function()
    local n_error_can, nMat, iMat, eMat, test;
    n_error_can:=0;
    nMat:=Length(ListMat);
    for iMat in [1..nMat]
    do
        eMat:=ListMat[iMat];
        Print("         iMat=", iMat, "/", nMat, " |eMat|=", Length(eMat), " m_error_can=", n_error_can, "\n");
        test:=test_canonicalization_function(eMat);
        if test=fail then
            n_error_can:=n_error_can+1;
        fi;
    od;
    Print("n_error_can=", n_error_can, "\n");
    return n_error_can;
end;


get_latt_automorphism_perm_group:=function(eMat)
    local LmatrGens, ListVect, ListPermGens, eGenMatr, ePerm, GRPperm, SetVect;
    LmatrGens:=get_latt_automorphism_group(eMat);
    if LmatrGens=fail then
        Print("get_latt_automorphism_perm_group, LmatrGens=fail\n");
        return fail;
    fi;
    ListVect:=get_fullrank_invariant_family(eMat, "fullrank");
    if is_error(ListVect) then
        return fail;
    fi;
    SetVect:=Set(ListVect);
    ListPermGens:=[];
    for eGenMatr in LmatrGens
    do
        ePerm:=SortingPerm(SetVect * eGenMatr);
        Add(ListPermGens, ePerm);
    od;
    GRPperm:=Group(ListPermGens);
    return GRPperm;
end;


test_canonicalization_function:=function(eMat)
    local GRPperm, ord_grp, n, GRP, iter, len, eP, eMat_B, GRPperm_B, ord_grp_b, LGen, i;
    GRPperm:=get_latt_automorphism_perm_group(eMat);
    if GRPperm=fail then
        Print("test_canonicalization_function, inconsistent GRPperm\n");
        return fail;
    fi;
    ord_grp:=Order(GRPperm);
    n:=Length(eMat);
    GRP:=GeneralLinearGroup(n, Integers);
    LGen:=GeneratorsOfGroup(GRP);
    for iter in [1..10]
    do
        Print("  |eMat|=", Length(eMat), " iter=", iter, "\n");
        len:=Random([1..2*n]);
        eP:=IdentityMat(n);
        for i in [1..len]
        do
            eP := eP * Random(LGen);
        od;
        eMat_B:=eP * eMat * TransposedMat(eP);
        GRPperm_B:=get_latt_automorphism_perm_group(eMat_B);
        if GRPperm_B=fail then
            Print("test_canonicalization_function, inconsistent GRPperm_B\n");
            return fail;
        fi;
        ord_grp_b:=Order(GRPperm_B);
        if ord_grp_b<>ord_grp then
            Print("test_canonicalization_function, inconsistent ord_grp\n");
            return fail;
        fi;
    od;
    return true;
end;


test_all_automs:=function()
    local n_error_aut, nMat, iMat, eMat, test;
    n_error_aut:=0;
    nMat:=Length(ListMat);
    for iMat in [1..nMat]
    do
        eMat:=ListMat[iMat];
        Print("         iMat=", iMat, "/", nMat, " |eMat|=", Length(eMat), " n_error_aut=", n_error_aut, "\n");
        test:=test_canonicalization_function(eMat);
        if test=fail then
            n_error_aut:=n_error_aut+1;
        fi;
    od;
    Print("n_error_aut=", n_error_aut, "\n");
    return n_error_aut;
end;


test_isomorphism_function:=function(eMat)
    local TheCan, n, GRP, iter, len, eP, eMat_B, TheCan_B, LGen, i, test_iso;
    n:=Length(eMat);
    GRP:=GeneralLinearGroup(n, Integers);
    LGen:=GeneratorsOfGroup(GRP);
    for iter in [1..10]
    do
        Print("  |eMat|=", Length(eMat), " iter=", iter, "\n");
        len:=Random([1..2*n]);
        eP:=IdentityMat(n);
        for i in [1..len]
        do
            eP := eP * Random(LGen);
        od;
        eMat_B:=eP * eMat * TransposedMat(eP);
        test_iso:=get_latt_isomorphism_test(eMat, eMat_B);
        if test_iso=fail or test_iso=false then
            Print("test_isomorphism_function, incorrect result\n");
            return fail;
        fi;
    od;
    return true;
end;


test_all_isoms:=function()
    local n_error_iso, nMat, iMat, eMat, test;
    n_error_iso:=0;
    nMat:=Length(ListMat);
    for iMat in [1..nMat]
    do
        eMat:=ListMat[iMat];
        Print("         iMat=", iMat, "/", nMat, " |eMat|=", Length(eMat), " n_error_iso=", n_error_iso, "\n");
        test:=test_isomorphism_function(eMat);
        if test=fail then
            n_error_iso:=n_error_iso+1;
        fi;
    od;
    Print("n_error_iso=", n_error_iso, "\n");
    return n_error_iso;
end;


test_all:=function()
    local n_error;
    n_error:=0;
#    n_error:=n_error + test_all_cans();
    n_error:=n_error + test_all_automs();
#    n_error:=n_error + test_all_isoms();
    return n_error;
end;


n_error:=test_all();

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
