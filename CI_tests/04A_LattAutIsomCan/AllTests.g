Read("../common.g");
Read("../access_points.g");
Print("Beginning TestCanonicalization\n");

ListMat:=[];

ClassicMatrices:=true;
ConwayExample:=true;
GlueLattices:=true;
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

# The lattice s Z^n + sum_g Z g for a list of glue vectors g, expressed in
# a basis. The glue vectors are chosen so that the minimal vectors are
# exactly the +-s e_i: they span the sublattice s Z^n whose index is the
# order of the glue group. Those lattices exercise the canonicalization of
# vector families of full rank that do not span the whole lattice.
get_glue_lattice_gram:=function(n, s, glues)
    local B, pivots, g, i, eLine;
    B:=[];
    pivots:=[];
    for g in glues
    do
        Add(B, ShallowCopy(g));
        for i in [1..n]
        do
            if g[i] mod s <> 0 then
                Add(pivots, i);
                break;
            fi;
        od;
    od;
    for i in [1..n]
    do
        if not i in pivots then
            eLine:=ListWithIdenticalEntries(n, 0);
            eLine[i]:=s;
            Add(B, eLine);
        fi;
    od;
    return B * TransposedMat(B);
end;

if GlueLattices then
    # index 2: 2 Z^10 glued by (1,1,1,1,1,0,0,0,0,0)
    Add(ListMat, get_glue_lattice_gram(10, 2, [[1,1,1,1,1,0,0,0,0,0]]));
    # index 3: 3 Z^10 glued by the all-one vector
    Add(ListMat, get_glue_lattice_gram(10, 3, [[1,1,1,1,1,1,1,1,1,1]]));
    # index 4: 2 Z^10 glued by two disjoint half-support vectors
    Add(ListMat, get_glue_lattice_gram(10, 2, [[1,1,1,1,1,0,0,0,0,0],[0,0,0,0,0,1,1,1,1,1]]));
fi;

if WellRoundedDim10 then
    FileSave:="../21B_ShortRealizability/ListGram_n10_rnk10";
    if IsExistingFile(FileSave) then
        ListGram:=ReadAsFunction(FileSave)();
        Append(ListMat, ListGram);
    fi;
fi;


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
    return eP * eMat * TransposedMat(eP);
end;


test_canonicalization_function:=function(eMat)
    local TheCan, iter, eMat_B, TheCan_B;
    TheCan:=get_latt_canonical_form(eMat);
    if is_error(TheCan) then
        return fail;
    fi;
    for iter in [1..10]
    do
        Print("  |eMat|=", Length(eMat), " iter=", iter, "\n");
        eMat_B:=get_random_conjugate(eMat);
        TheCan_B:=get_latt_canonical_form(eMat_B);
        if is_error(TheCan_B) then
            return fail;
        fi;
        if TheCan_B.eG<>TheCan.eG then
            Print("test_canonicalization_function, the canonical forms do not match\n");
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
        Print("         iMat=", iMat, "/", nMat, " |eMat|=", Length(eMat), " n_error_can=", n_error_can, "\n");
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


test_automorphism_function:=function(eMat)
    local GRPperm, ord_grp, iter, eMat_B, GRPperm_B, ord_grp_b;
    GRPperm:=get_latt_automorphism_perm_group(eMat);
    if GRPperm=fail then
        Print("test_automorphism_function, inconsistent GRPperm\n");
        return fail;
    fi;
    ord_grp:=Order(GRPperm);
    for iter in [1..10]
    do
        Print("  |eMat|=", Length(eMat), " iter=", iter, "\n");
        eMat_B:=get_random_conjugate(eMat);
        GRPperm_B:=get_latt_automorphism_perm_group(eMat_B);
        if GRPperm_B=fail then
            Print("test_automorphism_function, inconsistent GRPperm_B\n");
            return fail;
        fi;
        ord_grp_b:=Order(GRPperm_B);
        if ord_grp_b<>ord_grp then
            Print("test_automorphism_function, inconsistent ord_grp\n");
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
        test:=test_automorphism_function(eMat);
        if test=fail then
            n_error_aut:=n_error_aut+1;
        fi;
    od;
    Print("n_error_aut=", n_error_aut, "\n");
    return n_error_aut;
end;


test_isomorphism_function:=function(eMat)
    local iter, eMat_B, test_iso;
    for iter in [1..10]
    do
        Print("  |eMat|=", Length(eMat), " iter=", iter, "\n");
        eMat_B:=get_random_conjugate(eMat);
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


# Lattices whose automorphism group order is known, so that the tests above
# are not only comparing the code with itself. Checking that the order is
# invariant under conjugation catches an order that varies; it does not catch
# one that is systematically wrong, and a family of full rank that does not
# span the lattice did once give an order too large by a factor two.
#
# 2 (n+1)! for A_n, 2^n n! for Z^n and for D_n with n >= 5, 1152 for D_4
# where triality adds a factor 3, and the Weyl group orders for E_6, E_7 and
# E_8. The orthogonal lattice diag(1, 4, 9) has only the sign changes.
get_An:=function(n)
    local M, i;
    M:=NullMat(n, n);
    for i in [1..n]
    do
        M[i][i]:=2;
        if i < n then
            M[i][i+1]:=-1;
            M[i+1][i]:=-1;
        fi;
    od;
    return M;
end;

get_Dn:=function(n)
    local M, i;
    M:=NullMat(n, n);
    for i in [1..n]
    do
        M[i][i]:=2;
    od;
    for i in [1..n-2]
    do
        M[i][i+1]:=-1;
        M[i+1][i]:=-1;
    od;
    M[n-2][n]:=-1;
    M[n][n-2]:=-1;
    return M;
end;

get_En:=function(n)
    local M, i;
    M:=NullMat(n, n);
    for i in [1..n]
    do
        M[i][i]:=2;
    od;
    for i in [1..n-2]
    do
        M[i][i+1]:=-1;
        M[i+1][i]:=-1;
    od;
    M[3][n]:=-1;
    M[n][3]:=-1;
    return M;
end;

ListKnownOrder:=[
  rec(name:="A2", eG:=get_An(2), ord:=12),
  rec(name:="A3", eG:=get_An(3), ord:=48),
  rec(name:="A4", eG:=get_An(4), ord:=240),
  rec(name:="A5", eG:=get_An(5), ord:=1440),
  rec(name:="D4", eG:=get_Dn(4), ord:=1152),
  rec(name:="D5", eG:=get_Dn(5), ord:=3840),
  rec(name:="E6", eG:=get_En(6), ord:=103680),
  rec(name:="E7", eG:=get_En(7), ord:=2903040),
  rec(name:="E8", eG:=get_En(8), ord:=696729600),
  rec(name:="Z2", eG:=IdentityMat(2), ord:=8),
  rec(name:="Z3", eG:=IdentityMat(3), ord:=48),
  rec(name:="Z4", eG:=IdentityMat(4), ord:=384),
  rec(name:="diag149", eG:=DiagonalMat([1,4,9]), ord:=8)];

test_known_orders:=function()
    local n_error, eRec, GRPperm, ord;
    n_error:=0;
    for eRec in ListKnownOrder
    do
        GRPperm:=get_latt_automorphism_perm_group(eRec.eG);
        if GRPperm=fail then
            Print("test_known_orders, ", eRec.name, ": no group\n");
            n_error:=n_error+1;
        else
            ord:=Order(GRPperm);
            Print("  ", eRec.name, " |Aut|=", ord, " expected=", eRec.ord, "\n");
            if ord <> eRec.ord then
                Print("test_known_orders, ", eRec.name, ": wrong order\n");
                n_error:=n_error+1;
            fi;
        fi;
    od;
    Print("n_error_known_order=", n_error, "\n");
    return n_error;
end;

# Lattices of the same dimension that are not isometric. Testing only that a
# lattice is isomorphic to its own conjugates cannot catch a test that
# answers yes too often.
ListNonIsomorphic:=[
  rec(name1:="A2", eG1:=get_An(2), name2:="Z2", eG2:=IdentityMat(2)),
  rec(name1:="A3", eG1:=get_An(3), name2:="Z3", eG2:=IdentityMat(3)),
  rec(name1:="Z3", eG1:=IdentityMat(3), name2:="diag149", eG2:=DiagonalMat([1,4,9])),
  rec(name1:="D4", eG1:=get_Dn(4), name2:="Z4", eG2:=IdentityMat(4)),
  rec(name1:="D4", eG1:=get_Dn(4), name2:="A4", eG2:=get_An(4)),
  rec(name1:="D5", eG1:=get_Dn(5), name2:="A5", eG2:=get_An(5))];

test_non_isomorphic:=function()
    local n_error, eRec, test_iso;
    n_error:=0;
    for eRec in ListNonIsomorphic
    do
        test_iso:=get_latt_isomorphism_test(eRec.eG1, eRec.eG2);
        Print("  ", eRec.name1, " against ", eRec.name2, " iso=", test_iso, "\n");
        if test_iso=fail then
            Print("test_non_isomorphic: program failure\n");
            n_error:=n_error+1;
        else
            if test_iso<>false then
                Print("test_non_isomorphic: not isometric but found isomorphic\n");
                n_error:=n_error+1;
            fi;
        fi;
    od;
    Print("n_error_non_isomorphic=", n_error, "\n");
    return n_error;
end;


test_all:=function()
    local n_error;
    n_error:=0;
    n_error:=n_error + test_all_cans();
    n_error:=n_error + test_all_automs();
    n_error:=n_error + test_all_isoms();
    n_error:=n_error + test_known_orders();
    n_error:=n_error + test_non_isomorphic();
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
