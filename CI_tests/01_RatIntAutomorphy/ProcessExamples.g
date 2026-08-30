Read("../common.g");
Read("../access_points.g");

# Unified test of the rational automorphy codes, of the integral
# automorphy codes and of the double coset codes. It replaces the
# former separate CI sections 01_MatrixDoubleCosets,
# 02A_IntegralAutomorphy and 02B_RationalAutomorphy.
#
# The examples are the classification of simplices in dimension 5 to 7
# and the simple root systems of the reflective forms in dimension 4
# and 5. As in the former sections, the root systems are used for the
# integral tests only: the rational section had AppendReflectiveDim45
# set to false and the double coset tests enumerate the groups element
# by element, which is out of reach for those examples.

ListEXT_simplices:=[];
for n in [5..7]
do
    FileName:=Concatenation("../DATA/ClassificationSimplices", String(n));
    Append(ListEXT_simplices, ReadAsFunction(FileName)());
od;

ListEXT_reflective:=ReadAsFunction("ListSimpleRootSystem_4_56_X_5_47")();

ListEXT_integral:=Concatenation(ListEXT_simplices, ListEXT_reflective);


#
# Shared helpers
#

scramble_ext:=function(EXT)
    local n, n_vert, P;
    n:=Length(EXT[1]);
    n_vert:=Length(EXT);
    P:=RandomIntegralUnimodularMatrix(n);
    return Permuted(EXT * P, Random(SymmetricGroup(n_vert)));
end;


#
# The rational tests
#

TestCase_Rational_Automorphy:=function(EXT)
    local EXT_img, rec_gap1, rec_gap2;
    EXT_img:=scramble_ext(EXT);
    rec_gap1:=get_grp_automorphy(EXT);
    if is_error(rec_gap1) then
        return false;
    fi;
    rec_gap2:=get_grp_automorphy(EXT_img);
    if is_error(rec_gap2) then
        return false;
    fi;
    if Order(rec_gap1.GAPperm)<>Order(rec_gap2.GAPperm) then
        return false;
    fi;
    return true;
end;

TestCase_Rational_Canonic:=function(EXT)
    local EXT_img, EXT_can1, EXT_can2;
    EXT_img:=scramble_ext(EXT);
    EXT_can1:=get_canonic_form(EXT);
    if is_error(EXT_can1) then
        return false;
    fi;
    EXT_can2:=get_canonic_form(EXT_img);
    if is_error(EXT_can2) then
        return false;
    fi;
    if EXT_can1<>EXT_can2 then
        Print("Different canonical form\n");
        return false;
    fi;
    return true;
end;

TestCase_Rational_Isomorphy:=function(EXT)
    local EXT_img, result;
    EXT_img:=scramble_ext(EXT);
    result:=get_isomorphism_result(EXT, EXT_img);
    if is_error(result) then
        return false;
    fi;
    if result=fail then
        # Found to be non-isomorphic
        return false;
    fi;
    return true;
end;


#
# The integral tests
#

GeneratorsPreservePolytope:=function(TheGRP, EXT)
    local gens, BasisEXT, BasisIndices, M, n, i, ImageEXT, ImageIndices, Mb, Mi, A, g;
    BasisEXT:=[];
    BasisIndices:=[];
    n:=Length(EXT[1]);
    for i in [1..Length(EXT)] do
        M:=Concatenation(BasisEXT,[EXT[i]]);
        if RankMat(M)=Length(BasisEXT)+1 then
            Add(BasisEXT, EXT[i]);
            Add(BasisIndices, i);
        fi;
        if Length(BasisEXT)=n then
            break;
        fi;
    od;
    gens:=GeneratorsOfGroup(TheGRP);
    for g in gens do
        ImageIndices:=List(BasisIndices, i -> i^g);
        ImageEXT := List(ImageIndices, i -> EXT[i]);
        Mb:=TransposedMat(BasisEXT);
        if DeterminantMat(Mb) = 0 then
            return false;
        fi;
        Mi:=TransposedMat(ImageEXT);
        A:=Mi*Inverse(Mb);
        if not ForAll(Flat(A), IsInt) then
            return false;
        fi;
        for i in Difference([1..Length(EXT)], BasisIndices) do
            if not A*EXT[i] = EXT[i^g] then
                return false;
            fi;
        od;
    od;
    return true;
end;

TestCase_Integral_Automorphy:=function(EXT)
    local TheGRP;
    TheGRP:=get_grp_integral_automorphy(EXT);
    if is_error(TheGRP) then
        return false;
    fi;
    if not GeneratorsPreservePolytope(TheGRP.GAPperm, EXT) then
        return false;
    fi;
    Print("|TheGRP|=", Order(TheGRP.GAPperm), "\n");
    return true;
end;

TestCase_Integral_Isomorphy:=function(EXT)
    local EXT_img, result;
    EXT_img:=scramble_ext(EXT);
    result:=test_polytope_integral_equivalence(EXT, EXT_img);
    if is_error(result) then
        return false;
    fi;
    if result=fail then
        return false;
    fi;
    return true;
end;

TestCase_Integral_Automorphy_RightCoset:=function(EXT)
    local result;
    result:=get_linpolytopeintegral_aut_rightcoset(EXT);
    if is_error(result) then
        return false;
    fi;
    return true;
end;


#
# The double coset tests.
#
# We check that G_rat = \cup_i G_int g_i G_int is a partition of the
# rational group into double cosets of the integral group.
#

TestDoubleCosetPermDecomposition:=function(GRP, GRP_U, GRP_V, DoubleCosetsPerm)
    local ListSets, sumElt, eCos, eU, eV, eG, iSet, jSet, eInt, LGen, eGen, NewGen, GRP_v_conj, GRPinter, ListRightCos, eRightCos, fRightCos, Lentry, eEntry, eList_method1, eList_method2, eSet_method1, eSet_method2, method1;
    Print("|G|=", Order(GRP), " |GRP_U|=", Order(GRP_U), " |GRP_V|=", Order(GRP_V), " |LCos|=", Length(DoubleCosetsPerm), "\n");
    ListSets:=[];
    sumElt:=0;
    for eCos in DoubleCosetsPerm
    do
        LGen:=[];
        for eGen in GeneratorsOfGroup(GRP_V)
        do
            NewGen:=eCos * eGen * Inverse(eCos);
            Add(LGen, NewGen);
        od;
        GRP_v_conj:=Group(LGen);
        GRPinter:=Intersection(GRP_v_conj, GRP_U);
        ListRightCos:=RightCosets(GRP_v_conj, GRPinter);
        Lentry:=[];
        for eRightCos in ListRightCos
        do
            fRightCos:=Representative(eRightCos) * eCos;
            Add(Lentry, fRightCos);
        od;
        #
        eList_method2:=[];
        for eU in GRP_U
        do
            for eEntry in Lentry
            do
                eG:=eU * eEntry;
                Add(eList_method2, eG);
            od;
        od;
        eSet_method2:=Set(eList_method2);
        #
        method1:=false;
        if method1 then
            eList_method1:=[];
            for eU in GRP_U
            do
                for eV in GRP_V
                do
                    eG:=eU * eCos * eV;
                    Add(eList_method1, eG);
                od;
            od;
            eSet_method1:=Set(eList_method1);
            if eSet_method1<>eSet_method2 then
                Error("The computation methods are inconsistent");
            fi;
        fi;

        sumElt:=sumElt + Length(eSet_method2);
        Add(ListSets, eSet_method2);
    od;
    Print("Creation of the double cosets as raw sets done\n");
    for iSet in [1..Length(ListSets)]
    do
        for jSet in [iSet+1..Length(ListSets)]
        do
            eInt:=Intersection(ListSets[iSet], ListSets[jSet]);
            if Length(eInt) > 0 then
                Print("Intersection are not empty\n");
                return false;
            fi;
        od;
    od;
    Print("Test intersection done\n");
    if sumElt<>Order(GRP) then
        Print("The union of all cosets is not the full group\n");
        return false;
    fi;
    Print("All checks done\n");
    return true;
end;

# Common preliminary of the two double coset tests: the rational group
# and the integral group, the latter having to be a subgroup of the
# former.
GetRationalAndIntegralGroups:=function(EXT)
    local GRP_rat, GRP_V;
    GRP_rat:=get_grp_automorphy(EXT);
    if is_error(GRP_rat) then
        Print("GRP_rat=fail which we do not want\n");
        return fail;
    fi;
    GRP_V:=get_grp_integral_automorphy(EXT);
    if is_error(GRP_V) then
        Print("GRP_V=fail which we do not want\n");
        return fail;
    fi;
    if IsSubgroup(GRP_rat.GAPperm, GRP_V.GAPperm)=false then
        Print("|EXT|=", Length(EXT), " / ", Length(EXT[1]), " Det=", DeterminantMat(EXT), "\n");
        Print("EXT=\n");
        PrintArray(EXT);
        Print("GRP_rat=", GeneratorsOfGroup(GRP_rat.GAPperm), " |GRP_rat|=", Order(GRP_rat.GAPperm), "\n");
        Print("GRP_V=", GeneratorsOfGroup(GRP_V.GAPperm), " |GRP_V|=", Order(GRP_V.GAPperm), "\n");
        Print("The integral group should be a subgroup of the rational group\n");
        return fail;
    fi;
    return rec(GRP_rat:=GRP_rat, GRP_V:=GRP_V);
end;

TestCase_LinPolytopeIntegral_Automorphy_DoubleCoset:=function(EXT)
    local eRecGrp, GRP_rat, GRP_V, GRP_U, RecResult;
    Print("Begin TestCase_LinPolytopeIntegral_Automorphy_DoubleCoset, Det(BaseIntMat(EXT))=", DeterminantMat(BaseIntMat(EXT)), "\n");
    eRecGrp:=GetRationalAndIntegralGroups(EXT);
    if eRecGrp=fail then
        return false;
    fi;
    GRP_rat:=eRecGrp.GRP_rat;
    GRP_V:=eRecGrp.GRP_V;
    RecResult:=get_double_cosets_matrix_group(EXT, GRP_V.GAPperm);
    if is_error(RecResult) then
        return false;
    fi;
    GRP_U:=RecResult.GAPperm;
    if GRP_U<>GRP_V.GAPperm then
        Print("GRP_U and GRP_V should be equal\n");
        return false;
    fi;
    return TestDoubleCosetPermDecomposition(GRP_rat.GAPperm, GRP_U, GRP_V.GAPperm, RecResult.DoubleCosetsPerm);
end;

TestCase_LinearSpace_Stabilizer_DoubleCoset:=function(EXT)
    local dim, eRecGrp, GRP_rat, TheSpace, GRP_V, GRP_U, RecResult, DoubleCosetsPerm;
    dim:=Length(EXT[1]);
    Print("Begin TestCase_LinearSpace_Stabilizer_DoubleCoset, Det(BaseIntMat(EXT))=", DeterminantMat(BaseIntMat(EXT)), "\n");
    eRecGrp:=GetRationalAndIntegralGroups(EXT);
    if eRecGrp=fail then
        return false;
    fi;
    GRP_rat:=eRecGrp.GRP_rat;
    GRP_V:=eRecGrp.GRP_V;
    TheSpace:=IdentityMat(dim);
    RecResult:=get_linear_space_stabilizer_double_cosets(GRP_rat.GAPmatr, TheSpace, GRP_V.GAPmatr);
    if is_error(RecResult) then
        return false;
    fi;
    Print("We have RecResult\n");
    DoubleCosetsPerm:=List(RecResult.ListCos, xMat->PermList(List(EXT, xVert->Position(EXT, xVert*xMat))));
    GRP_U:=Group(List(GeneratorsOfGroup(RecResult.GRPmatr), xMat->PermList(List(EXT, xVert->Position(EXT, xVert*xMat)))));
    return TestDoubleCosetPermDecomposition(GRP_rat.GAPperm, GRP_U, GRP_V.GAPperm, DoubleCosetsPerm);
end;


#
# Running the test suites
#

RunTestSuite:=function(SuiteName, ListTests, ListEXT)
    local n_err, nEXT, iEXT, EXT, eTest, test;
    n_err:=0;
    nEXT:=Length(ListEXT);
    for iEXT in [1..nEXT]
    do
        EXT:=ListEXT[iEXT];
        Print("---------------------------------------- ", SuiteName, " ", iEXT, "/", nEXT, " |EXT|=", Length(EXT), "/", Length(EXT[1]), " n_err=", n_err, " ----------------------------------------\n");
        for eTest in ListTests
        do
            Print("---- ", eTest[1], "\n");
            test:=eTest[2](EXT);
            if test=false then
                Print("ERROR in ", eTest[1], "\n");
                n_err:=n_err + 1;
            fi;
        od;
    od;
    Print(SuiteName, " suite done, n_err=", n_err, "\n");
    return n_err;
end;

CI_Decision_Reset();
n_error:=0;

n_error:=n_error + RunTestSuite("Rational",
    [["Rational_Automorphy", TestCase_Rational_Automorphy],
     ["Rational_Canonic", TestCase_Rational_Canonic],
     ["Rational_Isomorphy", TestCase_Rational_Isomorphy]],
    ListEXT_simplices);

n_error:=n_error + RunTestSuite("Integral",
    [["Integral_Automorphy", TestCase_Integral_Automorphy],
     ["Integral_Isomorphy", TestCase_Integral_Isomorphy],
     ["Integral_Automorphy_RightCoset", TestCase_Integral_Automorphy_RightCoset]],
    ListEXT_integral);

n_error:=n_error + RunTestSuite("DoubleCoset",
    [["LinPolytopeIntegral_Automorphy_DoubleCoset", TestCase_LinPolytopeIntegral_Automorphy_DoubleCoset],
     ["LinearSpace_Stabilizer_DoubleCoset", TestCase_LinearSpace_Stabilizer_DoubleCoset]],
    ListEXT_simplices);

Print("-------------------------------------------------------\n");
Print("n_error=", n_error, "\n");
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
