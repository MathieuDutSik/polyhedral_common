Read("../common.g");

ListEXT:=[];
for n in [5..7]
do
    FileName:=Concatenation("../DATA/ClassificationSimplices", String(n));
    PartListEXT:=ReadAsFunction(FileName)();
    Append(ListEXT, PartListEXT);
od;

get_grp_automorphy:=function(EXT)
    local FileI, FileO, arith, eProg, TheCommand, TheGRP_rat;
    FileI:=Filename(DirectoryTemporary(), "Test.in");
    FileO:=Filename(DirectoryTemporary(), "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:="../../src_group/GRP_LinPolytope_Automorphism";
    TheCommand:=Concatenation(eProg, " rational ", FileI, " RecGAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return fail;
    fi;
    TheGRP_rat:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    Print("|TheGRP_rat|=", Order(TheGRP_rat.GAPperm), "\n");
    return TheGRP_rat;
end;

get_grp_integral_automorphy:=function(EXT)
    local FileI, FileO, arith, eProg, TheCommand, TheGRP_int;
    FileI:=Filename(DirectoryTemporary(), "Test.in");
    FileO:=Filename(DirectoryTemporary(), "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:="../../src_group/GRP_LinPolytopeIntegral_Automorphism";
    TheCommand:=Concatenation(eProg, " rational ", FileI, " RecGAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return fail;
    fi;
    TheGRP_int:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    Print("|TheGRP_int|=", Order(TheGRP_int.GAPperm), "\n");
    return TheGRP_int;
end;

# Compute the double cosets of G_rat = \cup_i G_int g_i G_int.


TestDoubleCosetPermDecomposition:=function(GRP, GRP_U, GRP_V, DoubleCosetsPerm)
    local ListSets, sumElt, eCos, eU, eV, eG, eList, eSet, iSet, jSet, eInt, LGen, eGen, NewGen, GRP_v_conj, GRPinter, ListRightCos, eRightCos, fRightCos, Lentry, eEntry, eList_method1, eList_method2, eSet_method1, eSet_method2, method1;
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




TestCase_LinPolytopeIntegral_Automorphy_DoubleCoset:=function(EXT)
    local FileI, FileO, FileGRP_V, arith, OutFormat, eProg, TheCommand, GRP_rat, GRP_V, GRP_U, RecResult;
    FileI:=Filename(DirectoryTemporary(), "Test.in");
    FileO:=Filename(DirectoryTemporary(), "Test.out");
    FileGRP_V:=Filename(DirectoryTemporary(), "Test.grp_V");
    WriteMatrixFile(FileI, EXT);
    Print("Begin TestCase_Automorphy_DoubleCoset, Det(BaseIntMat(EXT))=", DeterminantMat(BaseIntMat(EXT)), "\n");
    GRP_rat:=get_grp_automorphy(EXT);
    if GRP_rat=fail then
        Print("GRP_rat=fail which we do not want\n");
        return false;
    fi;
    GRP_V:=get_grp_integral_automorphy(EXT);
    if GRP_V=fail then
        Print("GRP_V=fail which we do not want\n");
        return false;
    fi;
    if IsSubgroup(GRP_rat.GAPperm, GRP_V.GAPperm)=false then
        Print("|EXT|=", Length(EXT), " / ", Length(EXT[1]), " Det=", DeterminantMat(EXT), "\n");
        Print("EXT=\n");
        PrintArray(EXT);
        Print("GRP_rat=", GeneratorsOfGroup(GRP_rat.GAPperm), " |GRP_rat|=", Order(GRP_rat.GAPperm), "\n");
        Print("GRP_V=", GeneratorsOfGroup(GRP_V.GAPperm), " |GRP_V|=", Order(GRP_V.GAPperm), "\n");
        Print("The integral group should be a subgroup of the rational group\n");
        return false;
    fi;
    WriteGroupFile(FileGRP_V, Length(EXT), GRP_V.GAPperm);
    eProg:="../../src_group/GRP_LinPolytopeIntegral_Automorphism_DoubleCoset";
    TheCommand:=Concatenation(eProg, " rational ", FileI, " ", FileGRP_V, " RecGAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("The FileO does not exist\n");
        return false;
    fi;
    RecResult:=ReadAsFunction(FileO)();
    Print("We have RecResult\n");
    GRP_U:=RecResult.GAPperm;
    if GRP_U<>GRP_V.GAPperm then
        Print("GRP_U and GRP_V should be equal\n");
        return false;
    fi;
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileGRP_V);
    return TestDoubleCosetPermDecomposition(GRP_rat.GAPperm, GRP_U, GRP_V.GAPperm, RecResult.DoubleCosetsPerm);
end;


TestCase_LinearSpace_Stabilizer_DoubleCoset:=function(EXT)
    local dim, FileGRP, FileSPA, FileGRP_V, FileO, arith, OutFormat, eProg, TheCommand, GRP_rat, TheSpace, GRP_V, GRP_U, RecResult, DoubleCosetsPerm;
    dim:=Length(EXT[1]);
    FileGRP:=Filename(DirectoryTemporary(), "Test.grp");
    FileSPA:=Filename(DirectoryTemporary(), "Test.space");
    FileGRP_V:=Filename(DirectoryTemporary(), "Test.grp_V");
    FileO:=Filename(DirectoryTemporary(), "Test.out");
    Print("Begin TestCase_Automorphy_DoubleCoset, Det(BaseIntMat(EXT))=", DeterminantMat(BaseIntMat(EXT)), "\n");
    GRP_rat:=get_grp_automorphy(EXT);
    if GRP_rat=fail then
        Print("GRP_rat=fail which we do not want\n");
        return false;
    fi;
    TheSpace:=IdentityMat(dim);
    GRP_V:=get_grp_integral_automorphy(EXT);
    if GRP_V=fail then
        Print("GRP_V=fail which we do not want\n");
        return false;
    fi;
    if IsSubgroup(GRP_rat.GAPperm, GRP_V.GAPperm)=false then
        Print("|EXT|=", Length(EXT), " / ", Length(EXT[1]), " Det=", DeterminantMat(EXT), "\n");
        Print("EXT=\n");
        PrintArray(EXT);
        Print("GRP_rat=", GeneratorsOfGroup(GRP_rat.GAPperm), " |GRP_rat|=", Order(GRP_rat.GAPperm), "\n");
        Print("GRP_V=", GeneratorsOfGroup(GRP_V.GAPperm), " |GRP_V|=", Order(GRP_V.GAPperm), "\n");
        Print("The integral group should be a subgroup of the rational group\n");
        return false;
    fi;
    WriteListMatrixFile(FileGRP, GRP_rat.GAPmatr);
    WriteMatrixFile(FileSPA, TheSpace);
    WriteListMatrixFile(FileGRP_V, GRP_V.GAPmatr);
    eProg:="../../src_group/GRP_LinearSpace_Stabilizer_DoubleCoset";
    TheCommand:=Concatenation(eProg, " ", FileGRP, " ", FileSPA, " ", FileGRP_V, " ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("The FileO does not exist\n");
        return false;
    fi;
    RecResult:=ReadAsFunction(FileO)();
    RemoveFile(FileGRP);
    RemoveFile(FileSPA);
    RemoveFile(FileGRP_V);
    RemoveFile(FileO);
    Print("We have RecResult\n");
    DoubleCosetsPerm:=List(RecResult.ListCos, xMat->PermList(List(EXT, xVert->Position(EXT, xVert*xMat))));
    GRP_U:=Group(List(GeneratorsOfGroup(RecResult.GRPmatr), xMat->PermList(List(EXT, xVert->Position(EXT, xVert*xMat)))));
    return TestDoubleCosetPermDecomposition(GRP_rat.GAPperm, GRP_U, GRP_V.GAPperm, DoubleCosetsPerm);
end;




#

CI_Decision_Reset();
pos:=0;
n_error:=0;

for EXT in ListEXT
do
    Print("     pos=", pos, "/", Length(ListEXT), " |EXT|=", Length(EXT), "/", Length(EXT[1]), "\n");
    test:=TestCase_LinPolytopeIntegral_Automorphy_DoubleCoset(EXT);
    if test=false then
        Error("Stop here");
        n_error:=n_error + 1;
    fi;
    test:=TestCase_LinearSpace_Stabilizer_DoubleCoset(EXT);
    if test=false then
        Error("Stop here");
        n_error:=n_error + 1;
    fi;
    pos:=pos + 1;
od;

Print("n_error=", n_error, "\n");
if n_error=0 then
    Print("Normal case\n");
    CI_Write_Ok();
else
    Print("Error case\n");
fi;
