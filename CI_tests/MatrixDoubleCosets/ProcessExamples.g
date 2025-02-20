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
    TheGRP_int:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    Print("|TheGRP_int|=", Order(TheGRP_int.GAPperm), "\n");
    return TheGRP_int;
end;

# Compute the double cosets of G_rat = \cup_i G_int g_i G_int.

TestCase_Automorphy_DoubleCoset:=function(EXT)
    local FileI, FileO, FileGRP_V, arith, OutFormat, eProg, TheCommand, GRP_rat, GRP_V, GRP_U, RecResult, ListSets, sumElt, eCos, eU, eV, eG, eList, eSet, iSet, jSet, eInt;
    FileI:=Filename(DirectoryTemporary(), "Test.in");
    FileO:=Filename(DirectoryTemporary(), "Test.out");
    FileGRP_V:=Filename(DirectoryTemporary(), "Test.grp_V");
    WriteMatrixFile(FileI, EXT);
    Print("Begin TestCase_Automorphy_DoubleCoset, Det(BaseIntMat(EXT))=", DeterminantMat(BaseIntMat(EXT)), "\n");
    GRP_rat:=get_grp_automorphy(EXT);
    GRP_V:=get_grp_integral_automorphy(EXT);
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
    RecResult:=ReadAsFunction(FileO)();
    Print("We have RecResult\n");
    GRP_U:=RecResult.GAPperm;
    Print("|G_rat|=", Order(GRP_rat.GAPperm), " |GRP_U|=", Order(GRP_U), " |GRP_V|=", Order(GRP_V.GAPperm), " |LCos|=", Length(RecResult.DoubleCosetsPerm), "\n");
    ListSets:=[];
    sumElt:=0;
    if GRP_U<>GRP_V.GAPperm then
        Print("GRP_U and GRP_V should be equal\n");
        return false;
    fi;
    for eCos in RecResult.DoubleCosetsPerm
    do
        eList:=[];
        for eU in GRP_U
        do
            for eV in GRP_V.GAPperm
            do
                eG:=eU * eCos * eV;
                Add(eList, eG);
            od;
        od;
        eSet:=Set(eList);
        sumElt:=sumElt + Length(eSet);
        Add(ListSets, eSet);
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
    if sumElt<>Order(GRP_rat.GAPperm) then
        Print("The union of all cosets is not the full group\n");
        return false;
    fi;
    Print("All checks done\n");
    return true;
end;

#

CI_Decision_Reset();
pos:=0;
n_error:=0;

for EXT in ListEXT
do
    Print("     pos=", pos, "/", Length(ListEXT), " |EXT|=", Length(EXT), "/", Length(EXT[1]), "\n");
    test:=TestCase_Automorphy_DoubleCoset(EXT);
    if test=false then
        n_error:=n_error + 1;
    fi;
    pos:=pos + 1;
od;

if n_error=0 then
    Print("Normal case\n");
    CI_Write_Ok();
else
    Print("Error case\n");
fi;
