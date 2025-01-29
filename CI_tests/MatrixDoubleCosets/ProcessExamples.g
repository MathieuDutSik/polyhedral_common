Read("../common.g");

TmpDir:=DirectoryTemporary();

ListEXT:=[];
for n in [5..7]
do
    FileName:=Concatenation("../DATA/ClassificationSimplices", String(n));
    PartListEXT:=ReadAsFunction(eFile)();
    Append(ListEXT, PartListEXT);
od;

get_grp_automorphy:=function(EXT)
    local FileI, FileO, arith, eProg, TheCommand, TheGRP_rat;
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:="../../src_group/GRP_LinPolytope_Automorphism";
    TheCommand:=Concatenation(eProg, " rational ", FileI, " GAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    TheGRP_rat:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    Print("|TheGRP_rat|=", Order(TheGRP_rat), "\n");
    return TheGRP_rat;
end;

get_grp_integral_automorphy:=function(EXT)
    local FileI, FileO, arith, eProg, TheCommand, TheGRP_int;
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:="../../src_group/GRP_LinPolytopeIntegral_Automorphism";
    TheCommand:=Concatenation(eProg, " rational ", FileI, " GAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    TheGRP_int:=ReadAsFunction(FileO)();
    Print("|TheGRP_int|=", Order(TheGRP_int), "\n");
    return TheGRP_int;
end;

# Compute the double cosets of G_rat = \cup_i G_int g_i G_int.

TestCase_Automorphy_DoubleCoset:=function(EXT)
    local FileI, FileO, arith, OutFormat, eProg, TheCommand, TheGRP;
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileI, EXT);
    GRP_rat:=get_grp_automorphy(EXT);
    GRP_V:=get_grp_integral_automorphy(EXT);
    eProg:="../../src_group/GRP_LinPolytopeIntegral_Automorphism_DoubleCoset";
    TheCommand:=Concatenation(eProg, " rational ", FileI, " RecGAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    RecResult:=ReadAsFunction(FileO)();
    GRP_U:=RecResult.GAPperm;
    ListSets:=[];
    sumElt:=0;
    for eCos in RecResult.DoubleCosetsPerm
    do
        eList:=[];
        for eU in GRP_U
        do
            for eV in GRP_V
            do
                eG:=eU * eCos * eV;
                Add(eList, eG);
            od;
        od;
        eSet:=Set(eList);
        sumElt:=sumElt + Length(eSet);
        Add(ListSets, eSet);
    od;
    for iSet in [1..Length(ListSets)]
    do
        for jSet in [iSet+1..Length(ListSet)]
        do
            eInt:=Intersection(ListSet[iSet], ListSet[jSet]);
            if Length(eInt) > 0 then
                Print("Intersection are not empty\n");
                return false;
            fi;
        od;
    od;
    if sumElt<>Order(GRP_rat) then
        Print("The union of all cosets is not the full group\n");
        return false;
    fi;
    return true;

end;

#

CI_Decision_Reset();
pos:=0;
n_error:=0;

for EXT in ListEXT
do
    Print("pos=", pos, "/", Length(ListEXT), " |EXT|=", Length(EXT), "/", Length(EXT[1]), "\n");
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
