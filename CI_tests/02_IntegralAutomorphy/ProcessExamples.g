Read("../common.g");

ListFiles:=[];
for n in [5..7]
do
    FileName:=Concatenation("../DATA/ClassificationSimplices", String(n));
    Add(ListFiles, FileName);
od;

TestCase_Automorphy:=function(EXT)
    local FileI, FileO, arith, eProg, TheCommand, TheGRP;
    FileI:=Filename(DirectoryTemporary(), "Test.in");
    FileO:=Filename(DirectoryTemporary(), "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:="../../src_group/GRP_LinPolytopeIntegral_Automorphism";
    TheCommand:=Concatenation(eProg, " rational ", FileI, " GAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("FileO does not exist\n");
        return false;
    fi;
    TheGRP:=ReadAsFunction(FileO)();
    Print("|TheGRP|=", Order(TheGRP), "\n");
    return true;
end;

TestCase_Automorphy_RightCoset:=function(EXT)
    local FileI, FileO, arith, OutFormat, eProg, TheCommand, TheGRP;
    FileI:=Filename(DirectoryTemporary(), "Test.in");
    FileO:=Filename(DirectoryTemporary(), "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:="../../src_group/GRP_LinPolytopeIntegral_Automorphism_RightCoset";
    TheCommand:=Concatenation(eProg, " rational ", FileI, " GAP ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("FileO does not exist\n");
        return false;
    fi;
    TheGRP:=ReadAsFunction(FileO)();
    Print("|TheGRP|=", Order(TheGRP), "\n");
    return true;
end;

#

n_error:=0;
for eFile in ListFiles
do
    ListEXT:=ReadAsFunction(eFile)();
    Print("eFile=", eFile, " |ListEXT|=", Length(ListEXT), "\n");
    for EXT in ListEXT
    do
        test:=TestCase_Automorphy(EXT);
        if test=false then
            n_error:=n_error+1;
        fi;
        test:=TestCase_Automorphy_RightCoset(EXT);
        if test=false then
            n_error:=n_error+1;
        fi;
    od;
od;

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
