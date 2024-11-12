Read("../common.g");
CI_Decision_Reset();

ListFiles:=[];
for n in [5..7]
do
    FileName:=Concatenation("ClassificationSimplices", String(n));
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
    TheGRP:=ReadAsFunction(FileO)();
    Print("|TheGRP|=", Order(TheGRP), "\n");
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
    TheGRP:=ReadAsFunction(FileO)();
    Print("|TheGRP|=", Order(TheGRP), "\n");
end;

#

for eFile in ListFiles
do
    ListEXT:=ReadAsFunction(eFile)();
    Print("eFile=", eFile, " |ListEXT|=", Length(ListEXT), "\n");
    for EXT in ListEXT
    do
        TestCase_Automorphy(EXT);
        TestCase_Automorphy_RightCoset(EXT);
    od;
od;

CI_Write_Ok();
