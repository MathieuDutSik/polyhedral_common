Read("../common.g");

ListFiles:=[];
for n in [5..7]
do
    FileName:=Concatenation("ClassificationSimplices", String(n));
    Add(ListFiles, FileName);
od;

TestCase:=function(EXT)
    local FileI, FileO, arith, OutFormat, eProg, TheCommand, TheGRP;
    FileI:="Test.in";
    FileO:="Test.out";
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    WriteMatrixFile(FileI, EXT);
    arith:="rational";
    OutFormat:="GAP";
    eProg:="../../src_group/GRP_LinPolytopeIntegral_Automorphism";
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " ", OutFormat, " ", FileO);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    TheGRP:=ReadAsFunction(FileO)();
    Print("|TheGRP|=", Order(TheGRP), "\n");
end;




for eFile in ListFiles
do
    ListEXT:=ReadAsFunction(eFile)();
    Print("eFile=", eFile, " |ListEXT|=", Length(ListEXT), "\n");
    for EXT in ListEXT
    do
        TestCase(EXT);
    od;
od;
