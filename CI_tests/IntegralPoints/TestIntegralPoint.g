Read("../common.g");
Print("Beginning TestIntegralPoint\n");

TestIntegralPoint:=function(FileEXT)
    local FileOut, FileFAC, eProg1, command, TheCommand1, eProg2, TheCommand2, EXT, answer;
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    FileFAC:=Filename(DirectoryTemporary(), "Test.fac");
    #
    eProg1:="../../src_poly/POLY_dual_description";
    command:="cdd";
    TheCommand1:=Concatenation(eProg1, " rational  ", command, " CPP ", FileEXT, " ", FileFAC);
    Exec(TheCommand1);
    if IsExistingFile(FileFAC)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    #
    eProg2:="../../src_poly/POLY_IntegralPoints";
    TheCommand2:=Concatenation(eProg2, " ", FileFAC, " ", FileEXT, " ", FileOut);
    Exec(TheCommand2);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    EXT:=ReadMatrixFile(FileEXT);
    answer:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    RemoveFile(FileFAC);
    if Length(EXT)<>Length(answer) then
        Print("For the Delaunay polytopes, the set of integer points should be equal to the vertices");
        return false;
    fi;
    return true;
end;

ListFile:=["ER35.ext", "G6.ext"];
for i in [1..27]
do
    eFile:=Concatenation("Perfect8_", String(i));
    Add(ListFile, eFile);
od;

FullTest:=function()
    local iFile, eFile, test;
    iFile:=0;
    for eFile in ListFile
    do
        Print("iFile=", iFile, " / ", Length(ListFile), "\n");
        test:=TestIntegralPoint(eFile);
        if test=false then
            return false;
        fi;
        iFile:=iFile + 1;
    od;
    return true;
end;

test:=FullTest();

CI_Decision_Reset();
if test=false then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

