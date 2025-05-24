Read("../common.g");
Print("Beginning TestLamination\n");

TestLamination:=function(eRec)
    local FileOut, eProg, TheCommand, answer, test1, test2;
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    #
    eProg:="../../src_poly/POLY_TwoLaminations";
    TheCommand:=Concatenation(eProg, " rational one ", eRec.FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    answer:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    test1:=answer=fail;
    test2:=eRec.answer=fail;
    if test1<>test2 then
        Print("incoherence of the result\n");
        return false;
    fi;
    return true;
end;

ListRec:=[];
Add(ListRec, rec(FileIn:="G6.ext", answer:=fail));
Add(ListRec, rec(FileIn:="24cell_poly.ext", answer:=fail));
Add(ListRec, rec(FileIn:="H3.ext", answer:=[]));

FullTest:=function()
    local iRec, eRec, test;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        test:=TestLamination(eRec);
        if test=false then
            return false;
        fi;
        iRec:=iRec + 1;
    od;
    return true;
end;

test:=FullTest();
Print("test=", test, "\n");

CI_Decision_Reset();
if test=false then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

