TestLamination:=function(eRec)
    local FileOut, eProg, TheCommand, the_volume;
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    RemoveFileIfExist(FileOut);
    #
    eProg:="../../src_poly/POLY_TwoLaminations";
    TheCommand:=Concatenation(eProg, " rational one ", eRec.FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    answer:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    test1:=answer=fail;
    test2:=eRec.answer=fail;
    if test1<>test2 then
        Print("incoherence of the result\n");
        return false;
    fi;
    return true;
end;

eRec1:=rec(FileIn:="G6.ext", answer:=fail);
eRec2:=rec(FileIn:="24cell_poly.ext", answer:=fail);
eRec3:=rec(FileIn:="H3.ext", answer:=[]);
ListRec:=[eRec1, eRec2, eRec3];

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
if test=false then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

