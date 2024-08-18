Read("../common.g");
Print("Beginning FindIsotropicVector\n");

SingleTest:=function(M, CritNorm, StrictIneq)
    local FileIn, FileOut, eProg, TheCommand, V, eNorm;
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    RemoveFileIfExist(FileIn);
    RemoveFileIfExist(FileOut);
    #
    WriteMatrixFile(FileIn, M);
    #
    eProg:="../../src_isotropy/LATT_FindPositiveVector";
    TheCommand:=Concatenation(eProg, " gmp ", FileIn, " ", String(CritNorm), " ", String(StrictIneq), " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    V:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    eNorm:=V * M * V;
    if StrictIneq then
        if eNorm <= CritNorm then
            Print("1: eNorm=", eNorm, " CritNorm=", CritNorm, " StrictIneq=", StrictIneq, "\n");
            return false;
        fi;
    else
        if eNorm < CritNorm then
            Print("2: eNorm=", eNorm, " CritNorm=", CritNorm, " StrictIneq=", StrictIneq, "\n");
            return false;
        fi;
    fi;
    return true;
end;

SeveralTest_A:=function(M, CritNorm)
    local test;
    test:=SingleTest(M, 0, true);
    if test=false then
        return false;
    fi;
    test:=SingleTest(M, 0, false);
    if test=false then
        return false;
    fi;
end;

SeveralTest_B:=function(M)
    local CritNorm, test;
    for CritNorm in [-1,0,1]
    do
        test:=SeveralTest_A(M, CritNorm);
        if test=false then
            return false;
        fi;
    od;
    return true;
end;

SeveralTest_C:=function(M)
    local test;
    test:=SeveralTest_B(M);
    if test=false then
        return false;
    fi;
    test:=SeveralTest_B(-M);
    if test=false then
        return false;
    fi;
    return true;
end;




ListM:=ReadAsFunction("IndefiniteForms")();;

FullTest:=function()
    local iRec, eM, reply;
    iRec:=0;
    for eM in ListM
    do
        Print("iRec=", iRec, " / ", Length(ListM), "\n");
        reply:=SeveralTest_C(eM);
        if reply=false then
            return false;
        fi;
        iRec:=iRec + 1;
    od;
    return true;
end;

is_correct:=FullTest();
if is_correct=false then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

