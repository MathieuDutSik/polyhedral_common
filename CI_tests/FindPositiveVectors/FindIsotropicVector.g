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
    Print("V=", V, " eNorm=", eNorm, " StrictIneq=", StrictIneq, " CritNorm=", CritNorm, "\n");
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
    test:=SingleTest(M, CritNorm, true);
    if test=false then
        return false;
    fi;
    test:=SingleTest(M, CritNorm, false);
    if test=false then
        return false;
    fi;
    return true;
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



#ListM:=ReadAsFunction("ListM1")();;
ListM:=ReadAsFunction("LGramReflect")();;

GetSqrMat:=function(eV)
    local dim, eM, pos, i, j;
    dim:=Sqrt(Length(eV));
    eM:=NullMat(dim, dim);
    pos:=0;
    for i in [1..dim]
    do
        for j in [1..dim]
        do
            pos:=pos+1;
            eM[i][j]:=eV[pos];
        od;
    od;
    return eM;
end;




FullTest:=function()
    local iRec, eV, eM, reply;
    iRec:=0;
    for eV in ListM
    do
        Print("iRec=", iRec, " / ", Length(ListM), "\n");
        eM:=GetSqrMat(eV);
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

