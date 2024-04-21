Read("../common.g");
Print("Beginning TestIsotropic\n");

TestIsotropic:=function(eRec)
    local n, FileIn, FileOut, eProg, TheCommand, U, V, eNorm;
    n:=Length(eRec.M);
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    RemoveFileIfExist(FileIn);
    RemoveFileIfExist(FileOut);
    #
    WriteMatrixFile(FileIn, eRec.M);
    #
    eProg:="../../src_indefinite/LATT_FindIsotropic";
    TheCommand:=Concatenation(eProg, " rational ", FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileOut);
    if eRec.has_isotropic<>U.has_isotropic then
        Print("The result of has_isotropic is inconsistent\n");
        return rec(is_correct:=false);
    fi;
    if U.has_isotropic then
        V:=U.V;
        eNorm:=V * eRec.M * V;
        if eNorm<>0 then
            Print("The vector is not isotropic\n");
            return rec(is_correct:=false);
        fi;
    fi;
    return rec(is_correct:=true);
end;

ListRec:=ReadAsFunction("IsotropicCases")();;


FullTest:=function()
    local n_error, iRec, eRec, RecReply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        RecReply:=TestIsotropic(eRec);
        if RecReply.is_correct=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;

n_error:=FullTest();
if n_error > 0 then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

