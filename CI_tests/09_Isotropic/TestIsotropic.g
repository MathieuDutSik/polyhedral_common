Read("../common.g");
Print("Beginning TestIsotropic\n");

OnlyTestExistence:=false;

TestIsotropic:=function(eRec)
    local n, FileIn, FileOut, eProg, TheCommand, U, V, eNorm, eNormSqr;
    n:=Length(eRec.M);
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    RemoveFileIfExist(FileIn);
    RemoveFileIfExist(FileOut);
    #
    WriteMatrixFile(FileIn, eRec.M);
    #
    if OnlyTestExistence then
        eProg:="../../src_isotropy/LATT_TestIsotropic";
    else
        eProg:="../../src_isotropy/LATT_FindIsotropic";
    fi;
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
    if U.has_isotropic and OnlyTestExistence=false then
        V:=U.V;
        eNorm:=V * eRec.M * V;
        if eNorm<>0 then
            Print("The vector is not isotropic\n");
            return rec(is_correct:=false);
        fi;
        eNormSqr:=V * V;
        if eNormSqr=0 then
            Print("The vector is zero\n");
            return rec(is_correct:=false);
        fi;
    fi;
    return rec(is_correct:=true);
end;

ListRec:=ReadAsFunction("../DATA/IsotropicCases")();;
ListDim3:=Filtered(ListRec, x->Length(x.M) = 3);
ListDim4:=Filtered(ListRec, x->Length(x.M) = 4);
ListDim5:=Filtered(ListRec, x->Length(x.M) = 5);
ListDim5gt:=Filtered(ListRec, x->Length(x.M) > 5);
#ListRec:=Concatenation(ListDim3_B, ListDim4_B, ListDim5high);
#ListRec:=Concatenation(ListDim3_B, ListDim5high);


#ListRec:=ListDim3;
#ListRec:=ListDim4;
#ListRec:=ListDim5;
#ListRec:=ListDim5gt;
#ListRec:=ListDim3_A{[3271..3272]};
#ListRec:=ListDim4_A;
#ListRec:=ListDim5;


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
Print("n_error=", n_error, "\n");
CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

