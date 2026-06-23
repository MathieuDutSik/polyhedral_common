Read("../common.g");
Print("Beginning TestQuantizationDeformation\n");

# Deformation of the quantizer along Q + t H. We test the E6 root lattice with
# H = e1 e1^T, for which the exact answer is known:
#   SecMoment(t) = (15/28 + 5/6 t + 1/9 t^2) / (1 + 4/3 t)
#   SecMoment(0)=15/28, SecMoment'(0)=5/42, SecMoment''(0)=-2/21
TestDeformation:=function(eRec)
    local Q, H, TmpDir, FileQ, FileH, FileO, eProg, TheCommand, U, is_correct;
    Q:=eRec.Q;
    H:=eRec.H;
    TmpDir:=DirectoryTemporary();
    FileQ:=Filename(TmpDir, "Q.in");
    FileH:=Filename(TmpDir, "H.in");
    FileO:=Filename(TmpDir, "Out.g");
    WriteMatrixFile(FileQ, Q);
    WriteMatrixFile(FileH, H);
    eProg:=GetBinaryFilename("LATT_QuantizationDeformation");
    TheCommand:=Concatenation(eProg, " gmp ", FileQ, " ", FileH, " ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileQ);
    RemoveFile(FileH);
    RemoveFile(FileO);
    is_correct:=U.SecMoment0=eRec.SecMoment0
                and U.SecMoment1=eRec.SecMoment1
                and U.SecMoment2=eRec.SecMoment2
                and U.numerator=eRec.numerator
                and U.denominator=eRec.denominator;
    Print("name=", eRec.name, " is_correct=", is_correct, "\n");
    Print("  SecMoment0=", U.SecMoment0, " SecMoment1=", U.SecMoment1,
          " SecMoment2=", U.SecMoment2, "\n");
    return rec(is_correct:=is_correct);
end;

ListRec:=[
rec(name:="E6_e1",
    Q:=[ [ 2, 1, 1, 0, 1, 1 ], [ 1, 4, 1, 1, 1, 3 ], [ 1, 1, 2, 1, 1, 1 ],
         [ 0, 1, 1, 2, 1, 2 ], [ 1, 1, 1, 1, 2, 2 ], [ 1, 3, 1, 2, 2, 4 ] ],
    H:=[ [ 1, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0 ],
         [ 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0 ] ],
    SecMoment0:=15/28, SecMoment1:=5/42, SecMoment2:=-2/21,
    numerator:=[ 15/28, 5/6, 1/9 ],
    denominator:=[ 1, 4/3, 0, 0, 0, 0, 0 ])
];

FullTest:=function()
    local n_error, iRec, eRec, RecReply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), " name=", eRec.name, "\n");
        RecReply:=TestDeformation(eRec);
        if RecReply.is_correct=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    Print("FullTest: n_error=", n_error, "\n");
    return n_error;
end;

NestFunction:=function()
    local n_error;
    n_error:=FullTest();
    CI_Decision_Reset();
    if n_error > 0 then
        Print("Error case\n");
    else
        Print("Normal case\n");
        CI_Write_Ok();
    fi;
    CI_PrintExistConclusion();
end;

NestFunction();
