Read("../common.g");
Print("Beginning TestQuantizationHessian\n");

# Hessian of the normalized quantizer constant G at a lattice Q, computed as a
# QUERIES option of LATT_SerialComputeDelaunay (FileHessian). The method uses the
# derivative of the second-moment matrix M(Q) along rank-one directions v v^T
# (which span Sym^n by linearity of the moment derivative), assembled via the
# Aut(Q) orbit equivariance. We check the exact signature of the m=n(n+1)/2-1
# dimensional Hessian, together with the two internal consistency residuals
# (radial: beta Q = 0 at a critical point; cross-check: an independent scalar
# deformation reproduces the Hessian prediction on a held-out direction).
TestHessian:=function(eRec)
    local TmpDir, FileG, FileN, FileO, FileE, strOut, eProg, TheCommand, U,
          obtained, is_correct;
    TmpDir:=DirectoryTemporary();
    FileG:=Filename(TmpDir, "Gram.in");
    FileN:=Filename(TmpDir, "Hessian.nml");
    FileO:=Filename(TmpDir, "Hessian.out");
    FileE:=Filename(TmpDir, "Hessian.err");
    WriteMatrixFile(FileG, eRec.eG);
    #
    strOut:="&SYSTEM\n";
    strOut:=Concatenation(strOut, " OutFormat = \"nothing\"\n");
    strOut:=Concatenation(strOut, " OutFile = \"unset.out\"\n");
    strOut:=Concatenation(strOut, " max_runtime_second = 0\n");
    strOut:=Concatenation(strOut, "/\n\n");
    strOut:=Concatenation(strOut, "&DATA\n");
    strOut:=Concatenation(strOut, " arithmetic = \"gmp\"\n");
    strOut:=Concatenation(strOut, " GRAMfile = \"", FileG, "\"\n");
    strOut:=Concatenation(strOut, " SVRfile = \"unset.svr\"\n");
    strOut:=Concatenation(strOut, " CacheFile = \"none\"\n");
    strOut:=Concatenation(strOut, "/\n\n");
    strOut:=Concatenation(strOut, "&QUERIES\n");
    strOut:=Concatenation(strOut, " FileHessian = \"", FileO, "\"\n");
    strOut:=Concatenation(strOut, "/\n");
    WriteStringFile(FileN, strOut);
    #
    eProg:=GetBinaryFilename("LATT_SerialComputeDelaunay");
    TheCommand:=Concatenation(eProg, " ", FileN, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileG);
    RemoveFile(FileN);
    RemoveFile(FileO);
    obtained:=rec(nbPlus:=U.signature.nbPlus, nbMinus:=U.signature.nbMinus,
                  nbZero:=U.signature.nbZero);
    is_correct:=U.solved=true
                and U.radialResidual=0
                and U.checkResidual=0
                and obtained.nbPlus=eRec.nbPlus
                and obtained.nbMinus=eRec.nbMinus
                and obtained.nbZero=eRec.nbZero;
    Print("name=", eRec.name, " is_correct=", is_correct, "\n");
    Print("  signature=", obtained, " radialResidual=", U.radialResidual,
          " checkResidual=", U.checkResidual, " nbEval=", U.nbEval, "\n");
    return rec(is_correct:=is_correct);
end;

ListRec:=[
# A3 (FCC) is a degenerate critical point: 2 flat directions on the m=5 space.
rec(name:="A3", eG:=[ [ 2, -1, 0 ], [ -1, 2, -1 ], [ 0, -1, 2 ] ],
    nbPlus:=3, nbMinus:=0, nbZero:=2),
# E6 is a local minimum: positive definite Hessian on the m=20 space.
rec(name:="E6",
    eG:=[ [ 2, 1, 1, 0, 1, 1 ], [ 1, 4, 1, 1, 1, 3 ], [ 1, 1, 2, 1, 1, 1 ],
          [ 0, 1, 1, 2, 1, 2 ], [ 1, 1, 1, 1, 2, 2 ], [ 1, 3, 1, 2, 2, 4 ] ],
    nbPlus:=20, nbMinus:=0, nbZero:=0)
];

FullTest:=function()
    local n_error, iRec, eRec, RecReply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), " name=", eRec.name, "\n");
        RecReply:=TestHessian(eRec);
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
