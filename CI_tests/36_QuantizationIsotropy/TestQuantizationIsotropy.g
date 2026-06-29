Read("../common.g");
Print("Beginning TestQuantizationIsotropy\n");

# Isotropy (extremal / "white" quantizer) test, computed as a QUERIES option of
# LATT_SerialComputeDelaunay (FileIsotropy). A lattice is isotropic when the
# second-moment matrix M of its Voronoi cell is proportional to the inverse Gram
# matrix, equivalently GramMat * M = Lambda * I for a scalar Lambda. The test is
# exact: the defect DefectMat = GramMat*M - Lambda*I, with Lambda the average of
# the diagonal of GramMat*M, must be the zero matrix. We check the boolean and
# its internal consistency with the defect matrix.
TestIsotropy:=function(eRec)
    local TmpDir, FileG, FileN, FileO, FileE, strOut, eProg, TheCommand, U, n,
          defect_is_zero, is_correct;
    TmpDir:=DirectoryTemporary();
    FileG:=Filename(TmpDir, "Gram.in");
    FileN:=Filename(TmpDir, "Isotropy.nml");
    FileO:=Filename(TmpDir, "Isotropy.out");
    FileE:=Filename(TmpDir, "Isotropy.err");
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
    strOut:=Concatenation(strOut, " FileIsotropy = \"", FileO, "\"\n");
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
    n:=Length(eRec.eG);
    # The boolean must agree with the defect being exactly zero, and with the
    # expected isotropy of the lattice.
    defect_is_zero:=U.DefectMat=NullMat(n, n);
    is_correct:=U.IsIsotropic=defect_is_zero
                and U.IsIsotropic=eRec.isIsotropic;
    Print("name=", eRec.name, " is_correct=", is_correct, "\n");
    Print("  IsIsotropic=", U.IsIsotropic, " Lambda=", U.Lambda, "\n");
    return rec(is_correct:=is_correct);
end;

ListRec:=[
# A3 (FCC / body-centered cubic dual) is an isotropic quantizer.
rec(name:="A3", eG:=[ [ 2, -1, 0 ], [ -1, 2, -1 ], [ 0, -1, 2 ] ],
    isIsotropic:=true),
# E6 is an isotropic quantizer.
rec(name:="E6",
    eG:=[ [ 2, 1, 1, 0, 1, 1 ], [ 1, 4, 1, 1, 1, 3 ], [ 1, 1, 2, 1, 1, 1 ],
          [ 0, 1, 1, 2, 1, 2 ], [ 1, 1, 1, 1, 2, 2 ], [ 1, 3, 1, 2, 2, 4 ] ],
    isIsotropic:=true),
# A rectangular lattice with unequal sides is not isotropic.
rec(name:="rect12", eG:=[ [ 1, 0 ], [ 0, 2 ] ],
    isIsotropic:=false)
];

FullTest:=function()
    local n_error, iRec, eRec, RecReply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), " name=", eRec.name, "\n");
        RecReply:=TestIsotropy(eRec);
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
