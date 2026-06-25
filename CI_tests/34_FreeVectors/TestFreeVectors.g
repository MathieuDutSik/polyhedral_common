Read("../common.g");
Print("Beginning TestFreeVectors\n");

# Free vectors of a lattice, computed as a QUERIES option of
# LATT_SerialComputeDelaunay (FileFreeVectors). We test a few lattices for which
# the answer is known/stable.
TestFreeVectors:=function(eRec)
    local TmpDir, FileG, FileN, FileO, FileE, strOut, eProg, TheCommand, U,
          obtained, is_correct;
    TmpDir:=DirectoryTemporary();
    FileG:=Filename(TmpDir, "Gram.in");
    FileN:=Filename(TmpDir, "FreeVect.nml");
    FileO:=Filename(TmpDir, "FreeVect.out");
    FileE:=Filename(TmpDir, "FreeVect.err");
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
    strOut:=Concatenation(strOut, " FileFreeVectors = \"", FileO, "\"\n");
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
    # Stable signature: [nbRelevantVector, sorted [SubspaceDim,OrbitSize,nMatched]].
    obtained:=rec(nbRelevantVector:=U.nbRelevantVector,
                  nbFreeVectorOrbit:=U.nbFreeVectorOrbit,
                  Signature:=Set(List(U.ListFreeVectorOrbit,
                      x->[x.SubspaceDim, x.OrbitSize, x.nMatched])));
    is_correct:=obtained.nbRelevantVector=eRec.nbRelevantVector
                and obtained.nbFreeVectorOrbit=eRec.nbFreeVectorOrbit
                and obtained.Signature=eRec.Signature;
    Print("name=", eRec.name, " is_correct=", is_correct, "\n");
    Print("  obtained=", obtained, "\n");
    return rec(is_correct:=is_correct);
end;

ListRec:=[
rec(name:="A2", eG:=[ [ 2, 1 ], [ 1, 2 ] ],
    nbRelevantVector:=3, nbFreeVectorOrbit:=1,
    Signature:=Set([ [ 1, 3, 1 ] ])),
rec(name:="A4", eG:=[ [ 2, -1, 0, 0 ], [ -1, 2, -1, 0 ], [ 0, -1, 2, -1 ],
                      [ 0, 0, -1, 2 ] ],
    nbRelevantVector:=10, nbFreeVectorOrbit:=2,
    Signature:=Set([ [ 1, 5, 6 ], [ 1, 10, 4 ] ])),
rec(name:="E6",
    eG:=[ [ 2, 1, 1, 0, 1, 1 ], [ 1, 4, 1, 1, 1, 3 ], [ 1, 1, 2, 1, 1, 1 ],
          [ 0, 1, 1, 2, 1, 2 ], [ 1, 1, 1, 1, 2, 2 ], [ 1, 3, 1, 2, 2, 4 ] ],
    nbRelevantVector:=36, nbFreeVectorOrbit:=1,
    Signature:=Set([ [ 1, 27, 20 ] ]))
];

FullTest:=function()
    local n_error, iRec, eRec, RecReply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), " name=", eRec.name, "\n");
        RecReply:=TestFreeVectors(eRec);
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
