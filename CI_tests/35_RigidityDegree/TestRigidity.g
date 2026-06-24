Read("../common.g");
Print("Beginning TestRigidityDegree\n");

# Rigidity degree of a lattice, computed as a QUERIES option of
# LATT_SerialComputeDelaunay (FileRigidityDegree). Mirrors GAP's
# MyPolyhedral/lib/LatticeDelaunays.g::GetRigidityDegree.
TestRigidity:=function(eRec)
    local TmpDir, FileG, FileN, FileO, FileE, strOut, eProg, TheCommand,
          obtained, is_correct;
    TmpDir:=DirectoryTemporary();
    FileG:=Filename(TmpDir, "Gram.in");
    FileN:=Filename(TmpDir, "Rigid.nml");
    FileO:=Filename(TmpDir, "Rigid.out");
    FileE:=Filename(TmpDir, "Rigid.err");
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
    strOut:=Concatenation(strOut, " FileRigidityDegree = \"", FileO, "\"\n");
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
    obtained:=ReadAsFunction(FileO)();
    RemoveFile(FileG);
    RemoveFile(FileN);
    RemoveFile(FileO);
    is_correct:=obtained=eRec.rigidity;
    Print("name=", eRec.name, " obtained=", obtained,
          " expected=", eRec.rigidity, " is_correct=", is_correct, "\n");
    return rec(is_correct:=is_correct);
end;

ListRec:=[
rec(name:="Z2", eG:=[ [ 2, 0 ], [ 0, 2 ] ], rigidity:=2),
rec(name:="A2", eG:=[ [ 2, 1 ], [ 1, 2 ] ], rigidity:=3),
rec(name:="Z3", eG:=[ [ 2, 0, 0 ], [ 0, 2, 0 ], [ 0, 0, 2 ] ], rigidity:=3),
rec(name:="D4",
    eG:=[ [ 2, -1, 0, 0 ], [ -1, 2, -1, -1 ], [ 0, -1, 2, 0 ],
          [ 0, -1, 0, 2 ] ],
    rigidity:=1),
rec(name:="D5",
    eG:=[ [ 2, -1, 0, 0, 0 ], [ -1, 2, -1, 0, 0 ], [ 0, -1, 2, -1, -1 ],
          [ 0, 0, -1, 2, 0 ], [ 0, 0, -1, 0, 2 ] ],
    rigidity:=1),
rec(name:="E6",
    eG:=[ [ 2, -1, 0, 0, 0, 0 ], [ -1, 2, -1, 0, 0, 0 ],
          [ 0, -1, 2, -1, 0, -1 ], [ 0, 0, -1, 2, -1, 0 ],
          [ 0, 0, 0, -1, 2, 0 ], [ 0, 0, -1, 0, 0, 2 ] ],
    rigidity:=1)
];

FullTest:=function()
    local n_error, iRec, eRec, RecReply;
    n_error:=0;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), " name=", eRec.name, "\n");
        RecReply:=TestRigidity(eRec);
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
