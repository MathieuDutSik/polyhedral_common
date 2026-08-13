Read("../common.g");
Print("Beginning TestPeriodicDelaunay\n");

# Exercise of the periodic point set support of src_delaunay.
#
# The five TEST binaries are self-checking: each one confronts the periodic
# code with the defining properties of the objects (empty circumspheres,
# direct enumeration beyond a wall, hand-computed group orders) and exits
# nonzero on any failure. LATT_SerialPeriodic_IsoDelaunayDomain is then run
# on the point set Z^2 + {(0,0), (1/3,0)} in the full symmetric space, for
# which the number of orbits of periodic iso-Delaunay domains is 1, a value
# cross-checked against the GAP Periodic_EnumerationProcedureLtype.

# GAP's Exec drops the exit status; recover it through the shell.
Exec_GetReturnValue:=function(TheCommand)
    local TmpDir, FileRC, list_lines;
    TmpDir:=DirectoryTemporary();
    FileRC:=Filename(TmpDir, "rc");
    Exec(Concatenation("( ", TheCommand, " ) ; echo $? > ", FileRC));
    list_lines:=ReadTextFile(FileRC);
    if Length(list_lines)=0 then
        return 1;
    fi;
    return Int(list_lines[1]);
end;

RunSelfCheckingTest:=function(eProgName)
    local eProg, TmpDir, FileErr, TheCommand, TheResult;
    eProg:=GetBinaryFilename(eProgName);
    TmpDir:=DirectoryTemporary();
    FileErr:=Filename(TmpDir, "test.err");
    TheCommand:=Concatenation(eProg, " 2> ", FileErr);
    TheResult:=Exec_GetReturnValue(TheCommand);
    if TheResult<>0 then
        Print("The test ", eProgName, " failed, its output:\n");
        Exec(Concatenation("cat ", FileErr));
        return 1;
    fi;
    Print("The test ", eProgName, " passed\n");
    return 0;
end;

RunPeriodicIsoDelaunay:=function()
    local eProg, TmpDir, FileCosets, FileNml, FileRes, FileErr, strNml, TheCommand, TheResult, TheRec;
    eProg:=GetBinaryFilename("LATT_SerialPeriodic_IsoDelaunayDomain");
    TmpDir:=DirectoryTemporary();
    FileCosets:=Filename(TmpDir, "Cosets");
    FileNml:=Filename(TmpDir, "Periodic.nml");
    FileRes:=Filename(TmpDir, "Result");
    FileErr:=Filename(TmpDir, "run.err");
    WriteMatrixFile(FileCosets, [ [0,0], [1/3,0] ]);
    strNml:="&SYSTEM\n";
    strNml:=Concatenation(strNml, " max_runtime_second = 0\n");
    strNml:=Concatenation(strNml, " ApplyStdUnitbuf = T\n");
    strNml:=Concatenation(strNml, " Saving = F\n");
    strNml:=Concatenation(strNml, " Prefix = \"/irrelevant/\"\n");
    strNml:=Concatenation(strNml, " OutFile = \"", FileRes, "\"\n");
    strNml:=Concatenation(strNml, " OutFormat = \"NumberGAP\"\n");
    strNml:=Concatenation(strNml, "/\n\n");
    strNml:=Concatenation(strNml, "&DATA\n");
    strNml:=Concatenation(strNml, " arithmetic = \"gmp\"\n");
    strNml:=Concatenation(strNml, " FileDualDescription = \"unset\"\n");
    strNml:=Concatenation(strNml, " FileCosets = \"", FileCosets, "\"\n");
    strNml:=Concatenation(strNml, "/\n\n");
    strNml:=Concatenation(strNml, "&TSPACE\n");
    strNml:=Concatenation(strNml, " TypeTspace = \"Classic\"\n");
    strNml:=Concatenation(strNml, " ClassicDim = 2\n");
    strNml:=Concatenation(strNml, "/\n");
    WriteStringFile(FileNml, strNml);
    TheCommand:=Concatenation(eProg, " ", FileNml, " 2> ", FileErr);
    TheResult:=Exec_GetReturnValue(TheCommand);
    if TheResult<>0 then
        Print("LATT_SerialPeriodic_IsoDelaunayDomain failed, its output:\n");
        Exec(Concatenation("cat ", FileErr));
        return 1;
    fi;
    TheRec:=ReadAsFunction(FileRes)();
    Print("Periodic iso-Delaunay domains: nb=", TheRec.nb, "\n");
    if TheRec.nb<>1 then
        Print("The number of orbits should be 1, the GAP-verified value\n");
        return 1;
    fi;
    return 0;
end;

# The A4 + cosets family (c generating A4*/A4 = Z/5 in the A4 basis) and
# HE6 = E6 + {0, c6} (c6 of order 3 in E6*/E6). The
# expected orbit data (vertex count, stabilizer order) was computed
# independently by the GAP Periodic_DelaunayDescriptionStandard and by
# LATT_SerialPeriodicDelaunay, in exact agreement (August 2026).
RunA4Family:=function()
    local eProg, TmpDir, GramA4, c, c2, ListCase, n_error, eCase, FileG, FileC, FileNml, FileRes, FileErr, strNml, TheCommand, TheResult, TheRec, ListSiz;
    eProg:=GetBinaryFilename("LATT_SerialPeriodicDelaunay");
    TmpDir:=DirectoryTemporary();
    GramA4:=[[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,2]];
    c :=[4/5, 3/5, 2/5, 1/5];
    c2:=[3/5, 1/5, 4/5, 2/5];
    GramE6:=[[2,-1,0,0,0,0],[-1,2,-1,0,0,0],[0,-1,2,-1,0,-1],
             [0,0,-1,2,-1,0],[0,0,0,-1,2,0],[0,0,-1,0,0,2]];
    c6:=[1/3, 2/3, 0, 1/3, 2/3, 0];
    ListCase:=[
      rec(name:="A4+c",   gram:=GramA4, cosets:=[[0,0,0,0], c],     expected:=[ [5,120], [20,240] ]),
      rec(name:="A4+2c",  gram:=GramA4, cosets:=[[0,0,0,0], c2],    expected:=[ [5,120], [8,48], [10,240] ]),
      rec(name:="A4+c2c", gram:=GramA4, cosets:=[[0,0,0,0], c, c2], expected:=[ [5,8], [5,120], [8,48] ]),
      rec(name:="HE6",    gram:=GramE6, cosets:=[[0,0,0,0,0,0], c6], expected:=[ [54,103680] ])
    ];
    n_error:=0;
    for eCase in ListCase
    do
        # The file writers of common.g append, so the names are made
        # unique per case.
        FileG:=Filename(TmpDir, Concatenation("Gram_", eCase.name));
        FileC:=Filename(TmpDir, Concatenation("Cosets_", eCase.name));
        FileNml:=Filename(TmpDir, Concatenation("Run_", eCase.name, ".nml"));
        FileRes:=Filename(TmpDir, Concatenation("Result_", eCase.name));
        FileErr:=Filename(TmpDir, Concatenation("run_", eCase.name, ".err"));
        WriteMatrixFile(FileG, eCase.gram);
        WriteMatrixFile(FileC, eCase.cosets);
        strNml:="&SYSTEM\n";
        strNml:=Concatenation(strNml, " max_runtime_second = 0\n");
        strNml:=Concatenation(strNml, " ApplyStdUnitbuf = T\n");
        strNml:=Concatenation(strNml, " Saving = F\n");
        strNml:=Concatenation(strNml, " Prefix = \"/irrelevant/\"\n");
        strNml:=Concatenation(strNml, " OutFile = \"", FileRes, "\"\n");
        strNml:=Concatenation(strNml, " OutFormat = \"SummaryGAP\"\n");
        strNml:=Concatenation(strNml, "/\n\n&DATA\n");
        strNml:=Concatenation(strNml, " arithmetic = \"gmp\"\n");
        strNml:=Concatenation(strNml, " GRAMfile = \"", FileG, "\"\n");
        strNml:=Concatenation(strNml, " FileCosets = \"", FileC, "\"\n");
        strNml:=Concatenation(strNml, "/\n");
        WriteStringFile(FileNml, strNml);
        TheCommand:=Concatenation(eProg, " ", FileNml, " 2> ", FileErr);
        TheResult:=Exec_GetReturnValue(TheCommand);
        if TheResult<>0 then
            Print("LATT_SerialPeriodicDelaunay failed on ", eCase.name, "\n");
            Exec(Concatenation("cat ", FileErr));
            n_error:=n_error+1;
        else
            TheRec:=ReadAsFunction(FileRes)();
            ListSiz:=Set(List(TheRec.ListRec, x->[x.nVert, x.ordStab]));
            Print(eCase.name, ": nb=", TheRec.nb, " orbits=", ListSiz, "\n");
            if ListSiz<>Set(eCase.expected) then
                Print("The orbit data differs from the GAP-verified value\n");
                n_error:=n_error+1;
            fi;
        fi;
    od;
    return n_error;
end;

FullTest:=function()
    local n_error, ListTest, eTest;
    n_error:=0;
    ListTest:=["TEST_PeriodicCVP", "TEST_PeriodicCosetSubgroup",
               "TEST_PeriodicFlipping", "TEST_PeriodicTspaceStabEqui",
               "TEST_PeriodicIsoDelaunay"];
    for eTest in ListTest
    do
        n_error:=n_error + RunSelfCheckingTest(eTest);
    od;
    n_error:=n_error + RunPeriodicIsoDelaunay();
    n_error:=n_error + RunA4Family();
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
