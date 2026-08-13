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
