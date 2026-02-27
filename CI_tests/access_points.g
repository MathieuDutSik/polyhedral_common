GetBinaryFilename:=function(FileName)
    local TmpDir, TmpFile, eProg, path, TmpFileB, command;
    # If accessible, use the path;
    TmpDir:=DirectoryTemporary();
    TmpFile:=Filename(TmpDir, "Test.in");
    Exec("which ", FileName, " > ", TmpFile);
    read_text_file:=function(TheFile)
        local list_lines;
        list_lines:=ReadTextFile(TheFile);
        if Length(list_lines)>1 then
            Print("list_lines=", list_lines, "\n");
            Error("Two programs available, odd case in my opinion");
        fi;
        if Length(list_lines)=1 then
            return list_lines[1];
        fi;
        return fail;
    end;
    eProg:=read_text_file(TmpFile);
    if eProg<fail then
        return eProg;
    fi;
    # Querying the environment variable
    if IsBound(GAPInfo.SystemEnvironment.POLYHEDRAL_COMMON_SOURCE_CODE) then
        path:=GAPInfo.SystemEnvironment.POLYHEDRAL_COMMON_SOURCE_CODE;
        TmpFileB:=Filename(TmpDir, "TestB.in");
        command:=Concatenation("(cd ", path, " && find . -name ", FileName, " 2> ", TmpFileB, ")");
        Exec(command);
        eProg:=read_text_file(TmpFileB);
        if eProg<>fail then
            return eProg;
        fi;
    fi;
    return fail;
end;


INDEF_FORM_Stabilizer:=function(eGram)
    local TmpDir, FileI, FileO, FileE, eProg, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=GetFreeFile("Err_indef_form_autom_");
    eProg:=GetBinaryFilename("INDEF_FORM_AutomorphismGroup");
    WriteMatrixFile(FileIn, eGram);
    TheCommand:=Concatenation(eProg, " gmp ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: INDEF_FORM_Stabilizer did not return anything, likely crash";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return Group(U);
end;

INDEF_FORM_TestEquivalence:=function(eGram1, eGram2)
    local TmpDir, eProg, FileM1, FileM2, FileO, FileE, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    eProg:=GetBinaryFileName("INDEF_FORM_TestEquivalence");
    FileM1:=Filename(TmpDir, "Mat1.in");
    FileM2:=Filename(TmpDir, "Mat2.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileM1, eGram1);
    WriteMatrixFile(FileM2, eGram2);
    TheCommand:=Concatenation(eProg, " gmp ", FileM1, " ", FileM2, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: Nothing returned by INDEF_FORM_TestEquivalence, maybe a crash";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileM1);
    RemoveFile(FileM2);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;


INDEF_FORM_GetOrbitRepresentative:=function(eGram, Xnorm)
    local TmpDir, eProg, FileM1, FileM2, FileO, FileE, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    eProg:=GetBinaryFileName("INDEF_FORM_GetOrbitRepresentative");
    FileM:=Filename(TmpDir, "Mat.in");
    FileOut:=Filename(TmpDir, "Test.out");
    FileErr:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileM, eGram);
    TheCommand:=Concatenation(eProg, " gmp ", FileM, " ", String(Xnorm), " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: Nothing returned by INDEF_FORM_GetOrbitRepresentative, likely crash or something";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileM);
    RemoveFile(FileO);
    return U;
end;

INDEF_FORM_GetOrbitIsotropic:=function(eGram, kDim, choice)
    local TmpDir, eProg, FileM, FileO, FileE, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    eProg:=GetBinaryFileName("INDEF_FORM_GetOrbit_IsotropicKplane");
    FileM:=Filename(TmpDir, "Mat.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileM, eGram);
    TheCommand:=Concatenation(eProg, " gmp ", FileM, " ", String(kDim), " ", choice, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: Nothing returned by INDEF_FORM_GetOrbitIsotropic, likely crash or something";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileM);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;





GetReflectivityInformation:=function(LorMat)
    local n, TmpDir, FileI, FileN, FileO, FileE, strOut, eProg, TheCommand, U;
    n:=Length(LorMat);
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileN:=Filename(TmpDir, "Test.nml");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    #
    WriteMatrixFile(FileIn, LorMat);
    #
    strOut:="&PROC\n";
    strOut:=Concatenation(strOut, " FileLorMat = \"", FileIn, "\"\n");
    strOut:=Concatenation(strOut, " OptionInitialVertex = \"isotropic_vinberg\"\n");
    strOut:=Concatenation(strOut, " OutFormat = \"GAP\"\n");
    strOut:=Concatenation(strOut, " FileOut = \"", FileOut, "\"\n");
    strOut:=Concatenation(strOut, " OptionNorms = \"all\"\n");
    strOut:=Concatenation(strOut, " EarlyTerminationIfNotReflective = T\n");
    strOut:=Concatenation(strOut, " ComputeAllSimpleRoots = T\n");
    strOut:=Concatenation(strOut, "/\n");
    #
    WriteStringFile(FileNml, strOut);
    #
    eProg:=GetBinaryFilename("LORENTZ_FundDomain_AllcockEdgewalk");
    TheCommand:=Concatenation(eProg, " ", FileN, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: Nothing returned by LORENTZ_FundDomain_AllcockEdgewalk, likely crash or something";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileN);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

