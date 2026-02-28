get_grp_automorphy:=function(EXT)
    local TmpDir, FileI, FileO, arith, eProg, TheCommand, TheGRP_rat;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Automorphism");
    TheCommand:=Concatenation(eProg, " rational ", FileI, " RecGAP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytope_Automorphism did not return anything, likely crash";
    fi;
    TheGRP_rat:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    return TheGRP_rat;
end;

get_grp_integral_automorphy:=function(EXT)
    local TmpDir, FileI, FileO, arith, eProg, TheCommand, TheGRP_int;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Automorphism");
    TheCommand:=Concatenation(eProg, " rational ", FileI, " RecGAP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytopeIntegral_Automorphism did not return anything, likely crash";
    fi;
    TheGRP_int:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    return TheGRP_int;
end;

test_polytope_integral_equivalence:=function(EXT1, EXT2)
    local TmpDir, File1, File2, FileO, eProg, TheCommand, result;
    TmpDir:=DirectoryTemporary();
    File1:=Filename(TmpDir, "Test1.in");
    File2:=Filename(TmpDir, "Test2.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(File1, EXT1);
    WriteMatrixFile(File2, EXT2);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Isomorphism");
    TheCommand:=Concatenation(eProg, " mpz_class ", File1, " ", File2, " GAP ", FileO);
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: The GRP_LinPolytopeIntegral_Isomorphism has not created a file";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFileIfExist(File1);
    RemoveFileIfExist(File2);
    RemoveFileIfExist(FileO);
    return result;
end;

get_canonic_form:=function(EXT)
    local TmpDir, FileI, FileO, eProg, TheCommand, EXT_can;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Can.in");
    FileO:=Filename(TmpDir, "Can.out");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Canonic");
    TheCommand:=Concatenation(eProg, " mpq_class ", FileI, " GAP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: The GRP_LinPolytope_Canonic did not return anything";
    fi;
    EXT_can:=ReadAsFunction(FileO)();
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    return EXT_can;
end;

get_isomorphism_result:=function(EXT1, EXT2)
    local TmpDir, File1, File2, FileO, eProg, TheCommand, result;
    TmpDir:=DirectoryTemporary();
    File1:=Filename(TmpDir, "Test1.in");
    File2:=Filename(TmpDir, "Test2.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(File1, EXT1);
    WriteMatrixFile(File2, EXT2);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Isomorphism");
    TheCommand:=Concatenation(eProg, " rational ", File1, " ", File2, " GAP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: The GRP_LinPolytope_Isomorphism did not return anything";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFileIfExist(File1);
    RemoveFileIfExist(File2);
    RemoveFileIfExist(FileO);
    return result;
end;


get_saturated_space:=function(ListMat)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, ListMatRet;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Saturation.in");
    FileO:=Filename(TmpDir, "Saturation.out");
    FileE:=Filename(TmpDir, "Saturation.err");
    WriteListMatrixFile(FileI, ListMat);

    eProg:=GetBinaryFilename("TSPACE_IntegralSaturation");
    TheCommand:=Concatenation(eProg, " mpq_class ", FileI, " GAP ", FileO);
    Exec(TheCommand);

    ListMatRet:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListMatRet;
end;

get_latt_canonical_form:=function(eMat)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Can.in");
    FileO:=Filename(TmpDir, "Can.out");
    FileE:=Filename(TmpDir, "Can.err");
    WriteMatrixFile(FileI, eMat);
    #
    eProg:=GetBinaryFilename("LATT_Canonicalize");
    TheCommand:=Concatenation(eProg, " gmp ", FileI, " GAP_full ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: The LATT_Canonicalize has failed";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

get_latt_automorphism_group:=function(eMat)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Aut.in");
    FileO:=Filename(TmpDir, "Aut.out");
    FileE:=Filename(TmpDir, "Aut.err");
    WriteListMatrixFile(FileI, [eMat]);
    #
    eProg:=GetBinaryFilename("LATT_Automorphism");
    TheCommand:=Concatenation(eProg, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("get_latt_automorphism_group, no return file\n");
        return fail;
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

get_fullrank_invariant_family:=function(eMat, method)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Fullrank.in");
    FileO:=Filename(TmpDir, "Fullrank.out");
    FileE:=Filename(TmpDir, "Fullrank.err");
    WriteMatrixFile(FileI, eMat);
    eProg:=GetBinaryFilename("LATT_GenerateCharacteristicVectorSet");
    TheCommand:=Concatenation(eProg, " mpq_class mpz_class ", method, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    if IsExistingFile(FileO)=false then
        return "program failure: LATT_GenerateCharacteristicVectorSet did not return anything";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

get_latt_isomorphism_test:=function(eMat1, eMat2)
    local TmpDir, FileIn1, FileIn2, FileOut, FileErr, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileIn1:=Filename(TmpDir, "Isom1.in");
    FileIn2:=Filename(TmpDir, "Isom2.in");
    FileOut:=Filename(TmpDir, "Aut.out");
    FileErr:=Filename(TmpDir, "Aut.err");
    WriteListMatrixFile(FileIn1, [eMat1]);
    WriteListMatrixFile(FileIn2, [eMat2]);
    #
    eProg:=GetBinaryFilename("LATT_Isomorphism");
    TheCommand:=Concatenation(eProg, " ", FileIn1, " ", FileIn2, " GAP ", FileOut, " 2> ", FileErr);
    if IsExistingFile(FileOut)=false then
        Print("get_latt_isomorphism_test, no return file\n");
        return fail;
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn1);
    RemoveFile(FileIn2);
    RemoveFile(FileOut);
    RemoveFile(FileErr);
    return U;
end;




INDEF_FORM_Stabilizer:=function(eGram)
    local TmpDir, FileI, FileO, FileE, eProg, U, TheCommand;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=GetFreeFile("Err_indef_form_autom_");
    WriteMatrixFile(FileI, eGram);
    eProg:=GetBinaryFilename("INDEF_FORM_AutomorphismGroup");
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
    FileM1:=Filename(TmpDir, "Mat1.in");
    FileM2:=Filename(TmpDir, "Mat2.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileM1, eGram1);
    WriteMatrixFile(FileM2, eGram2);
    eProg:=GetBinaryFilename("INDEF_FORM_TestEquivalence");
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
    local TmpDir, eProg, FileM, FileO, FileE, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileM:=Filename(TmpDir, "Mat.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileM, eGram);
    eProg:=GetBinaryFilename("INDEF_FORM_GetOrbitRepresentative");
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
    FileM:=Filename(TmpDir, "Mat.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileM, eGram);
    eProg:=GetBinaryFilename("INDEF_FORM_GetOrbit_IsotropicKplane");
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
    WriteMatrixFile(FileI, LorMat);
    #
    strOut:="&PROC\n";
    strOut:=Concatenation(strOut, " FileLorMat = \"", FileI, "\"\n");
    strOut:=Concatenation(strOut, " OptionInitialVertex = \"isotropic_vinberg\"\n");
    strOut:=Concatenation(strOut, " OutFormat = \"GAP\"\n");
    strOut:=Concatenation(strOut, " FileOut = \"", FileO, "\"\n");
    strOut:=Concatenation(strOut, " OptionNorms = \"all\"\n");
    strOut:=Concatenation(strOut, " EarlyTerminationIfNotReflective = T\n");
    strOut:=Concatenation(strOut, " ComputeAllSimpleRoots = T\n");
    strOut:=Concatenation(strOut, "/\n");
    #
    WriteStringFile(FileN, strOut);
    #
    eProg:=GetBinaryFilename("LORENTZ_FundDomain_AllcockEdgewalk");
    TheCommand:=Concatenation(eProg, " ", FileN, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
#        Print(NullMat(5));
        return "program failure: Nothing returned by LORENTZ_FundDomain_AllcockEdgewalk, likely crash or something";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileN);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;


get_matrix_group_mod_information:=function(ListGen, p_val)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteListMatrixFile(FileI, ListGen);
    eProg:=GetBinaryFilename("LATT_ResolveModAction");
    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " ", String(p_val), " GAP ", FileO);
#    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " ", String(p_val), " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: Nothing returned by LATT_ResolveModAction, likely crash or something";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;


get_grp_size_matrix_group_mod_action:=function(ListGen, p_val)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteListMatrixFile(FileI, ListGen);
    eProg:=GetBinaryFilename("LATT_ComputeGroupModAction");
    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " ", String(p_val), " GAP ", FileO);
#    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " ", String(p_val), " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: Nothing returned by LATT_ResolveModAction, likely crash or something";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;
