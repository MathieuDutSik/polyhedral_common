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
    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " RecGAP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytopeIntegral_Automorphism did not return anything, likely crash";
    fi;
    TheGRP_int:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    return TheGRP_int;
end;

get_double_cosets_matrix_group:=function(EXT, GRPperm)
    local TmpDir, FileI, FileO, FileGRP_V, eProg, TheCommand, RecResult;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileGRP_V:=Filename(TmpDir, "Test.grp_V");
    WriteMatrixFile(FileI, EXT);
    WriteGroupFile(FileGRP_V, Length(EXT), GRPperm);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Automorphism_DoubleCoset");
    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " ", FileGRP_V, " RecGAP ", FileO);
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytopeIntegral_Automorphism_DoubleCoset should have created a file";
    fi;
    RecResult:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileGRP_V);
    return RecResult;
end;

get_linear_space_stabilizer_double_cosets:=function(GRPmatr, TheSpace, GRP_Vmatr)
    local TmpDir, FileGRP, FileSPA, FileGRP_V, FileO, eProg, TheCommand, RecResult;
    TmpDir:=DirectoryTemporary();
    FileGRP:=Filename(TmpDir, "Test.grp");
    FileSPA:=Filename(TmpDir, "Test.space");
    FileGRP_V:=Filename(TmpDir, "Test.grp_V");
    FileO:=Filename(TmpDir, "Test.out");
    WriteListMatrixFile(FileGRP, GRPmatr);
    WriteMatrixFile(FileSPA, TheSpace);
    WriteListMatrixFile(FileGRP_V, GRP_Vmatr);
    eProg:=GetBinaryFilename("GRP_LinearSpace_Stabilizer_DoubleCoset");

    TheCommand:=Concatenation(eProg, " ", FileGRP, " ", FileSPA, " ", FileGRP_V, " ", FileO);
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinearSpace_Stabilizer_DoubleCoset should have created a file";
    fi;
    RecResult:=ReadAsFunction(FileO)();
    RemoveFile(FileGRP);
    RemoveFile(FileSPA);
    RemoveFile(FileGRP_V);
    RemoveFile(FileO);
    return RecResult;
end;


get_linpolytopeintegral_aut_rightcoset:=function(EXT)
    local TmpDir, FileI, FileO, eProg, TheCommand, result;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Automorphism_RightCoset");
    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " GAP ", FileO);
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytopeIntegral_Automorphism_RightCoset should have created a file";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    return result;
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

get_fullrank_invariant_family:=function(eG, method)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Fullrank.in");
    FileO:=Filename(TmpDir, "Fullrank.out");
    FileE:=Filename(TmpDir, "Fullrank.err");
    WriteMatrixFile(FileI, eG);
    eProg:=GetBinaryFilename("LATT_GenerateCharacteristicVectorSet");
    TheCommand:=Concatenation(eProg, " mpq_class mpz_class ", method, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: LATT_GenerateCharacteristicVectorSet did not return anything";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

get_ext_volume:=function(EXT)
    local TmpDir, FileI, FileO, eProg, TheCommand, the_volume;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("POLY_lrs_volume");
    TheCommand:=Concatenation(eProg, " rational ", FileI, " GAP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_lrs_volume should have created a file";
    fi;
    the_volume:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    return the_volume;
end;

get_interior_point:=function(FAC)
    local TmpDir, FileI, FileO, eProg, TheCommand, TheV;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Interior.in");
    FileO:=Filename(TmpDir, "Interior.out");
    WriteMatrixFile(FileI, FAC);
    eProg:=GetBinaryFilename("POLY_GeometricallyUniqueInteriorPoint");
    TheCommand:=Concatenation(eProg, " rational ", FileI, " GAP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_GeometricallyUniqueInteriorPoint failed to create an output file";
    fi;
    TheV:=ReadAsFunction(FileO)();
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    return TheV;
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
        return "program failure: LATT_Isomorphism failed to create a file";
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn1);
    RemoveFile(FileIn2);
    RemoveFile(FileOut);
    RemoveFile(FileErr);
    return U;
end;



get_group_skeleton:=function(EXT, GRP, LevSearch)
    local TmpDir, FileEXT, FileGRP, FileInputNml, FileGrpOut, FileResultOut, strOut, output, eProg, TheCommand, GRPout, result;
    TmpDir:=DirectoryTemporary();
    FileEXT:=Filename(TmpDir, "Test.ext");
    FileGRP:=Filename(TmpDir, "Test.grp");
    FileInputNml:=Filename(TmpDir, "Test.nml");
    FileGrpOut:=Filename(TmpDir, "Test.grp");
    FileResultOut:=Filename(TmpDir, "Test.result");
    WriteMatrixFile(FileEXT, EXT);
    WriteGroupFile(FileGRP, Length(EXT), GRP);
    #
    strOut:="&PROC\n";
    strOut:=Concatenation(strOut, " FACfile = \"", FileEXT, "\"\n");
    strOut:=Concatenation(strOut, " EXTfile = \"unset.ext\"\n");
    strOut:=Concatenation(strOut, " GRPfile = \"", FileGRP, "\"\n");
    strOut:=Concatenation(strOut, " OutFile =\"", FileResultOut, "\"\n");
    strOut:=Concatenation(strOut, " ComputeTotalNumberFaces = T\n");
    strOut:=Concatenation(strOut, " method_spann = \"LinearProgramming\"\n");
    strOut:=Concatenation(strOut, " method_final = \"all\"\n");
    strOut:=Concatenation(strOut, " Arithmetic = \"rational\"\n");
    strOut:=Concatenation(strOut, " LevSearch = ", String(LevSearch), "\n");
    strOut:=Concatenation(strOut, "/\n");
    strOut:=Concatenation(strOut, "\n");
    strOut:=Concatenation(strOut, "&GROUP\n");
    strOut:=Concatenation(strOut, " ComputeAutGroup = T\n");
    strOut:=Concatenation(strOut, " OutFormat = \"GAP\"\n");
    strOut:=Concatenation(strOut, " FileGroup = \"", FileGrpOut, "\"\n");
    strOut:=Concatenation(strOut, "/\n");
    #
    WriteStringFile(FileInputNml, strOut);
    #
    eProg:=GetBinaryFilename("POLY_FaceLatticeGen");
    TheCommand:=Concatenation(eProg, " ", FileInputNml);
    Exec(TheCommand);
    if IsExistingFile(FileGrpOut)=false then
        return "program failure: POLY_FaceLatticeGen failed to create the group file";
    fi;
    if IsExistingFile(FileResultOut)=false then
        return "program failure: POLY_FaceLatticeGen failed to create the result file";
    fi;
    GRPout:=ReadAsFunction(FileGrpOut)();
    result:=ReadAsFunction(FileResultOut)();
    RemoveFile(FileEXT);
    RemoveFile(FileGRP);
    RemoveFile(FileInputNml);
    RemoveFile(FileGrpOut);
    RemoveFile(FileResultOut);
    return [GRPout, result];
end;


sample_facet_polytope:=function(EXT, method)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, result;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    #
    eProg:=GetBinaryFilename("POLY_sampling_facets");
    TheCommand:=Concatenation(eProg, " rational ", method, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_sampling_facets failed to create the output";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return result;
end;


test_matrix_group_finiteness:=function(GRPmatr)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Finiteness_test.input");
    FileO:=Filename(TmpDir, "Finiteness_test.result");
    FileE:=Filename(TmpDir, "Finiteness_test.err");
    WriteListMatrixFile(FileI, GRPmatr);
    #
    eProg:=GetBinaryFilename("GRP_TestFiniteness");
    TheCommand:=Concatenation(eProg, " gmp ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    #
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_TestFiniteness failed to create the output";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;



get_integral_interior_point:=function(FAC, method)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, EXTint;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.fac");
    FileO:=Filename(TmpDir, "Test.vertint");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, FAC);
    eProg:=GetBinaryFilename("POLY_IntegralPoints");
    TheCommand:=Concatenation(eProg, " gmp ", method, " ", FileI, " CPP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_IntegralPoints failed to create the file";
    fi;
    EXTint:=ReadMatrixFile(FileO);
    RemoveFile(FileI);
    RemoveFile(FileO);
    return EXTint;
end;

get_dual_desc:=function(EXT)
    local TmpDir, FileI, FileO, FileE, eProg, command, TheCommand, FAC;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.out");
    FileO:=Filename(TmpDir, "Test.fac");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("POLY_dual_description");
    command:="cdd";
    TheCommand:=Concatenation(eProg, " rational  ", command, " CPP ", FileI, " ", FileO, " 2>", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_dual_description failed to create the file";
    fi;
    FAC:=ReadMatrixFile(FileO);
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return FAC;
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
