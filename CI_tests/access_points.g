get_grp_automorphy:=function(EXT)
    local TmpDir, FileI, FileO, FileE, arith, eProg, TheCommand, TheGRP_rat;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Automorphism");
    TheCommand:=Concatenation(eProg, " rational ", FileI, " RecGAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytope_Automorphism did not return anything, likely crash";
    fi;
    TheGRP_rat:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return TheGRP_rat;
end;

get_grp_integral_automorphy:=function(EXT)
    local TmpDir, FileI, FileO, FileE, arith, eProg, TheCommand, TheGRP_int;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Automorphism");
    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " RecGAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytopeIntegral_Automorphism did not return anything, likely crash";
    fi;
    TheGRP_int:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return TheGRP_int;
end;

get_double_cosets_matrix_group:=function(EXT, GRPperm)
    local TmpDir, FileI, FileO, FileE, FileGRP_V, eProg, TheCommand, RecResult;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    FileGRP_V:=Filename(TmpDir, "Test.grp_V");
    WriteMatrixFile(FileI, EXT);
    WriteGroupFile(FileGRP_V, Length(EXT), GRPperm);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Automorphism_DoubleCoset");
    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " ", FileGRP_V, " RecGAP ", FileO, " 2> ", FileE);
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytopeIntegral_Automorphism_DoubleCoset should have created a file";
    fi;
    RecResult:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    RemoveFile(FileGRP_V);
    return RecResult;
end;

get_linear_space_stabilizer_double_cosets:=function(GRPmatr, TheSpace, GRP_Vmatr)
    local TmpDir, FileGRP, FileSPA, FileGRP_V, FileO, FileE, eProg, TheCommand, RecResult;
    TmpDir:=DirectoryTemporary();
    FileGRP:=Filename(TmpDir, "Test.grp");
    FileSPA:=Filename(TmpDir, "Test.space");
    FileGRP_V:=Filename(TmpDir, "Test.grp_V");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteListMatrixFile(FileGRP, GRPmatr);
    WriteMatrixFile(FileSPA, TheSpace);
    WriteListMatrixFile(FileGRP_V, GRP_Vmatr);
    eProg:=GetBinaryFilename("GRP_LinearSpace_Stabilizer_DoubleCoset");
    TheCommand:=Concatenation(eProg, " ", FileGRP, " ", FileSPA, " ", FileGRP_V, " ", FileO, " 2> ", FileE);
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
    RemoveFile(FileE);
    return RecResult;
end;


get_linpolytopeintegral_aut_rightcoset:=function(EXT)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, result;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Automorphism_RightCoset");
    TheCommand:=Concatenation(eProg, " mpz_class ", FileI, " GAP ", FileO, " 2> ", FileE);
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytopeIntegral_Automorphism_RightCoset should have created a file";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return result;
end;





test_polytope_integral_equivalence:=function(EXT1, EXT2)
    local TmpDir, File1, File2, FileO, FileE, eProg, TheCommand, result;
    TmpDir:=DirectoryTemporary();
    File1:=Filename(TmpDir, "Test1.in");
    File2:=Filename(TmpDir, "Test2.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(File1, EXT1);
    WriteMatrixFile(File2, EXT2);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Isomorphism");
    TheCommand:=Concatenation(eProg, " mpz_class ", File1, " ", File2, " GAP ", FileO, " 2> ", FileE);
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: The GRP_LinPolytopeIntegral_Isomorphism has not created a file";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFileIfExist(File1);
    RemoveFileIfExist(File2);
    RemoveFileIfExist(FileO);
    RemoveFileIfExist(FileE);
    return result;
end;

get_canonic_form:=function(EXT)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, EXT_can;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Can.in");
    FileO:=Filename(TmpDir, "Can.out");
    FileE:=Filename(TmpDir, "Can.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Canonic");
    TheCommand:=Concatenation(eProg, " mpq_class ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: The GRP_LinPolytope_Canonic did not return anything";
    fi;
    EXT_can:=ReadAsFunction(FileO)();
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    RemoveFileIfExist(FileE);
    return EXT_can;
end;

get_isomorphism_result:=function(EXT1, EXT2)
    local TmpDir, File1, File2, FileO, FileE, eProg, TheCommand, result;
    TmpDir:=DirectoryTemporary();
    File1:=Filename(TmpDir, "Test1.in");
    File2:=Filename(TmpDir, "Test2.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(File1, EXT1);
    WriteMatrixFile(File2, EXT2);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Isomorphism");
    TheCommand:=Concatenation(eProg, " rational ", File1, " ", File2, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: The GRP_LinPolytope_Isomorphism did not return anything";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFile(File1);
    RemoveFile(File2);
    RemoveFile(FileO);
    RemoveFile(FileE);
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
    TheCommand:=Concatenation(eProg, " mpq_class ", FileI, " GAP ", FileO, " 2> ", FileE);
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
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, the_volume;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("POLY_lrs_volume");
    TheCommand:=Concatenation(eProg, " rational ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_lrs_volume should have created a file";
    fi;
    the_volume:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return the_volume;
end;

get_interior_point:=function(FAC)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, TheV;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Interior.in");
    FileO:=Filename(TmpDir, "Interior.out");
    WriteMatrixFile(FileI, FAC);
    eProg:=GetBinaryFilename("POLY_GeometricallyUniqueInteriorPoint");
    TheCommand:=Concatenation(eProg, " rational ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_GeometricallyUniqueInteriorPoint failed to create an output file";
    fi;
    TheV:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
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
    local TmpDir, FileEXT, FileGRP, FileInputNml, FileGrpOut, FileResultOut, FileE, strOut, output, eProg, TheCommand, GRPout, result;
    TmpDir:=DirectoryTemporary();
    FileEXT:=Filename(TmpDir, "Test.ext");
    FileGRP:=Filename(TmpDir, "Test.grp");
    FileInputNml:=Filename(TmpDir, "Test.nml");
    FileGrpOut:=Filename(TmpDir, "Test.grp");
    FileResultOut:=Filename(TmpDir, "Test.result");
    FileE:=Filename(TmpDir, "Test.err");
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
    TheCommand:=Concatenation(eProg, " ", FileInputNml, " 2> ", FileE);
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
    RemoveFile(FileE);
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
    RemoveFile(FileE);
    return EXTint;
end;

get_dual_desc:=function(EXT, method)
    local TmpDir, FileI, FileO, FileE, eProg, command, TheCommand, FAC;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.out");
    FileO:=Filename(TmpDir, "Test.fac");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("POLY_dual_description");
    TheCommand:=Concatenation(eProg, " rational  ", method, " CPP ", FileI, " ", FileO, " 2>", FileE);
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


test_shortest_realizability:=function(SHV)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, eRec;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    #
    WriteMatrixFile(FileI, SHV);
    eProg:=GetBinaryFilename("SHORT_TestRealizability");
    TheCommand:=Concatenation(eProg, " gmp ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: SHORT_TestRealizability did not return a file";
    fi;
    #
    eRec:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return eRec;
end;


quadratic_form_isotropic_test_generic:=function(eProg, mat)
    local TmpDir, FileI, FileO, FileE, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    RemoveFileIfExist(FileE);
    #
    WriteMatrixFile(FileI, mat);
    #
    TheCommand:=Concatenation(eProg, " rational ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: Isotropic program errors";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

quadratic_form_is_isotropic:=function(mat)
    local eProg, result;
    eProg:=GetBinaryFilename("LATT_TestIsotropic");
    result:=quadratic_form_isotropic_test_generic(eProg, mat);
    if is_error(result) then
        return result;
    fi;
    return result.has_isotropic;
end;


quadratic_form_isotropic_vector:=function(mat)
    local eProg, result;
    eProg:=GetBinaryFilename("LATT_FindIsotropic");
    result:=quadratic_form_isotropic_test_generic(eProg, mat);
    if is_error(result) then
        return result;
    fi;
    if result.has_isotropic then
        return result.V;
    else
        return fail;
    fi;
end;

get_lattice_covering:=function(eMat)
    local TmpDir, FileI, FileN, FileO, FileE, output, strOut, eProg, TheCommand, U, UseMpi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "DelaunayEnum.in");
    FileN:=Filename(TmpDir, "DelaunayEnum.nml");
    FileO:=Filename(TmpDir, "DelaunayEnum.out");
    FileE:=Filename(TmpDir, "DelaunayEnum.err");
    #
    WriteMatrixFile(FileI, eMat);
    Print("TestEnumeration, FileIn created\n");
    #
    UseMpi:=false;
    if UseMpi then
        strOut:="&DATA\n";
        strOut:=Concatenation(strOut, " arithmetic_T = \"gmp_rational\"\n");
        strOut:=Concatenation(strOut, " arithmetic_Tint = \"gmp_integer\"\n");
        strOut:=Concatenation(strOut, " GRAMfile = \"", FileI, "\"\n");
        strOut:=Concatenation(strOut, " SVRfile = \"unset.svr\"\n");
        strOut:=Concatenation(strOut, " OutFormat = \"GAP\"\n");
        strOut:=Concatenation(strOut, " OutFile = \"", FileO, "\"\n");
        strOut:=Concatenation(strOut, " max_runtime_second = 10800\n");
        strOut:=Concatenation(strOut, "/\n");
        strOut:=Concatenation(strOut, "\n");
        strOut:=Concatenation(strOut, "&STORAGE\n");
        strOut:=Concatenation(strOut, " Saving = F\n");
        strOut:=Concatenation(strOut, " Prefix = \"DATA/\"\n");
        strOut:=Concatenation(strOut, "/\n");
        #
        WriteStringFile(FileN, strOut);
        Print("TestEnumeration, FileNml created\n");
        #
        eProg:=GetBinaryFilename("LATT_MPI_ComputeDelaunay");
        TheCommand:=Concatenation(eProg, " ", FileN);
    else
        eProg:=GetBinaryFilename("LATT_SerialComputeDelaunay");
        TheCommand:=Concatenation(eProg, " gmp ", FileI, " GAP_Covering ", FileO, " 2> ", FileE);
    fi;
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    #
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return "program failure: get gap covering";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileN);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

get_robust_lattice_covering_density:=function(eMat)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, the_result;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "RobustEnum.in");
    FileO:=Filename(TmpDir, "RobustEnum.out");
    FileE:=Filename(TmpDir, "RobustEnum.err");
    #
    WriteMatrixFile(FileI, eMat);
    Print("TestEnumeration, FileIn created\n");

    eProg:=GetBinaryFilename("Robust_ExactRobustCoveringDensity");
    TheCommand:=Concatenation(eProg, " gmp ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    the_result:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return the_result;
end;


find_positive_vectors:=function(M, CritNorm, StrictIneq)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, V;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    #
    WriteMatrixFile(FileI, M);
    #
    eProg:=GetBinaryFilename("LATT_FindPositiveVector");
    TheCommand:=Concatenation(eProg, " gmp ", FileI, " ", String(CritNorm), " ", String(StrictIneq), " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: LATT_FindPositiveVector did not return anything, likely crash";
    fi;
    V:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return V;
end;


get_orbit_shortest:=function(eG)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, eG);
    #
    eProg:=GetBinaryFilename("LATT_ComputeShortestOrbits");
    TheCommand:=Concatenation(eProg, " gmp ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: LATT_ComputeShortestOrbits did not return anything, likely crash";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
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
    RemoveFile(FileE);
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



test_copositivity:=function(eMat)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, eMat);
    #
    eProg:=GetBinaryFilename("CP_TestCopositivity");
    TheCommand:=Concatenation(eProg, " gmp ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: CP_TestCopositivity did not return anything";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

test_complete_positivity:=function(eMat)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, eMat);
    #
    eProg:=GetBinaryFilename("CP_TestCompletePositivity");
    TheCommand:=Concatenation(eProg, " gmp ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        return "program failure: CP_TestCompletePositivity did not return anything";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
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


test_equivalent_lorentzian_matrices:=function(mat1, mat2)
    local TmpDir, File1, File2, FileO, FileE, eProg, TheCommand, U;
    #
    TmpDir:=DirectoryTemporary();
    File1:=Filename(TmpDir, "Indef_Equi.mat1");
    File2:=Filename(TmpDir, "Indef_Equi.mat2");
    FileO:=Filename(TmpDir, "Indef_Equi.out");
    FileE:=Filename(TmpDir, "Indef_Equi.err");
    WriteMatrixFile(File1, mat1);
    WriteMatrixFile(File2, mat2);
    #
    eProg:=GetBinaryFilename("LORENTZ_PERF_Isomorphism");
    TheCommand:=Concatenation(eProg, " ", File1, " ", File2, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    #
    if IsExistingFile(FileO)=false then
        return "program failure: LORENTZ_PERF_Isomorphism failure for perfect matrices";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(File1);
    RemoveFile(File2);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

stabilizer_lorentzian_matrix:=function(mat)
    local TmpDir, FileI, FileO, FileE, eProg, TheCommand, U;
    #
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Indef_Stab.mat");
    FileO:=Filename(TmpDir, "Indef_Stab.out");
    FileE:=Filename(TmpDir, "Indef_Stab.err");
    #
    WriteMatrixFile(FileI, mat);
    #
    eProg:=GetBinaryFilename("LORENTZ_PERF_Automorphism");
    TheCommand:=Concatenation(eProg, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    #
    if IsExistingFile(FileO)=false then
        return "program failure: LORENTZ_PERF_Automorphism failed to create a file";
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



GenerateTspaceDescription_classic:=function(n, only_well_rounded)
    local CacheFile;
    CacheFile:=Concatenation("/tmp/Perfect_classic_", String(n));
    return rec(nature:="classic",
               n:=n,
               only_well_rounded:=only_well_rounded,
               CacheFile:=CacheFile);
end;

GenerateTspaceDescription_imag_quad:=function(n, d, only_well_rounded)
    local info, CacheFile;
    info:=GetFundamentalInfo(d);
    CacheFile:=Concatenation("/tmp/Perfect_imag_quad_", String(n), "_", String(d));
    return rec(nature:="imag_quad",
               info:=info,
               n:=n,
               only_well_rounded:=only_well_rounded,
               CacheFile:=CacheFile);
end;

FORTRAN_logical:=function(test)
    if test=true then
        return "T";
    fi;
    if test=false then
        return "F";
    fi;
    Error("Wrong input to FORTRAN_logical");
end;


write_t_space:=function(output, desc)
    if desc.nature = "classic" then
        AppendTo(output, "&TSPACE\n");
        AppendTo(output, " TypeTspace = \"Classic\"\n");
        AppendTo(output, " FileLinSpa = \"unset.linspa\"\n");
        AppendTo(output, " SuperMatMethod = \"NotNeeded\"\n");
        AppendTo(output, " ListComm = \"Use_realimag\"\n");
        AppendTo(output, " PtGroupMethod = \"Trivial\"\n");
        AppendTo(output, " FileListSubspaces = \"unset\"\n");
        AppendTo(output, " ClassicDim = ", desc.n, "\n");
        AppendTo(output, "/\n");
        AppendTo(output, "\n");
        return;
    fi;
    if desc.nature = "imag_quad" then
        AppendTo(output, "&TSPACE\n");
        AppendTo(output, " TypeTspace = \"", desc.info.type_tspace, "\"\n");
        AppendTo(output, " FileLinSpa = \"unset.linspa\"\n");
        AppendTo(output, " SuperMatMethod = \"NotNeeded\"\n");
        AppendTo(output, " ListComm = \"Use_realimag\"\n");
#        AppendTo(output, " PtGroupMethod = \"Trivial\"\n");
        AppendTo(output, " PtGroupMethod = \"Compute\"\n");
        AppendTo(output, " FileListSubspaces = \"unset\"\n");
        AppendTo(output, " RealImagDim = ", desc.n, "\n");
        AppendTo(output, " RealImagSum = ", desc.info.eSum, "\n");
        AppendTo(output, " RealImagProd = ", desc.info.eProd, "\n");
        AppendTo(output, "/\n");
        AppendTo(output, "\n");
        return;
    fi;
    Print("nature=", desc.nature, "\n");
    Print("Supported types: classic and imag_quad\n");
    Error("Unsupposed type");
end;



__PERFCOMP_Write_t_space:=function(output, desc)
    AppendTo(output, "&DATA\n");
    AppendTo(output, " arithmetic_T = \"gmp_rational\"\n");
    AppendTo(output, " arithmetic_Tint = \"gmp_integer\"\n");
    AppendTo(output, " OnlyWellRounded = ", FORTRAN_logical(desc.only_well_rounded), "\n");
    AppendTo(output, " ComputeBoundary = T\n");
    AppendTo(output, " ComputeContractingHomotopy = T\n");
#    AppendTo(output, " CacheFile = \"", desc.CacheFile, "\"\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    #
    write_t_space(output, desc);
end;

PERFCOMP_group_generators:=function(desc)
    local TmpDir, FileN, FileO, FileE, output, ListGen, binary, cmd;
    TmpDir:=DirectoryTemporary();
    FileN:=Filename(TmpDir, "PerfComp.nml");
    FileO:=Filename(TmpDir, "PerfComp.out");
    FileE:=Filename(TmpDir, "PerfComp.err");
    output:=OutputTextFile(FileN, true);
    __PERFCOMP_Write_t_space(output, desc);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileGroupGenerators = \"", FileO, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileN, " 2> ", FileE);
    Exec(cmd);
    #
    ListGen:=ReadAsFunction(FileO)();
    RemoveFile(FileN);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListGen;
end;



get_rec_tspace:=function(desc)
    local TmpDir, FileN, FileO, FileE, binary, output, cmd, tspace;
    TmpDir:=DirectoryTemporary();
    FileN:=Filename(TmpDir, "PerfComp.nml");
    FileO:=Filename(TmpDir, "PerfComp.out");
    FileE:=Filename(TmpDir, "PerfComp.err");

    output:=OutputTextFile(FileN, true);
    write_t_space(output, desc);
    CloseStream(output);

    binary:=GetBinaryFilename("TSPACE_FileFormatConversion");
    cmd:=Concatenation(binary, " ", FileN, " GAP ", FileO, " 2> ", FileE);
    Exec(cmd);

    if IsExistingFile(FileO)=false then
        Error("The output file is not existing. That qualifies as a fail");
    fi;

    tspace:=ReadAsFunction(FileO)();
    RemoveFile(FileN);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return tspace;
end;







__PERFCOMP_face_query:=function(desc, ListVect, query)
    local TmpDir, FileN, FileI, FileO, FileE, output, binary, cmd, ListResult;
    TmpDir:=DirectoryTemporary();
    FileN:=Filename(TmpDir, "PerfComp.nml");
    FileI:=Filename(TmpDir, "PerfComp.mat");
    FileO:=Filename(TmpDir, "PerfComp.mat.output");
    FileE:=Filename(TmpDir, "PerfComp.err");
    WriteListMatrixFile(FileI, [ListVect]);
    #
    output:=OutputTextFile(FileN, true);
    __PERFCOMP_Write_t_space(output, desc);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " ", query, " = \"", FileI, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileN, " 2> ", FileE);
    Exec(cmd);
    #
    ListResult:=ReadAsFunction(FileO)();
    RemoveFile(FileN);
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListResult[1];
end;

PERFCOMP_dimension:=function(desc, ListVect)
    return __PERFCOMP_face_query(desc, ListVect, "FileDimensions");
end;



PERFCOMP_is_face:=function(desc, ListVect)
    return __PERFCOMP_face_query(desc, ListVect, "FileIsFace");
end;


PERFCOMP_is_well_rounded:=function(desc, ListVect)
    return __PERFCOMP_face_query(desc, ListVect, "FileIsWellRounded");
end;


PERFCOMP_face_search:=function(desc, ListVect)
    return __PERFCOMP_face_query(desc, ListVect, "FileFaceSearch");
end;

PERFCOMP_stabilizer:=function(desc, ListVect)
    return __PERFCOMP_face_query(desc, ListVect, "FileStabilizerQueries");
end;

PERFCOMP_lower_boundary_cell:=function(desc, ListVect)
    return __PERFCOMP_face_query(desc, ListVect, "FileCellLowerBoundary");
end;

PERFCOMP_upper_boundary_cell:=function(desc, ListVect)
    return __PERFCOMP_face_query(desc, ListVect, "FileCellUpperBoundary");
end;




PERFCOMP_test_equivalence:=function(desc, ListVect1, ListVect2)
    local TmpDir, FileN, FileI, FileO, FileE, output, binary, cmd, ListEquiv;
    TmpDir:=DirectoryTemporary();
    FileN:=Filename(TmpDir, "PerfComp.nml");
    FileI:=Filename(TmpDir, "PerfComp.mat");
    FileO:=Filename(TmpDir, "PerfComp.mat.output");
    FileE:=Filename(TmpDir, "PerfComp.err");
    WriteListMatrixFile(FileI, [ListVect1, ListVect2]);
    #
    output:=OutputTextFile(FileN, true);
    __PERFCOMP_Write_t_space(output, desc);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileEquivalenceQueries = \"", FileI, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileN, " 2> ", FileE);
    Exec(cmd);
    #
    ListEquiv:=ReadAsFunction(FileO)();
    RemoveFile(FileN);
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListEquiv[1];
end;



PERFCOMP_get_cells:=function(desc, index)
    local TmpDir, FileN, FileI, FileO, FileE, output, binary, cmd, ListCells;
    TmpDir:=DirectoryTemporary();
    FileN:=Filename(TmpDir, "PerfComp.nml");
    FileO:=Filename(TmpDir, "PerfComp.out");
    FileE:=Filename(TmpDir, "PerfComp.err");
    #
    output:=OutputTextFile(FileN, true);
    __PERFCOMP_Write_t_space(output, desc);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileCells = \"", FileO, "\"\n");
    AppendTo(output, " IndexCell = ", index, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileN, " 2> ", FileE);
    Exec(cmd);
    #
    ListCells:=ReadAsFunction(FileO)();
    RemoveFile(FileN);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListCells;
end;



PERFCOMP_get_lower_cells:=function(desc, index)
    local TmpDir, FileN, FileI, FileO, FileE, output, binary, cmd, ListLower;
    TmpDir:=DirectoryTemporary();
    FileN:=Filename(TmpDir, "PerfComp.nml");
    FileO:=Filename(TmpDir, "PerfComp.out");
    FileE:=Filename(TmpDir, "PerfComp.err");
    #
    output:=OutputTextFile(FileN, true);
    __PERFCOMP_Write_t_space(output, desc);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileListLowerBoundary = \"", FileO, "\"\n");
    AppendTo(output, " IndexLowerBoundary = ", index, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileN, " 2> ", FileE);
    Exec(cmd);
    #
    ListLower:=ReadAsFunction(FileO)();
    RemoveFile(FileN);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListLower;
end;


PERFCOMP_get_upper_cells:=function(desc, index)
    local TmpDir, FileN, FileI, FileO, FileE, output, binary, cmd, ListUpper;
    TmpDir:=DirectoryTemporary();
    FileN:=Filename(TmpDir, "PerfComp.nml");
    FileO:=Filename(TmpDir, "PerfComp.out");
    FileE:=Filename(TmpDir, "PerfComp.err");
    #
    output:=OutputTextFile(FileN, true);
    __PERFCOMP_Write_t_space(output, desc);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " FileListUpperBoundary = \"", FileO, "\"\n");
    AppendTo(output, " IndexUpperBoundary = ", index, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileN, " 2> ", FileE);
    Exec(cmd);
    #
    ListUpper:=ReadAsFunction(FileO)();
    RemoveFile(FileN);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListUpper;
end;















WriteChainStream:=function(output, the_chain)
    local len_chain, entry;
    len_chain:=Length(the_chain);
    AppendTo(output, len_chain, "\n");
    for entry in the_chain
    do
        AppendTo(output, entry.iOrb - 1, "\n");
        AppendTo(output, entry.value, "\n");
        WriteMatrix(output, entry.M);
    od;
end;


__PERFCOMP_chain_query:=function(desc, index, the_chain, query)
    local TmpDir, FileN, FileI, FileO, FileE, output, binary, cmd, ListDim;
    TmpDir:=DirectoryTemporary();
    FileN:=Filename(TmpDir, "PerfComp.nml");
    FileI:=Filename(TmpDir, "PerfComp.chain");
    FileO:=Filename(TmpDir, "PerfComp.chain.output");
    FileE:=Filename(TmpDir, "PerfComp.err");
    output:=OutputTextFile(FileI, true);
    AppendTo(output, index, "\n");
    AppendTo(output, "1\n");
    WriteChainStream(output, the_chain);
    CloseStream(output);
    #
    output:=OutputTextFile(FileN, true);
    __PERFCOMP_Write_t_space(output, desc);
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " ", query, " = \"", FileI, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
#    cmd:=Concatenation(binary, " ", FileN, " 2> ", FileE);
    cmd:=Concatenation(binary, " ", FileN);
    Exec(cmd);
    #
    ListDim:=ReadAsFunction(FileO)();
    RemoveFile(FileN);
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListDim[1];
end;






PERFCOMP_chain_contracting_homotopy:=function(desc, index, the_chain)
    return __PERFCOMP_chain_query(desc, index, the_chain, "FileChainContractingHomotopy");
end;

PERFCOMP_chain_boundary:=function(desc, index, the_chain)
    return __PERFCOMP_chain_query(desc, index, the_chain, "FileChainBoundary");
end;

PERFCOMP_chain_simplification:=function(desc, index, the_chain)
    return __PERFCOMP_chain_query(desc, index, the_chain, "FileChainSimplification");
end;
