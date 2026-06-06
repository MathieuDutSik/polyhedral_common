is_error:=function(input)
    local test;
    if IsString(input)=false then
        return false;
    fi;
    test:=starts_with(input, "program failure")<>fail;
    if test then
        Print("Something went wrong with error=", input, "\n");
    fi;
    return test;
end;

extract_runtime_from_log:=function(FileE)
    local lines, line, prefix, rest, idx;
    if IsExistingFile(FileE)=false then
        return "unknown";
    fi;
    lines:=SplitString(StringFile(FileE), "\n");
    prefix:="runtime = ";
    for line in lines do
        if Length(line) >= Length(prefix) and line{[1..Length(prefix)]}=prefix then
            rest:=line{[Length(prefix)+1..Length(line)]};
            idx:=PositionSublist(rest, " timeanddate=");
            if idx <> fail then
                rest:=rest{[1..idx-1]};
            fi;
            return rest;
        fi;
    od;
    return "unknown";
end;

get_grp_automorphy:=function(arg)
    local EXT, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, TheGRP_rat, runtime_str;
    EXT:=arg[1];
    arith:="rational";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Automorphism");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " RecGAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " command=GRP_LinPolytope_Automorphism runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytope_Automorphism did not return anything, likely crash";
    fi;
    TheGRP_rat:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return TheGRP_rat;
end;

get_grp_integral_automorphy:=function(arg)
    local EXT, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, TheGRP_int, runtime_str;
    EXT:=arg[1];
    arith:="mpz_class";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Automorphism");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " RecGAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " command=GRP_LinPolytopeIntegral_Automorphism runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytopeIntegral_Automorphism did not return anything, likely crash";
    fi;
    TheGRP_int:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return TheGRP_int;
end;

get_double_cosets_matrix_group:=function(arg)
    local EXT, GRPperm, options, arith, print_info, TmpDir, FileI, FileO, FileE, FileGRP_V, eProg, TheCommand, RecResult, runtime_str;
    EXT:=arg[1];
    GRPperm:=arg[2];
    arith:="mpz_class";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    FileGRP_V:=Filename(TmpDir, "Test.grp_V");
    WriteMatrixFile(FileI, EXT);
    WriteGroupFile(FileGRP_V, Length(EXT), GRPperm);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Automorphism_DoubleCoset");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " ", FileGRP_V, " RecGAP ", FileO, " 2> ", FileE);
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " command=GRP_LinPolytopeIntegral_Automorphism_DoubleCoset runtime=", runtime_str, "\n");
    fi;
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

get_linear_space_stabilizer_double_cosets:=function(arg)
    local GRPmatr, TheSpace, GRP_Vmatr, options, print_info, TmpDir, FileGRP, FileSPA, FileGRP_V, FileO, FileE, eProg, TheCommand, RecResult, runtime_str;
    GRPmatr:=arg[1];
    TheSpace:=arg[2];
    GRP_Vmatr:=arg[3];
    print_info:=false;
    if Length(arg) >= 4 then
        options:=arg[4];
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
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
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  TheSpace=", Length(TheSpace), "x", Length(TheSpace[1]), " command=GRP_LinearSpace_Stabilizer_DoubleCoset runtime=", runtime_str, "\n");
    fi;
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


get_linpolytopeintegral_aut_rightcoset:=function(arg)
    local EXT, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, result, runtime_str;
    EXT:=arg[1];
    arith:="mpz_class";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Automorphism_RightCoset");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " command=GRP_LinPolytopeIntegral_Automorphism_RightCoset runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_LinPolytopeIntegral_Automorphism_RightCoset should have created a file";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return result;
end;





test_polytope_integral_equivalence:=function(arg)
    local EXT1, EXT2, options, arith, print_info, TmpDir, File1, File2, FileO, FileE, eProg, TheCommand, result, runtime_str;
    EXT1:=arg[1];
    EXT2:=arg[2];
    arith:="mpz_class";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    File1:=Filename(TmpDir, "Test1.in");
    File2:=Filename(TmpDir, "Test2.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(File1, EXT1);
    WriteMatrixFile(File2, EXT2);
    eProg:=GetBinaryFilename("GRP_LinPolytopeIntegral_Isomorphism");
    TheCommand:=Concatenation(eProg, " ", arith, " ", File1, " ", File2, " GAP ", FileO, " 2> ", FileE);
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT1=", Length(EXT1), "x", Length(EXT1[1]), " EXT2=", Length(EXT2), "x", Length(EXT2[1]), " arith=", arith, " command=GRP_LinPolytopeIntegral_Isomorphism runtime=", runtime_str, "\n");
    fi;
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

get_canonic_form:=function(arg)
    local EXT, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, EXT_can, runtime_str;
    EXT:=arg[1];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Can.in");
    FileO:=Filename(TmpDir, "Can.out");
    FileE:=Filename(TmpDir, "Can.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Canonic");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " command=GRP_LinPolytope_Canonic runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: The GRP_LinPolytope_Canonic did not return anything";
    fi;
    EXT_can:=ReadAsFunction(FileO)();
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    RemoveFileIfExist(FileE);
    return EXT_can;
end;

get_isomorphism_result:=function(arg)
    local EXT1, EXT2, options, arith, print_info, TmpDir, File1, File2, FileO, FileE, eProg, TheCommand, result, runtime_str;
    EXT1:=arg[1];
    EXT2:=arg[2];
    arith:="rational";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    File1:=Filename(TmpDir, "Test1.in");
    File2:=Filename(TmpDir, "Test2.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(File1, EXT1);
    WriteMatrixFile(File2, EXT2);
    eProg:=GetBinaryFilename("GRP_LinPolytope_Isomorphism");
    TheCommand:=Concatenation(eProg, " ", arith, " ", File1, " ", File2, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT1=", Length(EXT1), "x", Length(EXT1[1]), " EXT2=", Length(EXT2), "x", Length(EXT2[1]), " arith=", arith, " command=GRP_LinPolytope_Isomorphism runtime=", runtime_str, "\n");
    fi;
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


get_saturated_space:=function(arg)
    local ListMat, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, ListMatRet, runtime_str;
    ListMat:=arg[1];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Saturation.in");
    FileO:=Filename(TmpDir, "Saturation.out");
    FileE:=Filename(TmpDir, "Saturation.err");
    WriteListMatrixFile(FileI, ListMat);

    eProg:=GetBinaryFilename("TSPACE_IntegralSaturation");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  |ListMat|=", Length(ListMat), " arith=", arith, " command=TSPACE_IntegralSaturation runtime=", runtime_str, "\n");
    fi;

    ListMatRet:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListMatRet;
end;

get_latt_canonical_form:=function(arg)
    local eMat, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, U, runtime_str;
    eMat:=arg[1];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Can.in");
    FileO:=Filename(TmpDir, "Can.out");
    FileE:=Filename(TmpDir, "Can.err");
    WriteMatrixFile(FileI, eMat);
    #
    eProg:=GetBinaryFilename("LATT_Canonicalize");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP_full ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eMat=", Length(eMat), "x", Length(eMat[1]), " arith=", arith, " command=LATT_Canonicalize runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: The LATT_Canonicalize has failed";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

get_latt_automorphism_group:=function(arg)
    local eMat, options, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, U, runtime_str;
    eMat:=arg[1];
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Aut.in");
    FileO:=Filename(TmpDir, "Aut.out");
    FileE:=Filename(TmpDir, "Aut.err");
    WriteListMatrixFile(FileI, [eMat]);
    #
    eProg:=GetBinaryFilename("LATT_Automorphism");
    TheCommand:=Concatenation(eProg, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eMat=", Length(eMat), "x", Length(eMat[1]), " command=LATT_Automorphism runtime=", runtime_str, "\n");
    fi;
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

get_fullrank_invariant_family:=function(arg)
    local eG, method, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, U, runtime_str;
    eG:=arg[1];
    method:=arg[2];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Fullrank.in");
    FileO:=Filename(TmpDir, "Fullrank.out");
    FileE:=Filename(TmpDir, "Fullrank.err");
    WriteMatrixFile(FileI, eG);
    eProg:=GetBinaryFilename("LATT_GenerateCharacteristicVectorSet");
    TheCommand:=Concatenation(eProg, " ", arith, " ", method, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eG=", Length(eG), "x", Length(eG[1]), " arith=", arith, " method=", method, " command=LATT_GenerateCharacteristicVectorSet runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: LATT_GenerateCharacteristicVectorSet did not return anything";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

get_ext_volume:=function(arg)
    local EXT, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, the_volume, runtime_str;
    EXT:=arg[1];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("POLY_lrs_volume");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " command=POLY_lrs_volume runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_lrs_volume should have created a file";
    fi;
    the_volume:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return the_volume;
end;

get_triangulation_of_polytope:=function(arg)
    local EXT, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, one_triangulation, runtime_str;
    EXT:=arg[1];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("POLY_lrs_triangulation");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " command=POLY_lrs_triangulation runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_lrs_volume should have created a file";
    fi;
    one_triangulation:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return one_triangulation;
end;


get_polygen_difference:=function(arg)
    local gp1, gp2, options, arith, print_info, TmpDir, File1, File2, FileO, FileE, eProg, TheCommand, result, runtime_str;
    gp1:=arg[1];
    gp2:=arg[2];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    File1:=Filename(TmpDir, "Test.in1");
    File2:=Filename(TmpDir, "Test.in2");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteListMatrixFile(File1, gp1);
    WriteListMatrixFile(File2, gp2);
    eProg:=GetBinaryFilename("PolyGen_Difference");
    TheCommand:=Concatenation(eProg, " ", arith, " ", File1, " ", File2, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  |gp1|=", Length(gp1), " |gp2|=", Length(gp2), " arith=", arith, " command=PolyGen_Difference runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: PolyGen_Difference should have created a file";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFile(File1);
    RemoveFile(File2);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return result;
end;

get_polygen_vertices:=function(arg)
    local gp, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, result, runtime_str;
    gp:=arg[1];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteListMatrixFile(FileI, gp);
    eProg:=GetBinaryFilename("PolyGen_vertices");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  |gp|=", Length(gp), " arith=", arith, " command=PolyGen_vertices runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: PolyGen_vertices should have created a file";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return result;
end;

get_interior_point:=function(arg)
    local FAC, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, TheV, runtime_str;
    FAC:=arg[1];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Interior.in");
    FileO:=Filename(TmpDir, "Interior.out");
    WriteMatrixFile(FileI, FAC);
    eProg:=GetBinaryFilename("POLY_GeometricallyUniqueInteriorPoint");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  FAC=", Length(FAC), "x", Length(FAC[1]), " arith=", arith, " command=POLY_GeometricallyUniqueInteriorPoint runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_GeometricallyUniqueInteriorPoint failed to create an output file";
    fi;
    TheV:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return TheV;
end;




get_latt_isomorphism_test:=function(arg)
    local eMat1, eMat2, options, print_info, TmpDir, FileIn1, FileIn2, FileOut, FileErr, eProg, TheCommand, U, runtime_str;
    eMat1:=arg[1];
    eMat2:=arg[2];
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
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
    if print_info then
        runtime_str:=extract_runtime_from_log(FileErr);
        Print("  eMat1=", Length(eMat1), "x", Length(eMat1[1]), " eMat2=", Length(eMat2), "x", Length(eMat2[1]), " command=LATT_Isomorphism runtime=", runtime_str, "\n");
    fi;
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



get_group_skeleton:=function(arg)
    local EXT, GRP, LevSearch, options, arith, print_info, TmpDir, FileEXT, FileGRP, FileInputNml, FileGrpOut, FileResultOut, FileE, strOut, output, eProg, TheCommand, GRPout, result, runtime_str;
    EXT:=arg[1];
    GRP:=arg[2];
    LevSearch:=arg[3];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 4 then
        options:=arg[4];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
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
    strOut:=Concatenation(strOut, " Arithmetic = \"", arith, "\"\n");
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
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " LevSearch=", LevSearch, " command=POLY_FaceLatticeGen runtime=", runtime_str, "\n");
    fi;
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


sample_facet_polytope:=function(arg)
    local EXT, method, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, result, runtime_str;
    EXT:=arg[1];
    method:=arg[2];
    arith:="rational";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    #
    eProg:=GetBinaryFilename("POLY_sampling_facets");
    TheCommand:=Concatenation(eProg, " ", arith, " ", method, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " method=", method, " command=POLY_sampling_facets runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_sampling_facets failed to create the output";
    fi;
    result:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return result;
end;


test_matrix_group_finiteness:=function(arg)
    local GRPmatr, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, U, runtime_str;
    GRPmatr:=arg[1];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Finiteness_test.input");
    FileO:=Filename(TmpDir, "Finiteness_test.result");
    FileE:=Filename(TmpDir, "Finiteness_test.err");
    WriteListMatrixFile(FileI, GRPmatr);
    #
    eProg:=GetBinaryFilename("GRP_TestFiniteness");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    #
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  |GRPmatr|=", Length(GRPmatr), " arith=", arith, " command=GRP_TestFiniteness runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: GRP_TestFiniteness failed to create the output";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;



get_integral_interior_point:=function(arg)
    local FAC, method, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, EXTint, runtime_str;
    FAC:=arg[1];
    method:=arg[2];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.fac");
    FileO:=Filename(TmpDir, "Test.vertint");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, FAC);
    eProg:=GetBinaryFilename("POLY_IntegralPoints");
    TheCommand:=Concatenation(eProg, " ", arith, " ", method, " ", FileI, " CPP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  FAC=", Length(FAC), "x", Length(FAC[1]), " arith=", arith, " method=", method, " command=POLY_IntegralPoints runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: POLY_IntegralPoints failed to create the file";
    fi;
    EXTint:=ReadMatrixFile(FileO);
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return EXTint;
end;

get_dual_desc_kernel:=function(arg)
    local choice, EXT, method, options, arith, print_info,
          TmpDir, FileI, FileO, FileE, eProg, TheCommand, FAC, runtime_str;
    choice:=arg[1];
    EXT:=arg[2];
    method:=arg[3];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 4 then
        options:=arg[4];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("POLY_dual_description");
    TheCommand:=Concatenation(eProg, " ", arith, " ", method, " ", choice, " ", FileI, " ", FileO, " 2>", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " command=", method, " runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        RemoveFile(FileI);
        RemoveFile(FileE);
        return "program failure: POLY_dual_description failed to create the file";
    fi;
    FAC:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return FAC;
end;

get_dual_desc:=function(arg)
    return CallFuncList(get_dual_desc_kernel, Concatenation(["GAP"], arg));
end;

get_dual_desc_incidence:=function(arg)
    return CallFuncList(get_dual_desc_kernel, Concatenation(["GAPincidence"], arg));
end;


# Run POLY_cdd_skeletons (cdd::DualDescriptionAdjacencies).
# Returns rec(FAC, nbVertSkel, SkelEdges, nbVertRidge, RidgeEdges) on
# success, or a string starting with "program failure: ..." on error.
get_cdd_skeletons:=function(arg)
    local EXT, options, arith, print_info,
          TmpDir, FileI, FileO, FileE, eProg, TheCommand, eRec, runtime_str;
    EXT:=arg[1];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, EXT);
    eProg:=GetBinaryFilename("POLY_cdd_skeletons");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO,
                              " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  EXT=", Length(EXT), "x", Length(EXT[1]), " arith=", arith, " command=POLY_cdd_skeletons runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        RemoveFile(FileI);
        RemoveFile(FileE);
        return "program failure: POLY_cdd_skeletons did not create the output";
    fi;
    eRec:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return eRec;
end;


# Run POLY_LinearDetermineByInequalities.  Given an H-representation
# matrix FAC, returns a matrix whose rows span the affine span of the
# polytope (so a 4x4 identity when the polytope is full-dimensional in
# 4 columns, fewer rows when the system has implicit linearities).
# Returns a string starting with "program failure: ..." on error.
get_linear_determine_by_inequalities:=function(arg)
    local FAC, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, Mat, runtime_str;
    FAC:=arg[1];
    arith:="mpq_class";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, FAC);
    eProg:=GetBinaryFilename("POLY_LinearDetermineByInequalities");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO,
                              " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  FAC=", Length(FAC), "x", Length(FAC[1]), " arith=", arith, " command=POLY_LinearDetermineByInequalities runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        RemoveFile(FileI);
        RemoveFile(FileE);
        return Concatenation("program failure: ",
                             "POLY_LinearDetermineByInequalities ",
                             "did not create the output");
    fi;
    Mat:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return Mat;
end;


test_shortest_realizability:=function(arg)
    local SHV, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, eRec, runtime_str;
    SHV:=arg[1];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    #
    WriteMatrixFile(FileI, SHV);
    eProg:=GetBinaryFilename("SHORT_TestRealizability");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  SHV=", Length(SHV), "x", Length(SHV[1]), " arith=", arith, " command=SHORT_TestRealizability runtime=", runtime_str, "\n");
    fi;
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


quadratic_form_isotropic_test_generic:=function(arg)
    local eProg, mat, options, arith, print_info, TmpDir, FileI, FileO, FileE, TheCommand, U, runtime_str;
    eProg:=arg[1];
    mat:=arg[2];
    arith:="rational";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
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
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  mat=", Length(mat), "x", Length(mat[1]), " arith=", arith, " command=", eProg, " runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: Isotropic program errors";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

quadratic_form_is_isotropic:=function(arg)
    local mat, eProg, result;
    mat:=arg[1];
    eProg:=GetBinaryFilename("LATT_TestIsotropic");
    result:=CallFuncList(quadratic_form_isotropic_test_generic, Concatenation([eProg], arg));
    if is_error(result) then
        return result;
    fi;
    return result.has_isotropic;
end;


quadratic_form_isotropic_vector:=function(arg)
    local mat, eProg, result;
    mat:=arg[1];
    eProg:=GetBinaryFilename("LATT_FindIsotropic");
    result:=CallFuncList(quadratic_form_isotropic_test_generic, Concatenation([eProg], arg));
    if is_error(result) then
        return result;
    fi;
    if result.has_isotropic then
        return result.V;
    else
        return fail;
    fi;
end;

get_lattice_covering:=function(arg)
    local eMat, options, arith, print_info, TmpDir, FileI, FileN, FileO, FileE, output, strOut, eProg, TheCommand, U, UseMpi, runtime_str;
    eMat:=arg[1];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
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
        strOut:=Concatenation(strOut, " arithmetic = \"gmp\"\n");
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
        TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP_Covering ", FileO, " 2> ", FileE);
    fi;
#    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eMat=", Length(eMat), "x", Length(eMat[1]), " arith=", arith, " command=LATT_SerialComputeDelaunay runtime=", runtime_str, "\n");
    fi;
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

get_robust_lattice_covering_density:=function(arg)
    local eMat, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, the_result, runtime_str;
    eMat:=arg[1];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "RobustEnum.in");
    FileO:=Filename(TmpDir, "RobustEnum.out");
    FileE:=Filename(TmpDir, "RobustEnum.err");
    #
    WriteMatrixFile(FileI, eMat);
    Print("TestEnumeration, FileIn created\n");

    eProg:=GetBinaryFilename("Robust_ExactRobustCoveringDensity");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eMat=", Length(eMat), "x", Length(eMat[1]), " arith=", arith, " command=Robust_ExactRobustCoveringDensity runtime=", runtime_str, "\n");
    fi;
    the_result:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return the_result;
end;


find_positive_vectors:=function(arg)
    local M, CritNorm, StrictIneq, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, V, runtime_str;
    M:=arg[1];
    CritNorm:=arg[2];
    StrictIneq:=arg[3];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 4 then
        options:=arg[4];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    #
    WriteMatrixFile(FileI, M);
    #
    eProg:=GetBinaryFilename("LATT_FindPositiveVector");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " ", String(CritNorm), " ", String(StrictIneq), " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  M=", Length(M), "x", Length(M[1]), " arith=", arith, " CritNorm=", CritNorm, " command=LATT_FindPositiveVector runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: LATT_FindPositiveVector did not return anything, likely crash";
    fi;
    V:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return V;
end;


get_orbit_shortest:=function(arg)
    local eG, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, U, runtime_str;
    eG:=arg[1];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, eG);
    #
    eProg:=GetBinaryFilename("LATT_ComputeShortestOrbits");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eG=", Length(eG), "x", Length(eG[1]), " arith=", arith, " command=LATT_ComputeShortestOrbits runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: LATT_ComputeShortestOrbits did not return anything, likely crash";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;







INDEF_FORM_Stabilizer:=function(arg)
    local eGram, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, U, TheCommand, runtime_str;
    eGram:=arg[1];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=GetFreeFile("Err_indef_form_autom_");
    WriteMatrixFile(FileI, eGram);
    eProg:=GetBinaryFilename("INDEF_FORM_AutomorphismGroup");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eGram=", Length(eGram), "x", Length(eGram[1]), " arith=", arith, " command=INDEF_FORM_AutomorphismGroup runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: INDEF_FORM_Stabilizer did not return anything, likely crash";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return Group(U);
end;

INDEF_FORM_TestEquivalence:=function(arg)
    local eGram1, eGram2, options, arith, print_info, TmpDir, eProg, FileM1, FileM2, FileO, FileE, TheCommand, U, runtime_str;
    eGram1:=arg[1];
    eGram2:=arg[2];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileM1:=Filename(TmpDir, "Mat1.in");
    FileM2:=Filename(TmpDir, "Mat2.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileM1, eGram1);
    WriteMatrixFile(FileM2, eGram2);
    eProg:=GetBinaryFilename("INDEF_FORM_TestEquivalence");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileM1, " ", FileM2, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eGram1=", Length(eGram1), "x", Length(eGram1[1]), " eGram2=", Length(eGram2), "x", Length(eGram2[1]), " arith=", arith, " command=INDEF_FORM_TestEquivalence runtime=", runtime_str, "\n");
    fi;
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


INDEF_FORM_GetOrbitRepresentative:=function(arg)
    local eGram, Xnorm, options, arith, print_info, TmpDir, eProg, FileM, FileO, FileE, TheCommand, U, runtime_str;
    eGram:=arg[1];
    Xnorm:=arg[2];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileM:=Filename(TmpDir, "Mat.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileM, eGram);
    eProg:=GetBinaryFilename("INDEF_FORM_GetOrbitRepresentative");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileM, " ", String(Xnorm), " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eGram=", Length(eGram), "x", Length(eGram[1]), " arith=", arith, " Xnorm=", Xnorm, " command=INDEF_FORM_GetOrbitRepresentative runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: Nothing returned by INDEF_FORM_GetOrbitRepresentative, likely crash or something";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileM);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

INDEF_FORM_GetOrbitIsotropic:=function(arg)
    local eGram, kDim, choice, options, arith, print_info, TmpDir, eProg, FileM, FileO, FileE, TheCommand, U, runtime_str;
    eGram:=arg[1];
    kDim:=arg[2];
    choice:=arg[3];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 4 then
        options:=arg[4];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileM:=Filename(TmpDir, "Mat.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileM, eGram);
    eProg:=GetBinaryFilename("INDEF_FORM_GetOrbit_IsotropicKplane");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileM, " ", String(kDim), " ", choice, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eGram=", Length(eGram), "x", Length(eGram[1]), " arith=", arith, " kDim=", kDim, " choice=", choice, " command=INDEF_FORM_GetOrbit_IsotropicKplane runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: Nothing returned by INDEF_FORM_GetOrbitIsotropic, likely crash or something";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileM);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;





GetReflectivityInformation:=function(arg)
    local LorMat, options, print_info, n, TmpDir, FileI, FileN, FileO, FileE, strOut, eProg, TheCommand, U, runtime_str;
    LorMat:=arg[1];
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
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
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  LorMat=", Length(LorMat), "x", Length(LorMat[1]), " command=LORENTZ_FundDomain_AllcockEdgewalk runtime=", runtime_str, "\n");
    fi;
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



test_copositivity:=function(arg)
    local eMat, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, U, runtime_str;
    eMat:=arg[1];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, eMat);
    #
    eProg:=GetBinaryFilename("CP_TestCopositivity");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eMat=", Length(eMat), "x", Length(eMat[1]), " arith=", arith, " command=CP_TestCopositivity runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: CP_TestCopositivity did not return anything";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;

test_complete_positivity:=function(arg)
    local eMat, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, U, runtime_str;
    eMat:=arg[1];
    arith:="gmp";
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteMatrixFile(FileI, eMat);
    #
    eProg:=GetBinaryFilename("CP_TestCompletePositivity");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  eMat=", Length(eMat), "x", Length(eMat[1]), " arith=", arith, " command=CP_TestCompletePositivity runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: CP_TestCompletePositivity did not return anything";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;



get_matrix_group_mod_information:=function(arg)
    local ListGen, p_val, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, U, runtime_str;
    ListGen:=arg[1];
    p_val:=arg[2];
    arith:="mpz_class";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteListMatrixFile(FileI, ListGen);
    eProg:=GetBinaryFilename("LATT_ResolveModAction");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " ", String(p_val), " GAP ", FileO);
#    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " ", String(p_val), " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  |ListGen|=", Length(ListGen), " arith=", arith, " p_val=", p_val, " command=LATT_ResolveModAction runtime=", runtime_str, "\n");
    fi;
    if IsExistingFile(FileO)=false then
        return "program failure: Nothing returned by LATT_ResolveModAction, likely crash or something";
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileI);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return U;
end;


test_equivalent_lorentzian_matrices:=function(arg)
    local mat1, mat2, options, print_info, TmpDir, File1, File2, FileO, FileE, eProg, TheCommand, U, runtime_str;
    mat1:=arg[1];
    mat2:=arg[2];
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
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
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  mat1=", Length(mat1), "x", Length(mat1[1]), " mat2=", Length(mat2), "x", Length(mat2[1]), " command=LORENTZ_PERF_Isomorphism runtime=", runtime_str, "\n");
    fi;
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

stabilizer_lorentzian_matrix:=function(arg)
    local mat, options, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, U, runtime_str;
    mat:=arg[1];
    print_info:=false;
    if Length(arg) >= 2 then
        options:=arg[2];
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
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
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  mat=", Length(mat), "x", Length(mat[1]), " command=LORENTZ_PERF_Automorphism runtime=", runtime_str, "\n");
    fi;
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



get_grp_size_matrix_group_mod_action:=function(arg)
    local ListGen, p_val, options, arith, print_info, TmpDir, FileI, FileO, FileE, eProg, TheCommand, U, runtime_str;
    ListGen:=arg[1];
    p_val:=arg[2];
    arith:="mpz_class";
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arith) then
            arith:=options.arith;
        fi;
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.in");
    FileO:=Filename(TmpDir, "Test.out");
    FileE:=Filename(TmpDir, "Test.err");
    WriteListMatrixFile(FileI, ListGen);
    eProg:=GetBinaryFilename("LATT_ComputeGroupModAction");
    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " ", String(p_val), " GAP ", FileO);
#    TheCommand:=Concatenation(eProg, " ", arith, " ", FileI, " ", String(p_val), " GAP ", FileO, " 2> ", FileE);
    Exec(TheCommand);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  |ListGen|=", Length(ListGen), " arith=", arith, " p_val=", p_val, " command=LATT_ComputeGroupModAction runtime=", runtime_str, "\n");
    fi;
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



__PERFCOMP_Write_t_space:=function(arg)
    local output, desc, options, arithmetic;
    output:=arg[1];
    desc:=arg[2];
    arithmetic:="gmp";
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.arithmetic) then
            arithmetic:=options.arithmetic;
        fi;
    fi;
    AppendTo(output, "&DATA\n");
    AppendTo(output, " arithmetic = \"", arithmetic, "\"\n");
    AppendTo(output, " OnlyWellRounded = ", FORTRAN_logical(desc.only_well_rounded), "\n");
    AppendTo(output, " ComputeBoundary = T\n");
    AppendTo(output, " ComputeContractingHomotopy = T\n");
#    AppendTo(output, " CacheFile = \"", desc.CacheFile, "\"\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    #
    write_t_space(output, desc);
end;

PERFCOMP_general_query:=function(arg)
    local desc, general_query, options, print_info, TmpDir, FileN, FileO, FileE, output, ListGen, binary, cmd, runtime_str;
    desc:=arg[1];
    general_query:=arg[2];
    print_info:=false;
    if Length(arg) >= 3 then
        options:=arg[3];
        if IsBound(options.print_info) and options.print_info then
            print_info:=true;
        fi;
    fi;
    TmpDir:=DirectoryTemporary();
    FileN:=Filename(TmpDir, "PerfComp.nml");
    FileO:=Filename(TmpDir, "PerfComp.out");
    FileE:=Filename(TmpDir, "PerfComp.err");
    output:=OutputTextFile(FileN, true);
    CallFuncList(__PERFCOMP_Write_t_space, Concatenation([output, desc], arg{[3..Length(arg)]}));
    AppendTo(output, "&QUERIES\n");
    AppendTo(output, " ", general_query, " = \"", FileO, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    binary:=GetBinaryFilename("PERF_SerialPerfectComputation");
    cmd:=Concatenation(binary, " ", FileN, " 2> ", FileE);
    Exec(cmd);
    if print_info then
        runtime_str:=extract_runtime_from_log(FileE);
        Print("  desc.nature=", desc.nature, " command=PERF_SerialPerfectComputation/FileGroupGenerators runtime=", runtime_str, "\n");
    fi;
    #
    ListGen:=ReadAsFunction(FileO)();
    RemoveFile(FileN);
    RemoveFile(FileO);
    RemoveFile(FileE);
    return ListGen;
end;


PERFCOMP_group_generators:=function(arg)
    local desc, general_query;
    desc:=arg[1];
    general_query:="FileGroupGenerators";
    return CallFuncList(PERFCOMP_general_query, Concatenation([desc, general_query], arg{[2..Length(arg)]}));
end;


PERFCOMP_list_number_orbit:=function(arg)
    local desc, general_query;
    desc:=arg[1];
    general_query:="FileListNumberOrbit";
    return CallFuncList(PERFCOMP_general_query, Concatenation([desc, general_query], arg{[2..Length(arg)]}));
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
