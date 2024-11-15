Read("../common.g");
Print("Beginning TestReflectivity\n");

GeneratePoincareFundamentalInput:=false;
GenerateIsotropicTestCases:=true;
iPoincare:=0;

LORENTZ_GetReflection:=function(LorMat, eRoot)
    local n, eNorm, TheMat, eLine, eScal, fLine;
    n:=Length(eRoot);
    eNorm:=eRoot * LorMat * eRoot;
    TheMat:=[];
    for eLine in IdentityMat(n)
    do
        eScal:=eLine * LorMat * eRoot;
        fLine:=eLine - 2 * (eScal / eNorm) * eRoot;
        Add(TheMat, fLine);
    od;
    if TheMat * LorMat * TransposedMat(TheMat) <> LorMat then
        Error("TheMat does not preserve the quadratic form");
    fi;
    if eRoot * TheMat<>-eRoot then
        Error("TheMat does not reflect eRoot");
    fi;
    if IsIntegralMat(TheMat)=false then
        Error("TheMat should be integral");
    fi;
    if TheMat * TheMat<>IdentityMat(n) then
        Error("The matrix TheMat *TheMat is not the identity");
    fi;
    return TheMat;
end;

AllAllowed:=true;

WritePoincareCase:=function(PrefixPoincare, ThePt, ListGen)
    local FilePoincare_Data, FilePoincare_Nml, output, strOut, eGen;
    FilePoincare_Data:=Concatenation(PrefixPoincare, ".data");
    FilePoincare_Nml:=Concatenation(PrefixPoincare, ".nml");
    output:=OutputTextFile(FilePoincare_Data, true);
    WriteVector(output, ThePt);
    AppendTo(output, Length(ListGen), "\n");
    for eGen in ListGen
    do
        WriteMatrix(output, eGen);
    od;
    CloseStream(output);
    #
    strOut:="&PROC\n";
    strOut:=Concatenation(strOut, " method_adjacent = \"linear_programming\"\n");
    strOut:=Concatenation(strOut, " eCommand_DD = \"glrs\"\n");
    strOut:=Concatenation(strOut, " MethodMissingI = \"Gen2\"\n");
    strOut:=Concatenation(strOut, " FileDataPoincare = \"", FilePoincare_Data, "\"\n");
    strOut:=Concatenation(strOut, " FileO = \"output.test\"\n");
    strOut:=Concatenation(strOut, " FileStepEnum = \"unset\"\n");
    strOut:=Concatenation(strOut, " Approach = \"IncrementallyAdd\"\n");
    strOut:=Concatenation(strOut, " n_iter_max = 1\n");
    strOut:=Concatenation(strOut, " n_expand = 0\n");
    strOut:=Concatenation(strOut, " Arithmetic = \"rational\"\n");
    strOut:=Concatenation(strOut, " ComputeStabilizerPermutation = T\n");
    strOut:=Concatenation(strOut, " ComputeGroupPresentation = T\n");
    strOut:=Concatenation(strOut, "/\n");
    #
    WriteStringFile(FilePoincare_Nml, strOut);
end;




TestReflectivity:=function(eRec)
    local n, FileIn, FileNml, FileOut, output, strOut, eProg, TheCommand, U, GRPmatr, ListVertNorm, isCocompact, PrefixPoincare, FilePoincare_Data, FilePoincare_Nml, ThePt, eRoot, eGen, ListGen, PrefixPoincareTot, ListGenTot, ListGenIsom, hasIsotropic, is_correct;
    n:=Length(eRec.LorMat);
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileNml:=Filename(DirectoryTemporary(), "Test.nml");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    Print("FileIn=", FileIn, "\n");
    Print("FileNml=", FileNml, "\n");
    Print("FileOut=", FileOut, "\n");
    RemoveFileIfExist(FileIn);
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileOut);
    #
    WriteMatrixFile(FileIn, eRec.LorMat);
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
    eProg:="../../src_lorentzian/LORENTZ_FundDomain_AllcockEdgewalk";
    TheCommand:=Concatenation(eProg, " ", FileNml);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileNml);
    RemoveFile(FileOut);
    ListVertNorm:=List(U.ListVertices, x->x.norm);
    isCocompact:=Maximum(ListVertNorm) < 0;
    Print("n=", n, " n_simple=", eRec.n_simple, " isCocompact=", isCocompact, " |GRPmatr|=", Order(U.GrpIsomCoxMatr), "\n");
    if GeneratePoincareFundamentalInput and isCocompact then
        iPoincare:=iPoincare + 1;
        PrefixPoincare:=Concatenation("Poincare_", String(iPoincare), "_-_", String(n), "_", String(U.n_simple));
        SaveDataToFile(PrefixPoincare, U);
        ListGen:=[];
        for eRoot in U.ListSimpleRoots
        do
            eGen:=LORENTZ_GetReflection(eRec.LorMat, eRoot);
            Add(ListGen, eGen);
        od;
        ThePt:=U.CentVect * Inverse(eRec.LorMat);
        #
        WritePoincareCase(PrefixPoincare, ThePt, ListGen);
        if Order(U.GrpIsomCoxMatr) > 1 then
            ListGenIsom:=Filtered(GeneratorsOfGroup(U.GrpIsomCoxMatr), x->x<>IdentityMat(n));
            ListGenTot:=Concatenation(ListGen, ListGenIsom);
            PrefixPoincareTot:=Concatenation(PrefixPoincare, "_Isom");
            WritePoincareCase(PrefixPoincareTot, ThePt, ListGenTot);
        fi;
        #
    fi;
    hasIsotropic:=not isCocompact;
    is_correct:=eRec.n_simple = U.n_simple;
    return rec(is_correct:=is_correct, hasIsotropic:=hasIsotropic);
end;

ListRec:=ReadAsFunction("ListReflect")();;
#ListRec:=ListRec{[1..100]};


FullTest:=function()
    local n_error, iRec, eRec, RecReply, ListIsotropicCases, eIsotropicCase;
    n_error:=0;
    iRec:=0;
    ListIsotropicCases:=[];
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        if iRec = 5345 or AllAllowed then
            RecReply:=TestReflectivity(eRec);
            if RecReply.is_correct=false then
                n_error:=n_error+1;
                return n_error;
            fi;
            if IsBound(RecReply.hasIsotropic) then
                eIsotropicCase:=rec(M:=eRec.LorMat, hasIsotropic:=RecReply.hasIsotropic);
                Add(ListIsotropicCases, eIsotropicCase);
            fi;
        fi;
        iRec:=iRec + 1;
    od;
    if GenerateIsotropicTestCases then
        SaveDataToFile("IsotropicCases", ListIsotropicCases);
    fi;
    return n_error;
end;

n_error:=FullTest();
Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

