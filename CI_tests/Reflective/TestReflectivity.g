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


WritePoincareCase:=function(PrefixPoincare, ThePt, ListGen)
    local FilePoincare_Data, FilePoincare_Nml, output, eGen;
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
    output:=OutputTextFile(FilePoincare_Nml, true);
    AppendTo(output, "&PROC\n");
    AppendTo(output, " method_adjacent = \"linear_programming\"\n");
    AppendTo(output, " eCommand_DD = \"glrs\"\n");
    AppendTo(output, " MethodMissingI = \"Gen2\"\n");
    AppendTo(output, " FileDataPoincare = \"", FilePoincare_Data, "\"\n");
    AppendTo(output, " FileO = \"output.test\"\n");
    AppendTo(output, " FileStepEnum = \"unset\"\n");
    AppendTo(output, " Approach = \"IncrementallyAdd\"\n");
    AppendTo(output, " n_iter_max = 1\n");
    AppendTo(output, " n_expand = 0\n");
    AppendTo(output, " Arithmetic = \"rational\"\n");
    AppendTo(output, " ComputeStabilizerPermutation = T\n");
    AppendTo(output, " ComputeGroupPresentation = T\n");
    AppendTo(output, "/\n");
    CloseStream(output);
end;




TestReflectivity:=function(eRec)
    local n, FileIn, FileNml, FileOut, output, eProg, TheCommand, U, GRPmatr, ListVertNorm, isCocompact, PrefixPoincare, FilePoincare_Data, FilePoincare_Nml, ThePt, eRoot, eGen, ListGen, PrefixPoincareTot, ListGenTot, ListGenIsom, hasIsotropic, is_correct;
    n:=Length(eRec.LorMat);
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileNml:=Filename(DirectoryTemporary(), "Test.nml");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    RemoveFileIfExist(FileIn);
    RemoveFileIfExist(FileNml);
    RemoveFileIfExist(FileOut);
    #
    WriteMatrixFile(FileIn, eRec.LorMat);
    #
    output:=OutputTextFile(FileNml, true);
    AppendTo(output, "&PROC\n");
    AppendTo(output, " FileLorMat = \"", FileIn, "\"\n");
    AppendTo(output, " OptionInitialVertex = \"vinberg\"\n");
    AppendTo(output, " OutFormat = \"GAP\"\n");
    AppendTo(output, " FileOut = \"", FileOut, "\"\n");
    AppendTo(output, " OptionNorms = \"all\"\n");
    AppendTo(output, " EarlyTerminationIfNotReflective = T\n");
    AppendTo(output, " ComputeAllSimpleRoots = T\n");
    AppendTo(output, "/\n");
    CloseStream(output);
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
        RecReply:=TestReflectivity(eRec);
        if RecReply.is_correct=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        if IsBound(RecReply.hasIsotropic) then
            eIsotropicCase:=rec(M:=eRec.LorMat, hasIsotropic:=RecReply.hasIsotropic);
            Add(ListIsotropicCases, eIsotropicCase);
        fi;
        iRec:=iRec + 1;
    od;
    if GenerateIsotropicTestCases then
        SaveDataToFile("IsotropicCases", ListIsotropicCases);
    fi;
    return n_error;
end;

n_error:=FullTest();
if n_error > 0 then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

