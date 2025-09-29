Read("../common.g");
Print("Beginning TestReflectivity\n");

GeneratePoincareFundamentalInput:=false;
GenerateIsotropicTestCases:=true;

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



GetReflectivityInformation:=function(eRec)
    local n, TmpDir, FileIn, FileNml, FileOut, strOut, eProg, TheCommand, U;
    n:=Length(eRec.LorMat);
    TmpDir:=DirectoryTemporary();
    FileIn:=Filename(TmpDir, "Test.in");
    FileNml:=Filename(TmpDir, "Test.nml");
    FileOut:=Filename(TmpDir, "Test.out");
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
    return U;
end;




TestReflectivity:=function(eRec)
    local U, ListVertNorm, isCocompact, n, PrefixPoincare, ListGen, ListGenTot, eRoot, ThePt, hasIsotropic, is_correct;
    U:=GetReflectivityInformation(eRec);
    ListVertNorm:=List(U.ListVertices, x->x.norm);
    isCocompact:=Maximum(ListVertNorm) < 0;
    n:=Length(eRec.LorMat);
    Print("n=", n, " n_simple=", eRec.n_simple, " isCocompact=", isCocompact, " |GRPmatr|=", Order(U.GrpIsomCoxMatr), "\n");
    hasIsotropic:=not isCocompact;
    is_correct:=eRec.n_simple = U.n_simple;
    return rec(is_correct:=is_correct, hasIsotropic:=hasIsotropic);
end;


GeneratePoincareInput:=function(ListRecInput)
    local iPoincare, eRec, U, ListVertNorm, isCocompact, n, PrefixPoincare, PrefixPoincareTot, ListGen, ThePt, ListGenTot, ListGenIsom;
    iPoincare:=0;
    for eRec in ListRecInput
    do
        U:=GetReflectivityInformation(eRec);
        ListVertNorm:=List(U.ListVertices, x->x.norm);
        isCocompact:=Maximum(ListVertNorm) < 0;
        n:=Length(eRec.LorMat);
        Print("n=", n, " n_simple=", eRec.n_simple, " isCocompact=", isCocompact, " |GRPmatr|=", Order(U.GrpIsomCoxMatr), "\n");
        if isCocompact then
            iPoincare:=iPoincare + 1;
            PrefixPoincare:=Concatenation("Poincare_", String(iPoincare), "_-_", String(n), "_", String(U.n_simple));
            SaveDataToFile(PrefixPoincare, U);
            ListGen:=List(U.ListSimpleRoots, x->LORENTZ_GetReflection(eRec.LorMat, x));
            ThePt:=U.CentVect * Inverse(eRec.LorMat);
            #
            WritePoincareCase(PrefixPoincare, ThePt, ListGen);
            if Order(U.GrpIsomCoxMatr) > 1 then
                ListGenIsom:=Filtered(GeneratorsOfGroup(U.GrpIsomCoxMatr), x->x<>IdentityMat(n));
                ListGenTot:=Concatenation(ListGen, ListGenIsom);
                PrefixPoincareTot:=Concatenation(PrefixPoincare, "_Isom");
                WritePoincareCase(PrefixPoincareTot, ThePt, ListGenTot);
            fi;
        fi;
    od;

end;



GenerateFamilyDomains:=function(ListRecInput)
    local ListEXT, ListDim, nRec, iRec, eRec, U, dim, n_simple, eColl, FileStr, IsFirst, ePair;
    ListEXT:=[];
    ListDim:=[];
    nRec:=Length(ListRecInput);
    for iRec in [1..nRec]
    do
        eRec:=ListRecInput[iRec];
        U:=GetReflectivityInformation(eRec);
        dim:=Length(eRec.LorMat);
        n_simple:=Length(U.ListSimpleRoots);
        Print("iRec=", iRec, " / ", nRec, " dim=", dim, " n_simple=", n_simple, "\n");
        Add(ListDim, dim);
        Add(ListEXT, U.ListSimpleRoots);
    od;
    eColl:=Collected(ListDim);
    FileStr:="ListSimpleRootSystem";
    IsFirst:=true;
    for ePair in eColl
    do
        if IsFirst = false then
            FileStr:=Concatenation(FileStr, "_X");
        fi;
        FileStr:=Concatenation(FileStr, "_", String(ePair[1]), "_", String(ePair[2]));
        IsFirst:=false;
    od;
    Print("eColl=", eColl, "\n");
    Print("FileStr=", FileStr, "\n");
    SaveDataToFile(FileStr, ListEXT);
end;





ListRec:=ReadAsFunction("ListReflect")();;
ListRec45:=Filtered(ListRec, x->Length(x.LorMat)>=4 and Length(x.LorMat) <=5);
#ListRec:=ListRec{[1..100]};


GenerateFamilyDomains(ListRec45);


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

