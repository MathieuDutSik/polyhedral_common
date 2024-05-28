Read("../common.g");

TestGroupSkeleton:=function(eRec)
    local n, FileIn, FileNml, FileOut, output, eProg, TheCommand, U, GRPmatr, ListVertNorm, isCocompact, PrefixPoincare, FilePoincare_Data, FilePoincare_Nml, ThePt, eRoot, eGen, ListGen, PrefixPoincareTot, ListGenTot, ListGenIsom, hasIsotropic, is_correct;
    nbRow:=Length(eRec.FAC);
    dim:=Length(eRec.FAC[1]);
    FileInputNml:=Filename(DirectoryTemporary(), "TestInput.nml");
    FileGrpOut:=Filename(DirectoryTemporary(), "TestGrp.out");
    FileResultOut:=Filename(DirectoryTemporary(), "TestResult.out");
    RemoveFileIfExist(FileInputNml);
    RemoveFileIfExist(FileGrpOut);
    RemoveFileIfExist(FileResultOut);
    #
    output:=OutputTextFile(FileNml, true);
    AppendTo(output, "&PROC\n");
    AppendTo(output, " FACfile = \"", eRec.FileExt, "\"\n");
    AppendTo(output, " EXTfile = \"unset.ext\"\n");
    AppendTo(output, " GRPfile = \"", eRec.FileGrp, "\"\n");
    AppendTo(output, " OUTfile =\"", FileResultOut, "\"\n");
    AppendTo(output, " ComputeTotalNumberFaces = T\n");
    AppendTo(output, " method_spann = \"LinearProgramming\"\n");
    AppendTo(output, " method_final = \"all\"\n");
    AppendTo(output, " Arithmetic = \"rational\"\n");
    AppendTo(output, " LevSearch = ", eRec.LevSearch, "\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&GROUP\n");
    AppendTo(output, " ComputeAutGroup = T\n");
    AppendTo(output, " OutFormat = \"CPP\"\n");
    AppendTo(output, " GileGroup = \"", FileGrpOut, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    eProg:="../../src_poly/POLY_FaceLatticeGen";
    TheCommand:=Concatenation(eProg, " ", FileNml);
    Exec(TheCommand);
    if IsExistingFile(FileGrpOut)=false or IsExistingFile(FileResultOut) then
        Print("The output file is not existing. That qualifies as a fail\n");
        return rec(is_correct:=false);
    fi;
    GRPout:=ReadAsFunction(FileGrpOut)();
    if Order(GRPout)<>eRec.GRPoutOrder then
        Print("The group order is not what we expect\n");
        return rec(is_correct:=false);
    fi;
    return rec(is_correct:=true);
end;

ListRec:=[
 rec(FileExt:="G6.ext", FileGrp:="G6.grp", LevSearch:=2, GRPoutOrder:=51840),
 rec(FileExt:="Pent.ext", FileGrp:="Pent.grp", LevSearch:=1, GRPoutOrder:=10)
];

FullTest:=function()
    local n_error, iRec, eRec, RecReply, ListIsotropicCases, eIsotropicCase;
    n_error:=0;
    iRec:=0;
    ListIsotropicCases:=[];
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        RecReply:=TestGroupSkeleton(eRec);
        if RecReply.is_correct=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;

n_error:=FullTest();
if n_error > 0 then
    # Error case
    GAP_EXIT_CODE(1);
else
    # No error case
    GAP_EXIT_CODE(0);
fi;

