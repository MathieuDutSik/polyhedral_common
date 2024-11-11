Read("../common.g");
Print("Doing the computation of skeleton and groups\n");

TestGroupSkeleton:=function(eRec)
    local FileInputNml, FileGrpOut, FileResultOut, strOut, output, eProg, TheCommand, GRPout;
    FileInputNml:=Filename(DirectoryTemporary(), "TestInput.nml");
    FileGrpOut:=Filename(DirectoryTemporary(), "TestGrp.out");
    FileResultOut:=Filename(DirectoryTemporary(), "TestResult.out");
    RemoveFileIfExist(FileInputNml);
    RemoveFileIfExist(FileGrpOut);
    RemoveFileIfExist(FileResultOut);
    #
    strOut:="&PROC\n";
    strOut:=Concatenation(strOut, " FACfile = \"", eRec.FileExt, "\"\n");
    strOut:=Concatenation(strOut, " EXTfile = \"unset.ext\"\n");
    strOut:=Concatenation(strOut, " GRPfile = \"", eRec.FileGrp, "\"\n");
    strOut:=Concatenation(strOut, " OUTfile =\"", FileResultOut, "\"\n");
    strOut:=Concatenation(strOut, " ComputeTotalNumberFaces = T\n");
    strOut:=Concatenation(strOut, " method_spann = \"LinearProgramming\"\n");
    strOut:=Concatenation(strOut, " method_final = \"all\"\n");
    strOut:=Concatenation(strOut, " Arithmetic = \"rational\"\n");
    strOut:=Concatenation(strOut, " LevSearch = ", String(eRec.LevSearch), "\n");
    strOut:=Concatenation(strOut, "/\n");
    strOut:=Concatenation(strOut, "\n");
    strOut:=Concatenation(strOut, "&GROUP\n");
    strOut:=Concatenation(strOut, " ComputeAutGroup = T\n");
    strOut:=Concatenation(strOut, " OutFormat = \"CPP\"\n");
    strOut:=Concatenation(strOut, " GileGroup = \"", FileGrpOut, "\"\n");
    strOut:=Concatenation(strOut, "/\n");
    #
    output:=OutputTextFile(FileInputNml, true);
    WriteAll(output, strOut);
    CloseStream(output);
    #
    eProg:="../../src_poly/POLY_FaceLatticeGen";
    TheCommand:=Concatenation(eProg, " ", FileInputNml);
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
CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

