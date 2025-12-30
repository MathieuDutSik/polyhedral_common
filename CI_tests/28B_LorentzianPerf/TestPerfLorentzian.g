Read("../common.g");
Print("Beginning TestPerfectEnumeration\n");

GetNumberPerfectLorentzian:=function(eRec, choice)
    local n, FileIn, FileOut, FileNml, eProg, TheCommand, U, V, eNorm, eNormSqr, strOut, LorMat;
    n:=Length(eRec.M);
    FileIn:=Filename(DirectoryTemporary(), "Test.in");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    FileNml:=Filename(DirectoryTemporary(), "Test.nml");
    RemoveFileIfExist(FileIn);
    RemoveFileIfExist(FileOut);
    #
    if choice = "isotropic" and eRec.has_isotropic = false then
        return fail;
    fi;
    #
    strOut:="&DATA\n";
    strOut:=Concatenation(strOut, "  LorMatFile = \"", FileIn, "\"\n");
    strOut:=Concatenation(strOut, "  Option = \"", choice, "\"\n");
    strOut:=Concatenation(strOut, "/\n");
    strOut:=Concatenation(strOut, "\n");
    strOut:=Concatenation(strOut, "&SYSTEM\n");
    strOut:=Concatenation(strOut, "  OutFormat = \"NumberGAP\"\n");
    strOut:=Concatenation(strOut, "  OutFile = \"", FileOut, "\"\n");
    strOut:=Concatenation(strOut, "/\n");
    #
    WriteStringFile(FileNml, strOut);
    #
    LorMat:= - eRec.M;
    WriteMatrixFile(FileIn, LorMat);
    #
    eProg:="../../src_lorentzian/LORENTZ_MPI_PerfectLorentzian";
    TheCommand:=Concatenation(eProg, " ", FileNml);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        Print(NullMat(5));
        return rec(is_correct:=false);
    fi;
    U:=ReadAsFunction(FileOut)();
    RemoveFile(FileIn);
    RemoveFile(FileNml);
    RemoveFile(FileOut);
    return U.nb;
end;




ListRec:=ReadAsFunction("../DATA/IsotropicCases")();;
ListDim3:=Filtered(ListRec, x->Length(x.M) = 3);
ListDim4:=Filtered(ListRec, x->Length(x.M) = 4);
ListDim5:=Filtered(ListRec, x->Length(x.M) = 5);
ListDim5gt:=Filtered(ListRec, x->Length(x.M) > 5);


ListWorkRec:=[];
Append(ListWorkRec, ListDim3{[1..10]});
Append(ListWorkRec, ListDim4{[1..10]});
Append(ListWorkRec, ListDim5{[1..10]});
Append(ListWorkRec, ListDim5gt{[1..10]});


ListChoices:=["isotropic", "total"];
#ListChoices:=["total"];
#ListChoices:=["isotropic"];


GetListResult:=function(ListRec)
    local ListEntry, nRec, iRec, eRec, choice, nPerf, eEntry;
    ListEntry:=[];
    nRec:=Length(ListRec);
    for iRec in [1..Length(ListRec)]
    do
        eRec:=ListRec[iRec];
        for choice in ListChoices
        do
            Print("iRec=", iRec, " / ", nRec, " choice=", choice, "\n");
            nPerf:=GetNumberPerfectLorentzian(eRec, choice);
            eEntry:=[eRec, choice, nPerf];
            Add(ListEntry, eEntry);
        od;
    od;
    return ListEntry;
end;

ListEntry:=GetListResult(ListWorkRec);
CI_Decision_Reset();

FileSave:="Result_Enumeration";
if IsExistingFile(FileSave)=false then
    SaveDataToFile(FileSave, ListEntry);
else
    ListResult:=ReadAsFunction(FileSave)();
    if ListResult<>ListEntry then
        Print("Error case\n");
    else
        Print("Normal case\n");
        CI_Write_Ok();
    fi;
fi;

Print("|ListDim3|=", Length(ListDim3), "\n");
Print("|ListDim4|=", Length(ListDim4), "\n");
Print("|ListDim5|=", Length(ListDim5), "\n");
Print("|ListDim5gt|=", Length(ListDim5gt), "\n");
