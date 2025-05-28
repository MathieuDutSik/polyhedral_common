Read("../common.g");
Print("Beginning TestPerfectEnumeration\n");

GetNumberPerfectLorentzian:=function(eRec, choice)
    local n, FileIn, FileOut, FileNml, eProg, TheCommand, U, V, eNorm, eNormSqr, output, LorMat;
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
    output:=OutputTextFile(FileNml, true);
    AppendTo(output, "&DATA\n");
    AppendTo(output, "  LorMatFile = \"", FileIn, "\"\n");
    AppendTo(output, "  Option = \"", choice, "\"\n");
    AppendTo(output, "  OutFormat = \"NumberGAP\"\n");
    AppendTo(output, "  OutFile = \"", FileOut, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);
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


ListWorkRec:=ListDim5{[1..20]};


#ListChoices:=["isotropic", "total"];
#ListChoices:=["total"];
ListChoices:=["isotropic"];


GetListPerf:=function(ListRec)
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

ListEntry:=GetListPerf(ListWorkRec);

