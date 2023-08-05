Print("Beginning TestSimpleDualDesc\n");

ReadMatrixFile:=function(eFile)
    local ListLines, TheMat, nRow, nCol, iRow;
    ListLines:=ReadVectorFile(eFile);
    TheMat:=ListLines{[2..Length(ListLines)]};
    nRow:=ListLines[1][1];
    nCol:=ListLines[1][2];
    if Length(TheMat)<>nRow then
        Error("Error in number of rows");
    fi;
    for iRow in [1..nRow]
    do
        if Length(TheMat[iRow]) <> nCol then
            Error("Error in number of cols");
        fi;
    od;
    return TheMat;
end;

WriteMatrixFile:=function(eFile, EXT)
    local output, eEXT;
    output:=OutputTextFile(eFile, true);
    AppendTo(output, Length(EXT), " ", Length(EXT[1]), "\n");
    for eEXT in EXT
    do
        WriteVector(output, eEXT);
    od;
    CloseStream(output);
end;


TestSimpleDD:=function(TheRec)
    local dim, FileI, FileO, arith, choice, eProg, TheCommand, FAC, eFAC, ListScal, ListIncd;
    dim:=Length(TheRec.EXT[1]);
    FileI:="Test.in";
    FileO:="Test.out";
    #
    WriteMatrixFile(FileIn, TheRec.EXT);
    #
    arith:="rational";
    choice:="CPP";
    eProg:="../../src_lorentzian/LORENTZ_FundDomain_AllcockEdgewalk";
    TheCommand:=Concatenation(eProg, " ", arith, " ", TheRec.command, " ", choice, " ", FileI, " ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    FAC:=ReadMatrixFile(FileO);
    RemoveFile(FileI);
    RemoveFile(FileO);
    if Length(FAC)<>TheRec.n_fac then
        Print("Incorrect number of facets. That qualifies as a fail\n");
        return false;
    fi;
    for eFAC in FAC
    do
        ListScal:=List(TheRec.EXT, x->x*eFAC);
        if Minimum(ListScal) < 0 then
            Print("Find a negative scalar product, a fail for sure\n");
            return false;
        fi;
        ListIncd:=Filtered([1..Length(EXT)], x->ListScal[x]=0);
        if RankMat(TheRec.EXT{ListIncd}) <>dim-1 then
            Print("The rank is not correct. A fail\n");
            return false;
        fi;
    od;
    return true;
end;


ListFiles:=["Example1_pd_lrs_1084_26"];


for eFile in ListFiles
do
    eRec:=ReadAsFunction(eFile)();
    test:=TestSimpleDD(eRec);
    if test=false then
        # Error case
        GAP_EXIT_CODE(1);
    fi;
od;
# No error case
GAP_EXIT_CODE(0);
