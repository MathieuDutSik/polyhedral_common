Read("../common.g");

# Runs ERDAHL_EnumeratePerfect and compares the perfect Delaunay polyhedra
# found with the expected ones, given by their number of vertex
# representatives and the rank of their isotropy lattice.
# eRec has the fields n, space ("full" or "centered"), FileCenter (for
# "centered"), method ("recursive" or "polytopes") and expected.
TestErdahl:=function(eRec)
    local FileNml, FileOut, eProg, output, TheCommand, answer, ListSignature, ePerf, eEXT, FileCenter;
    FileNml:=Filename(DirectoryTemporary(), "Test.nml");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    if IsBound(eRec.FileCenter) then
        FileCenter:=eRec.FileCenter;
    else
        FileCenter:="unset";
    fi;
    output:=OutputTextFile(FileNml, true);
    SetPrintFormattingStatus(output, false);
    AppendTo(output, "&DATA\n");
    AppendTo(output, " n = ", eRec.n, "\n");
    AppendTo(output, " space = \"", eRec.space, "\"\n");
    AppendTo(output, " FileCenter = \"", FileCenter, "\"\n");
    AppendTo(output, " method = \"", eRec.method, "\"\n");
    AppendTo(output, " FileStart = \"auto\"\n");
    AppendTo(output, " FileDualDesc = \"DualDescHeuristics.nml\"\n");
    AppendTo(output, " OutFile = \"", FileOut, "\"\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    eProg:=GetBinaryFilename("ERDAHL_EnumeratePerfect");
    TheCommand:=Concatenation(eProg, " ", FileNml);
    Print("TheCommand=", TheCommand, "\n");
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    answer:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    RemoveFile(FileNml);
    for ePerf in answer
    do
        for eEXT in ePerf.EXT
        do
            if Length(eEXT)<>eRec.n+1 or eEXT*ePerf.F*eEXT<>0 then
                Print("A vertex is not a zero of the function\n");
                return false;
            fi;
        od;
    od;
    ListSignature:=SortedList(List(answer, x->[Length(x.EXT), Length(x.L)]));
    Print("ListSignature=", ListSignature, "\n");
    if ListSignature<>SortedList(eRec.expected) then
        Print("The perfect Delaunay polyhedra are not the expected ones, expected=", eRec.expected, "\n");
        return false;
    fi;
    return true;
end;

ConcludeErdahl:=function(test)
    Print("test=", test, "\n");
    CI_Decision_Reset();
    if test=false then
        # Error case
        Print("Error case\n");
    else
        # No error case
        Print("Normal case\n");
        CI_Write_Ok();
    fi;
end;
