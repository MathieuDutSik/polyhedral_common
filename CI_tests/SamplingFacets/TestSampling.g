Read("../common.g");
Print("Beginning TestSampling\n");

TestSampling:=function(FileIn, command)
    local FileOut, TheLen, eProg, TheCommand, answer;
    Print("Running FileIn=", FileIn, " with command=", command, "\n");
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    #
    eProg:="../../src_poly/POLY_sampling_facets";
    TheCommand:=Concatenation(eProg, " rational ", command, " ", FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    answer:=ReadAsFunction(FileOut)();
    TheLen:=Length(answer);
    Print("TheLen=", TheLen, "\n");
    if TheLen = 0 then
        return false;
    fi;
    RemoveFile(FileOut);
    return true;
end;

TestRecord:=function(eRec)
    local EXT, command, test;
    EXT:=ReadMatrixFile(eRec.FileIn);
    Print("\n");
    Print("|EXT|=", Length(EXT), " / ", Length(EXT[1]), "\n");
    for command in eRec.l_command
    do
        test:=TestSampling(eRec.FileIn, command);
        if test=false then
            return false;
        fi;
    od;
    return true;
end;



Lcommand1:=["lp_cdd", "lp_cdd_min"];
Lcommand2:=["lp_cdd", "lrs_limited", "lp_cdd_min"];
Lcommand3:=["lp_cdd", "lrs_limited", "lp_cdd_min", "sampling"];


eRec1:=rec(FileIn:="Example_01_CUT_K333.ext", l_command:=Lcommand1);
eRec2:=rec(FileIn:="Example_02_MET_K333.ext", l_command:=Lcommand1);
eRec3:=rec(FileIn:="Example_03_CUT_K55.ext", l_command:=Lcommand1);
eRec4:=rec(FileIn:="Example_04_CUT_K144.ext", l_command:=Lcommand1);
eRec5:=rec(FileIn:="Example_05_MET_K144.ext", l_command:=Lcommand1);
eRec6:=rec(FileIn:="Example_06_Perfect_E7.ext", l_command:=Lcommand2);
eRec7:=rec(FileIn:="Example_07_CUT7.ext", l_command:=Lcommand2);
eRec8:=rec(FileIn:="Example_08_MET7.ext", l_command:=Lcommand3);
eRec9:=rec(FileIn:="Example_09_CUT8.ext", l_command:=Lcommand2);
ListRec:=[eRec1, eRec2, eRec3, eRec4, eRec5, eRec6, eRec7, eRec8, eRec9];

FullTest:=function()
    local iRec, eRec, test;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        test:=TestRecord(eRec);
        if test=false then
            return false;
        fi;
        iRec:=iRec + 1;
    od;
    return true;
end;

test:=FullTest();

CI_Decision_Reset();
if test=false then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

