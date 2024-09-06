TestSampling:=function(FileIn, command)
    local FileOut, eProg, TheCommand, the_volume;
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    RemoveFileIfExist(FileOut);
    #
    eProg:="../../src_poly/POLY_sampling_facets";
    TheCommand:=Concatenation(eProg, " rational ", command, " ", FileIn, " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    answer:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    return true;
end;

TestRecord:=function(eRec)
    for command in eRec.l_command
    do
        test:=TestSampling(eRec.FileIn, command);
        if test=false then
            return false;
        fi;
    od;
    return true;
end;



Lcommand1:=["lp_cdd", "lrs_limited", "lp_cdd_min"];
Lcommand2:=["lp_cdd", "lrs_limited", "lp_cdd_min", "sampling"];


eRec1:=rec(FileIn:="Example_01_CUT_K333.ext", l_command:=Lcommand1);
eRec2:=rec(FileIn:="Example_02_MET_K333.ext", l_command:=Lcommand1);
eRec3:=rec(FileIn:="Example_03_CUT_K55.ext", l_command:=Lcommand1);
eRec4:=rec(FileIn:="Example_04_CUT_K144.ext", l_command:=Lcommand1);
eRec5:=rec(FileIn:="Example_05_MET_K144.ext", l_command:=Lcommand1);
eRec6:=rec(FileIn:="Example_06_Perfect_E7.ext", l_command:=Lcommand2);
eRec7:=rec(FileIn:="Example_07_CUT7.ext", l_command:=Lcommand1);
eRec8:=rec(FileIn:="Example_08_MET7.ext", l_command:=Lcommand1);
eRec9:=rec(FileIn:="Example_09_CUT8.ext", l_command:=Lcommand1);
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
if test=false then
    # Error case
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;

