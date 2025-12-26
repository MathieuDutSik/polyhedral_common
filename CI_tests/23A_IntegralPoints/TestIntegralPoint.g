Read("../common.g");
Print("Beginning TestIntegralPoint\n");

get_integral_interior_point:=function(FAC, method)
    local TmpDir, FileI, FileO, eProg, TheCommand, EXTint;
    TmpDir:=DirectoryTemporary();
    FileI:=Filename(TmpDir, "Test.fac");
    FileO:=Filename(TmpDir, "Test.vertint");
    WriteMatrixFile(FileI, FAC);
    eProg:="../../src_poly/POLY_IntegralPoints";
    TheCommand:=Concatenation(eProg, " gmp ", method, " ", FileI, " CPP ", FileO);
    Exec(TheCommand);
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    EXTint:=ReadMatrixFile(FileO);
    RemoveFile(FileI);
    RemoveFile(FileO);
    return EXTint;
end;

get_dual_desc:=function(EXT)
    local TmpDir, FileEXT, FileFAC, eProg, command, TheCommand, FAC;
    TmpDir:=DirectoryTemporary();
    FileEXT:=Filename(TmpDir, "Test.out");
    FileFAC:=Filename(TmpDir, "Test.fac");
    WriteMatrixFile(FileEXT, EXT);
    eProg:="../../src_poly/POLY_dual_description";
    command:="cdd";
    TheCommand:=Concatenation(eProg, " rational  ", command, " CPP ", FileEXT, " ", FileFAC);
    Exec(TheCommand);
    if IsExistingFile(FileFAC)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    FAC:=ReadMatrixFile(FileFAC);
    RemoveFile(FileEXT);
    RemoveFile(FileFAC);
    return FAC;
end;



TestIntegralPoint:=function(FileEXT)
    local EXT, FAC, EXT_vert1, EXT_vert2, EXT_vert3, EXT_vert4;
    EXT:=ReadMatrixFile(FileEXT);
    FAC:=get_dual_desc(EXT);
    Print("|EXT|=", Length(EXT), " |FAC|=", Length(FAC), "\n");
    #
    EXT_vert1:=get_integral_interior_point(FAC, "LP_no_LLL");
    if EXT_vert1=false then
        Print("Failure for LP_no_LLL\n");
        return false;
    fi;
    #
    EXT_vert2:=get_integral_interior_point(FAC, "ITER_no_LLL");
    if EXT_vert2=false then
        Print("Failure for ITER_no_LLL\n");
        return false;
    fi;
    #
    EXT_vert3:=get_integral_interior_point(FAC, "ITER");
    if EXT_vert3=false then
        Print("Failure for ITER\n");
        return false;
    fi;
    #
    EXT_vert4:=get_integral_interior_point(FAC, "LP");
    if EXT_vert4=false then
        Print("Failure for LP\n");
        return false;
    fi;
    #
    if Set(EXT_vert1)<>Set(EXT_vert2) then
        Print("Inconsistency between LP_no_LLL and ITER_no_LLL\n");
        return false;
    fi;
    #
    if Set(EXT_vert2)<>Set(EXT_vert3) then
        Print("Inconsistency between ITER_no_LLL and ITER\n");
        return false;
    fi;
    #
    if Set(EXT_vert3)<>Set(EXT_vert4) then
        Print("Inconsistency between ITER and LP\n");
        return false;
    fi;
    return true;
end;

ListFile:=["G6.txt", "ER35.ext"];
for i in [1..27]
do
    eFile:=Concatenation("Perfect8_", String(i));
    Add(ListFile, eFile);
od;

FullTest:=function()
    local iFile, eFile, test;
    iFile:=0;
    for eFile in ListFile
    do
        Print("iFile=", iFile, " / ", Length(ListFile), "\n");
        test:=TestIntegralPoint(eFile);
        if test=false then
            return false;
        fi;
        iFile:=iFile + 1;
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

