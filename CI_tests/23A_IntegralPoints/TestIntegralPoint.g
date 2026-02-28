Read("../common.g");
Read("../access_points.g");
Print("Beginning TestIntegralPoint\n");


TestIntegralPoint:=function(FileEXT)
    local EXT, FAC, EXT_vert1, EXT_vert2, EXT_vert3, EXT_vert4;
    EXT:=ReadMatrixFile(FileEXT);
    FAC:=get_dual_desc(EXT);
    Print("|EXT|=", Length(EXT), " |FAC|=", Length(FAC), "\n");
    #
    EXT_vert1:=get_integral_interior_point(FAC, "LP_no_LLL");
    if is_error(EXT_vert1) then
        return false;
    fi;
    #
    EXT_vert2:=get_integral_interior_point(FAC, "ITER_no_LLL");
    if is_error(EXT_vert2) then
        return false;
    fi;
    #
    EXT_vert3:=get_integral_interior_point(FAC, "ITER");
    if is_error(EXT_vert3) then
        return false;
    fi;
    #
    EXT_vert4:=get_integral_interior_point(FAC, "LP");
    if is_error(EXT_vert4) then
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

