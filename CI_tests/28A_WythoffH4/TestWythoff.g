Read("../common.g");
Print("Beginning TestWythoff\n");

TestFCT:=function()
    local eProg, FileI, TheCommand, FileO, ListOrb;
    #
    eProg:="../../src_dualdesc/POLY_SerialDualDesc";
    FileI:="FacetsWythoffH4.nml";
    TheCommand:=Concatenation(eProg, " ", FileI);
    Exec(TheCommand);
    FileO:="ListOrbit";
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    ListOrb:=ReadAsFunction(FileO)();
    if Length(ListOrb)<>4 then
        Print("Incorrect number of facets. That qualifies as a fail\n");
        return false;
    fi;
    return true;
end;

test:=TestFCT();
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
