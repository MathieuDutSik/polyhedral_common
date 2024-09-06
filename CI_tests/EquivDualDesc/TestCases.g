TestDualDesc:=function(eRec)
    local prefix, FileOut, eProg, TheCommand, the_volume;
    prefix
    #
    eProg:="../../../src_poly/POLY_SerialDualDesc";
    TheCommand:=Concatenation("(cd ", prefix, " && ", eProg, " input.nml)");
    Exec(TheCommand);
    FileOut:=Concatenation(prefix, "/orbits");
    LOrb:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    if Length(LOrb)<>eRec.n_orbit then
        Print("incoherence of the result\n");
        return false;
    fi;
    return true;
end;

eRec1:=rec(FileIn:="24cell", n_orbit:=1);
eRec2:=rec(FileIn:="ContactE8", n_orbit:=2);
eRec3:=rec(FileIn:="CUT_K7", n_orbit:=11);
eRec4:=rec(FileIn:="CUT_K8", n_orbit:=147);
ListRec:=[eRec1, eRec2, eRec3, eRec4];

FullTest:=function()
    local iRec, eRec, test;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        test:=TestDualDesc(eRec);
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
