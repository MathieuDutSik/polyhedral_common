Read("../common.g");
Print("Beginning TestCases\n");

TestDualDesc:=function(eRec)
    local prefix, FileOut, eProg, TheCommand, LOrb;
    prefix:=eRec.prefix;
    #
    eProg:="../../../src_dualdesc/POLY_SerialDualDesc";
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

ListRec:=[];
Add(ListRec, rec(prefix:="24cell", n_orbit:=1));
Add(ListRec, rec(prefix:="ContactE8", n_orbit:=2));
Add(ListRec, rec(prefix:="CUT_K7", n_orbit:=11));
#Add(ListRec, rec(prefix:="CUT_K8", n_orbit:=147));

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
CI_Decision_Reset();
if test=false then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
