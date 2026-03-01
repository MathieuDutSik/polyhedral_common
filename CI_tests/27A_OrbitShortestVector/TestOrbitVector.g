Read("../common.g");
Read("../access_points.g");


TestOrbitShortest:=function(eRec)
    local GramMat, U;
    GramMat:=GetGramMatrixFromString(eRec.name);
    U:=get_orbit_shortest(GramMat);
    if is_error(U) then
        return false;
    fi;
    Print("|U.vf|=", Length(U.vf), " n_orbit=", eRec.n_orbit, "\n");
    if Length(U.vf)<>eRec.n_orbit then
        Print("The number of orbits is incorrect\n");
        return false;
    fi;
    return true;
end;

ListRec:=[];
Add(ListRec, rec(name:="A2", n_orbit:=1));
Add(ListRec, rec(name:="A3", n_orbit:=2));

FullTest:=function()
    local iRec, eRec, test;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        test:=TestOrbitShortest(eRec);
        if test=false then
            return false;
        fi;
        iRec:=iRec + 1;
    od;
    return true;
end;

result:=FullTest();
Print("result=", result, "\n");

CI_Decision_Reset();
if result = false then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;

