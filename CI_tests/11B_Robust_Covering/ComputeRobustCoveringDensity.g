Read("../common.g");
Read("../access_points.g");
Print("Beginning ComputeRobustCoveringDensity\n");


TestComputation:=function(eRec)
    local result;
    result:=get_robust_lattice_covering_density(eRec.eMat);
    if is_error(result) then
        return false;
    fi;
    if result.TheCov<>eRec.TheCov then
        return false;
    fi;
    Print("result=", result, "\n");
    return true;
end;

eRec1:=rec(name:="A2",
           eMat:=ClassicalSporadicLattices("A2"),
           TheCov:=8/3);
ListRec:=[eRec1];

FullTest:=function()
    local iRec, eRec, reply;
    iRec:=0;
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        reply:=TestComputation(eRec);
        Print("reply(B)=", reply, "\n");
        if reply=false then
            Print("Error occurred\n");
            return false;
        fi;
        iRec:=iRec + 1;
    od;
    Print("No error occurred\n");
    return true;
end;

NestFunction:=function()
    local result;
    result:=FullTest();
    Print("2: result=", result, "\n");
    CI_Decision_Reset();
    if result=false then
        # Error case
        Print("Error case\n");
    else
        # No error case
        Print("Normal case\n");
        CI_Write_Ok();
    fi;
    CI_PrintExistConclusion();
end;

NestFunction();


