Read("../common.g");
Read("../access_points.g");
Print("Doing the computation of skeleton and groups\n");

TestGroupSkeleton:=function(eRec)
    local GRP, pair_ans;
    GRP:=get_grp_automorphy(eRec.EXT);
    pair_ans:=get_group_skeleton(eRec.EXT, GRP.GAPperm, eRec.LevSearch);
    if is_error(pair_ans) then
        return false;
    fi;
    if Order(pair_ans[1]) <> eRec.GRPoutOrder then
        return false;
    fi;
    return true;
end;

EXT_G6:=ReadMatrixFile("G6.ext");
rec1:=rec(EXT:=EXT_G6, LevSearch:=2, GRPoutOrder:=51840);

EXT_Pent:=ReadMatrixFile("Pent.ext");
rec2:=rec(EXT:=EXT_Pent, LevSearch:=1, GRPoutOrder:=10);

ListRec:=[rec1, rec2];

FullTest:=function()
    local n_error, iRec, eRec, test, ListIsotropicCases, eIsotropicCase;
    n_error:=0;
    iRec:=0;
    ListIsotropicCases:=[];
    for eRec in ListRec
    do
        Print("iRec=", iRec, " / ", Length(ListRec), "\n");
        test:=TestGroupSkeleton(eRec);
        if test=false then
            n_error:=n_error+1;
            return n_error;
        fi;
        iRec:=iRec + 1;
    od;
    return n_error;
end;

n_error:=FullTest();
CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

