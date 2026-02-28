Read("../common.g");
Read("../access_points.g");
Print("Beginning TestSimpleDualDesc\n");


TestSimpleDD:=function(EXT, command, n_fac)
    local dim, FileI, FileO, arith, choice, eProg, TheCommand, FAC, eFAC, ListScal, ListIncd;
    dim:=Length(EXT[1]);
    FAC:=get_dual_desc(EXT, command);
    if is_error(FAC) then
        return false;
    fi;
    if Length(FAC)<>n_fac then
        Print("Incorrect number of facets. That qualifies as a fail\n");
        return false;
    fi;
    for eFAC in FAC
    do
        ListScal:=List(EXT, x->x*eFAC);
        if Minimum(ListScal) < 0 then
            Print("Find a negative scalar product, a fail for sure\n");
            return false;
        fi;
        ListIncd:=Filtered([1..Length(EXT)], x->ListScal[x]=0);
        if RankMat(EXT{ListIncd}) <> dim-1 then
            Print("The rank is not correct. A fail\n");
            return false;
        fi;
    od;
    return true;
end;

File1:="Example1_pd_lrs_1084_26";
File2:="Example2_lrs_cdd_27_99";
File3:="Example3_48_11432";
ListFiles:=[File1, File2, File3];

n_error:=0;
for iFile in [1..Length(ListFiles)]
do
    eFile:=ListFiles[iFile];
    eRec:=ReadAsFunction(eFile)();
    Print("iFile=", iFile, " |EXT|=", Length(eRec.EXT), "/", Length(eRec.EXT[1]), " n_fac=", eRec.n_fac, "\n");
    for i_command in [1..Length(eRec.commands)]
    do
        command:=eRec.commands[i_command];
        Print("   i_command=", i_command, " command=", command, "\n");
        test:=TestSimpleDD(eRec.EXT, command, eRec.n_fac);
        if test=false then
            n_error:=n_error + 1;
        fi;
    od;
od;
Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;

