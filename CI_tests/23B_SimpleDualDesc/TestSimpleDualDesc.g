Read("../common.g");
Read("../access_points.g");
Print("Beginning TestSimpleDualDesc\n");

#l_arith:=["mpq_class", "safe_rational", "cpp_rational", "mpq_rational"];
l_arith:=["mpq_class", "cpp_rational", "mpq_rational"];


TestSimpleDD:=function(EXT, command, n_fac)
    local dim, arith, options, choice, FAC, eFAC, ListScal, ListIncd;
    dim:=Length(EXT[1]);
    for arith in l_arith
    do
        options:=rec(print_info:=true, arith:=arith);
        FAC:=get_dual_desc(EXT, command, options);
        if is_error(FAC) then
            return false;
        fi;
        if Length(FAC)<>n_fac then
            Print("|FAC|=", Length(FAC), " n_fac=", n_fac, "\n");
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
    od;
    return true;
end;

File1:="Example1_pd_lrs_1084_26";
File2:="Example2_lrs_cdd_27_99";
File3:="Example3_48_11432";
File4:="Example4_cdd_lrs_CUTP6";
File5:="Example5_cdd_lrs_tsp6";
File6:="Example6_cdd_lrs_METP6";
File7:="Example7_cdd_lrs_mit41_16";
ListFiles:=[File1, File2, File3, File4, File5, File6, File7];

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

