Read("../common.g");
Read("../access_points.g");
Print("Beginning TestSimpleDualDesc\n");

l_arith1:=["mpq_class", "cpp_rational", "mpq_rational"];
l_arith2:=["mpq_class", "cpp_rational", "mpq_rational"];
l_arith3:=["mpq_class"];
# The integer ring arithmetic exercises the full ring computation
# (BB / LRS / DD / pd_lrs work over rings).
l_arith4:=["mpq_class", "mpz_class"];


TestSimpleDD:=function(EXT, command, n_fac)
    local dim, l_arith, arith, options, choice, FAC, eFAC, ListScal, ListIncd;
    dim:=Length(EXT[1]);
    l_arith:=l_arith4;
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
File2:="Example2_lrs_cdd_27_99_schlafli";
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
# The cross method comparison: on the same polytope, every method and
# every arithmetic must produce exactly the same facet incidence sets.
TestCrossMethod:=function(EXT, l_command, l_arith)
    local incd_ref, command, arith, vf, incd;
    incd_ref:=fail;
    for command in l_command
    do
        for arith in l_arith
        do
            vf:=get_dual_desc_incidence(EXT, command, rec(arith:=arith));
            if is_error(vf) then
                Print("Error for command=", command, " arith=", arith, "\n");
                return false;
            fi;
            incd:=Set(List(vf, Set));
            if incd_ref=fail then
                incd_ref:=incd;
            else
                if incd<>incd_ref then
                    Print("Incidence mismatch for command=", command, " arith=", arith, "\n");
                    return false;
                fi;
            fi;
        od;
    od;
    return true;
end;

GetRandomPolytope:=function(n_row, dim, amp)
    local EXT, i, j, eLine;
    while true
    do
        EXT:=[];
        for i in [1..n_row]
        do
            eLine:=[1];
            for j in [2..dim]
            do
                Add(eLine, Random([-amp..amp]));
            od;
            Add(EXT, eLine);
        od;
        if RankMat(EXT)=dim then
            return EXT;
        fi;
    od;
end;

l_command_general:=["cdd", "lrs", "pd_lrs", "beneath_beyond"];
for iter in [1..8]
do
    dim:=Random([4..6]);
    n_row:=Random([dim+3..dim+8]);
    EXT:=GetRandomPolytope(n_row, dim, 6);
    Print("Cross method iter=", iter, " |EXT|=", n_row, "/", dim, "\n");
    test:=TestCrossMethod(EXT, l_command_general, l_arith4);
    if test=false then
        n_error:=n_error + 1;
    fi;
od;

# The simplicial and near simplicial cases of small_polytopes.
l_command_small:=["small_polytopes", "lrs", "cdd"];
for iter in [1..8]
do
    dim:=Random([4..7]);
    n_row:=dim + Random([0, 1]);
    EXT:=GetRandomPolytope(n_row, dim, 6);
    Print("Small polytope iter=", iter, " |EXT|=", n_row, "/", dim, "\n");
    test:=TestCrossMethod(EXT, l_command_small, l_arith4);
    if test=false then
        n_error:=n_error + 1;
    fi;
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

