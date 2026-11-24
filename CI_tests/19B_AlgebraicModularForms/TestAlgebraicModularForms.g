Read("../common.g");
Read("../access_points.g");
Print("Beginning of TestAlgebraicModularForms\n");

# The Hecke matrices of GENUS_MPI_AlgebraicModularForms on the cases of
# cases.g. For every case: the result is the expected one, every row of a
# Hecke matrix sums to the number of neighbors, the matrices are
# self-adjoint for the inner product of the states, |O| / |Aut(L_j)| for
# the state (L_j, O), and the matrices at two primes commute. With the
# environment variable AMF_MPI=0 the runs on two MPI processes are skipped.

# The sparse matrix as the list of its rows, each a list of [column, value].
SparseRows:=function(entries, n_state)
    local rows, e;
    rows:=List([1..n_state], x->[]);
    for e in entries do
        Add(rows[e[1]], [e[2], e[3]]);
    od;
    return rows;
end;

SparseProduct:=function(A, B, n_state)
    local P, i, a, b, row, c;
    P:=[];
    for i in [1..n_state] do
        row:=[];
        for a in A[i] do
            for b in B[a[1]] do
                if IsBound(row[b[1]]) then
                    row[b[1]]:=row[b[1]] + a[2] * b[2];
                else
                    row[b[1]]:=a[2] * b[2];
                fi;
            od;
        od;
        P[i]:=[];
        for c in [1..Length(row)] do
            if IsBound(row[c]) and row[c] <> 0 then
                Add(P[i], [c, row[c]]);
            fi;
        od;
    od;
    return P;
end;

CheckAmfResult:=function(R)
    local weight, j, s, H, rows, sums, e, entry, i, k, A, B, dict;
    weight:=[];
    for j in [1..Length(R.ListAutOrder)] do
        for s in R.ListOrbitSize[j] do
            Add(weight, s / R.ListAutOrder[j]);
        od;
    od;
    if Length(weight) <> R.n_state then
        Print("The number of states is not the number of orbits\n");
        return false;
    fi;
    for H in R.ListHecke do
        sums:=ListWithIdenticalEntries(R.n_state, 0);
        dict:=NewDictionary([1, 1], true);
        for e in H.entries do
            sums[e[1]]:=sums[e[1]] + e[3];
            AddDictionary(dict, [e[1], e[2]], e[3]);
        od;
        if ForAny(sums, x->x <> H.n_neighbor) then
            Print("p=", H.prime, ": a row does not sum to the ", H.n_neighbor, " neighbors\n");
            return false;
        fi;
        for e in H.entries do
            entry:=LookupDictionary(dict, [e[2], e[1]]);
            if entry = fail then
                entry:=0;
            fi;
            if e[3] * weight[e[1]] <> entry * weight[e[2]] then
                Print("p=", H.prime, ": not self-adjoint at the states ", e[1], ", ", e[2], "\n");
                return false;
            fi;
        od;
    od;
    for i in [1..Length(R.ListHecke)] do
        for k in [i+1..Length(R.ListHecke)] do
            A:=SparseRows(R.ListHecke[i].entries, R.n_state);
            B:=SparseRows(R.ListHecke[k].entries, R.n_state);
            if SparseProduct(A, B, R.n_state) <> SparseProduct(B, A, R.n_state) then
                Print("The Hecke matrices at p=", R.ListHecke[i].prime, " and p=", R.ListHecke[k].prime, " do not commute\n");
                return false;
            fi;
        od;
    od;
    return true;
end;

FullTest:=function()
    local ListCase, use_mpi, eCase, ListGram, expected, R, R2;
    ListCase:=ReadAsFunction("cases.g")();
    use_mpi:=not (IsBound(GAPInfo.SystemEnvironment.AMF_MPI) and GAPInfo.SystemEnvironment.AMF_MPI = "0");
    for eCase in ListCase do
        Print("Case ", eCase.name, "\n");
        ListGram:=ReadListMatrixFile(eCase.classes);
        expected:=ReadAsFunction(Concatenation("expected_", eCase.name, ".g"))();
        R:=get_algebraic_modular_forms(ListGram, eCase.primes, eCase.level);
        if is_error(R) then
            return false;
        fi;
        if R <> expected then
            Print("The result of ", eCase.name, " is not the expected one\n");
            return false;
        fi;
        if CheckAmfResult(R) = false then
            Print("The checks of ", eCase.name, " failed\n");
            return false;
        fi;
        if eCase.mpi and use_mpi then
            R2:=get_algebraic_modular_forms(ListGram, eCase.primes, eCase.level, rec(np:=2));
            if is_error(R2) then
                return false;
            fi;
            if R2 <> R then
                Print("The result of ", eCase.name, " on two processes is not the same\n");
                return false;
            fi;
            Print("  the same on two MPI processes\n");
        fi;
        Print("  ", R.n_state, " states, correct\n");
    od;
    return true;
end;

result:=FullTest();
CI_Decision_Reset();
if result<>true then
    # Error case
    Print("Error case\n");
else
    # No error case
    Print("Normal case\n");
    CI_Write_Ok();
fi;
