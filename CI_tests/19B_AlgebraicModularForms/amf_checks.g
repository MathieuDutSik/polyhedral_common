# The checks of a result of GENUS_MPI_AlgebraicModularForms: every row of a
# Hecke matrix sums to the number of neighbors; T(s, t) w(s) = T'(t, s) w(t),
# w(s) = |O| / |Aut(L_j)| the weight of the state s = (L_j, O) and T' the
# adjoint of T (T itself over Z, the operator at the conjugate prime over a
# CM field, given by the field adjoint); the matrices at two primes commute.

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
    local weight, j, s, H, H_adj, rows, sums, e, entry, i, k, A, B, dict;
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
        for e in H.entries do
            sums[e[1]]:=sums[e[1]] + e[3];
        od;
        # the adjoint: the operator itself, or the one at the conjugate prime
        dict:=NewDictionary([1, 1], true);
        if IsBound(H.adjoint) then
            H_adj:=R.ListHecke[H.adjoint];
        else
            H_adj:=H;
        fi;
        for e in H_adj.entries do
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
                Print("p=", H.prime, ": not adjoint to its adjoint at the states ", e[1], ", ", e[2], "\n");
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
