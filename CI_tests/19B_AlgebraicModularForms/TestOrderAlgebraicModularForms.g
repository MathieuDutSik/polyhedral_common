Read("../common.g");
Print("Beginning of TestOrderAlgebraicModularForms\n");

# GENUS_AlgebraicModularForms over the ring of integers O_K of an imaginary
# quadratic field, on the standard unimodular Hermitian lattice I_m. The
# neighbors at the ramified prime alone must close the genus at the mass
# formula; with the split and inert primes as well, every row of a Hecke
# matrix sums to the number of neighbors, the matrices at a prime P with
# conj(P) = P (inert, ramified; not the split ones, written ell@r here,
# whose adjoint is T(conj P)) are self-adjoint for the weights 1 /
# |Aut(L_j)|, and two of them commute.

ListCase:=[
  rec(desc:="order_q7_3.txt", ramified:="7", primes:="7,2@0,3"),
  rec(desc:="order_q7_4.txt", ramified:="7", primes:="7,2@0,3"),
  rec(desc:="order_q11_3.txt", ramified:="11", primes:="11,3@0,2"),
  rec(desc:="order_q5_3.txt", ramified:="5", primes:="5,3@1,7@3")
];

RunOrderAmf:=function(desc, primes)
    local FileOut, prog, R;
    FileOut:=Filename(DirectoryTemporary(), "amf_order.g");
    prog:="../../src_genera/GENUS_AlgebraicModularForms";
    Exec(Concatenation(prog, " gmp ", desc, " 1 ", primes, " ", FileOut));
    if not IsExistingFile(FileOut) then
        Print("No output for ", desc, " at the primes ", primes, "\n");
        return fail;
    fi;
    R:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    return R;
end;

CheckHecke:=function(R)
    local n, H, M, i, j, k, A, B;
    n:=R.n_class;
    for H in R.ListHecke do
        M:=H.matrix;
        if ForAny(M, row->Sum(row) <> H.n_neighbor) then
            Print("p=", H.prime, ": a row does not sum to the ", H.n_neighbor, " neighbors\n");
            return false;
        fi;
        for i in [1..n] do
            for j in [1..n] do
                if not '@' in H.prime and M[i][j] / R.ListAutOrder[i] <> M[j][i] / R.ListAutOrder[j] then
                    Print("p=", H.prime, ": not self-adjoint at ", i, ", ", j, "\n");
                    return false;
                fi;
            od;
        od;
    od;
    for i in [1..Length(R.ListHecke)] do
        for k in [i+1..Length(R.ListHecke)] do
            A:=R.ListHecke[i].matrix;
            B:=R.ListHecke[k].matrix;
            if A * B <> B * A then
                Print("The Hecke matrices at ", R.ListHecke[i].prime, " and ", R.ListHecke[k].prime, " do not commute\n");
                return false;
            fi;
        od;
    od;
    return true;
end;

FullTest:=function()
    local eCase, R;
    for eCase in ListCase do
        Print("Case ", eCase.desc, "\n");
        R:=RunOrderAmf(eCase.desc, eCase.ramified);
        if R = fail or R.complete <> true then
            Print("The ramified prime ", eCase.ramified, " does not close the genus\n");
            return false;
        fi;
        Print("  ", R.n_class, " classes from the ramified prime ", eCase.ramified, "\n");
        R:=RunOrderAmf(eCase.desc, eCase.primes);
        if R = fail or R.complete <> true or CheckHecke(R) = false then
            Print("The checks at the primes ", eCase.primes, " failed\n");
            return false;
        fi;
        Print("  Hecke matrices at ", eCase.primes, " correct\n");
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
