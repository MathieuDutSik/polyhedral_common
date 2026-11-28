Read("../common.g");
Read("../access_points.g");
Print("Beginning of TestHermitianIndefinite\n");

# GENUS_HermitianIndefinite on the cases of hermitian_indefinite_cases.g:
# the number of classes of the genus is that of Hecke.jl, there are as many
# representatives, of the determinant of the given lattice (their local
# symbols at the primes of the neighbors are checked in the program with the
# SANITY_CHECK macros).

FullTest:=function()
    local ListCase, eCase, R, M, det;
    ListCase:=ReadAsFunction("hermitian_indefinite_cases.g")();
    for eCase in ListCase do
        Print("Case ", eCase.name, "\n");
        R:=get_hermitian_indefinite(eCase.d, eCase.h0, eCase.h1);
        if is_error(R) then
            return false;
        fi;
        if R.n_class <> eCase.n_class or Length(R.ListGram) <> eCase.n_class then
            Print("The genus of ", eCase.name, " has ", R.n_class, " classes and ",
                  Length(R.ListGram), " representatives instead of ", eCase.n_class, "\n");
            return false;
        fi;
        det:=DeterminantMat(R.ListGram[1]);
        for M in R.ListGram do
            if TransposedMat(M) <> M or DeterminantMat(M) <> det then
                Print("A representative of ", eCase.name, " is not symmetric of the determinant\n");
                return false;
            fi;
        od;
        Print("  ", R.n_class, " classes, correct\n");
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
