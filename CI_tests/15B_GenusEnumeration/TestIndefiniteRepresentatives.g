Read("../common.g");
Read("../access_points.g");
Print("Beginning of TestIndefiniteRepresentatives\n");

# GENUS_IndefiniteRepresentatives on the cases of indefinite_cases.g: the
# number of classes of the genus is the expected one, there are as many
# representatives, the first is the given lattice, and the representatives
# are symmetric of the determinant of the given lattice (their genus is
# checked in the program with the SANITY_CHECK macros).

FullTest:=function()
    local ListCase, eCase, R, M;
    ListCase:=ReadAsFunction("indefinite_cases.g")();
    for eCase in ListCase do
        Print("Case ", eCase.name, "\n");
        R:=get_indefinite_representatives(eCase.gram);
        if is_error(R) then
            return false;
        fi;
        if R.n_class <> eCase.n_class or Length(R.ListGram) <> eCase.n_class then
            Print("The genus of ", eCase.name, " has ", R.n_class, " classes and ",
                  Length(R.ListGram), " representatives instead of ", eCase.n_class, "\n");
            return false;
        fi;
        if R.ListGram[1] <> eCase.gram then
            Print("The first representative of ", eCase.name, " is not the given lattice\n");
            return false;
        fi;
        for M in R.ListGram do
            if TransposedMat(M) <> M or DeterminantMat(M) <> DeterminantMat(eCase.gram) then
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
