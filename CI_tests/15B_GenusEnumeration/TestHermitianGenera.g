Read("../common.g");
Read("../access_points.g");
Print("Beginning of TestHermitianGenera\n");

# GENUS_HermitianGenera on the cases of hermitian_genera_cases.g: the genera
# of the rank and determinant are those listed by Hecke.jl, one for one up to
# the equivalence of the local symbols (the program fails otherwise), and a
# lattice of each of them is built.

FullTest:=function()
    local ListCase, eCase, R;
    ListCase:=ReadAsFunction("hermitian_genera_cases.g")();
    for eCase in ListCase do
        Print("Case ", eCase.name, "\n");
        R:=get_hermitian_genera(eCase);
        if is_error(R) then
            return false;
        fi;
        if Length(R) <> Length(eCase.genera) then
            Print("The case ", eCase.name, " has ", Length(R), " genera instead of ",
                  Length(eCase.genera), "\n");
            return false;
        fi;
        Print("  ", Length(R), " genera, correct\n");
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
