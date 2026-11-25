Read("../common.g");
Read("../access_points.g");
Print("Beginning of TestHermitianEnumeration\n");

# The genera of Hermitian lattices of hermitian_enumeration_cases.g, by
# GENUS_HermitianEnumerate from their symbols: the enumeration must be
# complete, with the mass, the class number and the orders of the unitary
# groups of Hecke.jl. Over an abelian field the program computes the mass
# exactly; over a field given by a polynomial it is given, or, with
# compute_mass, computed by the program from numerical L-values (certified
# with Arb, the CI being built with flint).
FullTest:=function()
    local ListCase, eCase, mass_in, R;
    ListCase:=ReadAsFunction("hermitian_enumeration_cases.g")();
    for eCase in ListCase do
        Print("Case ", eCase.name, "\n");
        if IsBound(eCase.field.conductor) or IsBound(eCase.compute_mass) then
            mass_in:=fail;
        else
            mass_in:=eCase.mass;
        fi;
        R:=get_hermitian_enumeration(eCase.field, eCase.symbol, mass_in);
        if is_error(R) then
            return false;
        fi;
        if R.complete <> true then
            Print("The enumeration of ", eCase.name, " is not complete\n");
            return false;
        fi;
        if R.mass <> eCase.mass then
            Print("The mass of ", eCase.name, " is ", R.mass, " instead of ", eCase.mass, "\n");
            return false;
        fi;
        if IsBound(eCase.auts) and SortedList(R.ListAutOrder) <> SortedList(eCase.auts) then
            Print("The unitary group orders of ", eCase.name, " are ", SortedList(R.ListAutOrder), "\n");
            return false;
        fi;
        if IsBound(eCase.n_class) and Length(R.ListBasis) <> eCase.n_class then
            Print("The class number of ", eCase.name, " is ", Length(R.ListBasis), " instead of ", eCase.n_class, "\n");
            return false;
        fi;
        Print("  ", Length(R.ListBasis), " classes, correct\n");
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
