Read("../common.g");
Read("../access_points.g");
Read("amf_checks.g");
Print("Beginning of TestHermitianAlgebraicModularForms\n");

# The algebraic modular forms of genera of Hermitian lattices over CM
# fields: the classes by GENUS_HermitianEnumerate, the Hecke matrices by
# GENUS_MPI_AlgebraicModularForms (Ring = Hermitian) at every prime above
# the given rational primes, with or without a level structure, on the
# cases of hermitian_cases.g. For every case: the number of classes, the
# orders of the unitary groups, the number of states, and for each operator
# its prime ideal, its number of neighbors and its characteristic
# polynomial are the expected ones; the checks of amf_checks.g pass (at a
# split prime the adjoint is the operator at the conjugate prime). With
# mpi:=true the case is also run on two MPI processes (skipped with the
# environment variable AMF_MPI=0), which must give the same result.

DenseHecke:=function(R, H)
    local M, e;
    M:=NullMat(R.n_state, R.n_state);
    for e in H.entries do
        M[e[1]][e[2]]:=e[3];
    od;
    return M;
end;

FullTest:=function()
    local ListCase, use_mpi, eCase, R, R2, operators;
    ListCase:=ReadAsFunction("hermitian_cases.g")();
    use_mpi:=not (IsBound(GAPInfo.SystemEnvironment.AMF_MPI) and GAPInfo.SystemEnvironment.AMF_MPI = "0");
    for eCase in ListCase do
        Print("Case ", eCase.name, "\n");
        R:=get_hermitian_algebraic_modular_forms(eCase.field, eCase.symbol, eCase.primes, eCase.level, eCase.level_ideal);
        if is_error(R) then
            return false;
        fi;
        if R.n_class <> eCase.n_class or SortedList(R.ListAutOrder) <> eCase.auts or R.n_state <> eCase.n_state then
            Print("The classes or the states of ", eCase.name, " are not the expected ones: ", R.n_class, " classes, ",
                  SortedList(R.ListAutOrder), ", ", R.n_state, " states\n");
            return false;
        fi;
        operators:=List(R.ListHecke, H->[H.ideal, H.n_neighbor,
            CoefficientsOfUnivariatePolynomial(CharacteristicPolynomial(DenseHecke(R, H)))]);
        if operators <> eCase.operators then
            Print("The Hecke operators of ", eCase.name, " are not the expected ones: ", operators, "\n");
            return false;
        fi;
        if CheckAmfResult(R) = false then
            Print("The checks of ", eCase.name, " failed\n");
            return false;
        fi;
        if eCase.mpi and use_mpi then
            R2:=get_hermitian_algebraic_modular_forms(eCase.field, eCase.symbol, eCase.primes, eCase.level, eCase.level_ideal, rec(np:=2));
            if is_error(R2) then
                return false;
            fi;
            if R2 <> R then
                Print("The result of ", eCase.name, " on two processes is not the same\n");
                return false;
            fi;
            Print("  the same on two MPI processes\n");
        fi;
        Print("  ", R.n_class, " classes, ", R.n_state, " states, ", Length(R.ListHecke), " operators, correct\n");
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
