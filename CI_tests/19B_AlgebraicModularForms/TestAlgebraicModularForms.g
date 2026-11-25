Read("../common.g");
Read("../access_points.g");
Print("Beginning of TestAlgebraicModularForms\n");

# The Hecke matrices of GENUS_MPI_AlgebraicModularForms on the cases of
# cases.g. For every case: the result is the expected one, every row of a
# Hecke matrix sums to the number of neighbors, the matrices are
# self-adjoint for the inner product of the states, |O| / |Aut(L_j)| for
# the state (L_j, O), and the matrices at two primes commute. With the
# environment variable AMF_MPI=0 the runs on two MPI processes are skipped.

Read("amf_checks.g");

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
