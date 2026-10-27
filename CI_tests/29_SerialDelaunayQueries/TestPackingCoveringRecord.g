Read("../common.g");
Print("Beginning TestPackingCoveringRecord\n");

# Packing-covering constant of the dimension-5 record, through
# LATT_SerialPeriodicDelaunay.
#
# The record is a 2-periodic NON-LATTICE point set X = Z^5 + {0, c} whose
# packing-covering constant gamma = mu / rho is strictly below
#
#     gamma_5 = sqrt(3/2 + sqrt(13)/6) = 1.4494568681327953...,
#
# the best value over LATTICES, attained by Horvath's Ho_5. So the
# dimension-5 packing-covering problem is not solved by a lattice.
#
# What is checked here is the rational witness: an explicit point set with
# Gram matrix and coset of denominator 10^5. Being rational it is handled
# exactly, and being a legitimate point set in its own right the inequality
# it satisfies is a complete result -- no perturbation argument is needed to
# get back to the optimum. See PackingCoveringRecord/README.md for the
# algebraic optimum itself, which is a different matter.

# Runs LATT_SerialPeriodicDelaunay on the witness and returns the record of
# its PackingCoveringGAP output, or fail when there is no output.
RunRecord:=function()
    local TmpDir, FileN, FileO, FileE, strOut, eProg, TheCommand, U;
    TmpDir:=DirectoryTemporary();
    FileN:=Filename(TmpDir, "record.nml");
    FileO:=Filename(TmpDir, "record.out");
    FileE:=Filename(TmpDir, "record.err");
    #
    strOut:="&SYSTEM\n";
    strOut:=Concatenation(strOut, " OutFormat = \"PackingCoveringGAP\"\n");
    strOut:=Concatenation(strOut, " OutFile = \"", FileO, "\"\n");
    strOut:=Concatenation(strOut, " max_runtime_second = 0\n");
    strOut:=Concatenation(strOut, "/\n\n");
    strOut:=Concatenation(strOut, "&DATA\n");
    strOut:=Concatenation(strOut, " arithmetic = \"gmp\"\n");
    strOut:=Concatenation(strOut, " GRAMfile = \"PackingCoveringRecord/Gram_witness\"\n");
    strOut:=Concatenation(strOut, " FileCosets = \"PackingCoveringRecord/Cosets_witness\"\n");
    strOut:=Concatenation(strOut, "/\n");
    WriteStringFile(FileN, strOut);
    #
    eProg:=GetBinaryFilename("LATT_SerialPeriodicDelaunay");
    TheCommand:=Concatenation(eProg, " ", FileN, " 2> ", FileE);
    Exec(TheCommand);
    RemoveFile(FileN);
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return fail;
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileO);
    return U;
end;

# gamma^2 is an exact rational, so the comparison with gamma_5 is exact:
#   gamma < gamma_5  <=>  6 gamma^2 - 9 < sqrt(13)  <=>  (6 gamma^2 - 9)^2 < 13,
# both sides being positive.
TestRecord:=function()
    local U, g2, mu2, lam2, lhs, is_correct;
    U:=RunRecord();
    if U=fail then
        return false;
    fi;
    g2:=U.SquarePackingCoveringConstant;
    mu2:=U.SquareCoveringRadius;
    lam2:=U.SquareMinimalDistance;
    lhs:=6*g2 - 9;
    Print("orbits of Delaunay cells  : ", U.nb, "\n");
    Print("mu^2                      : ", Float(mu2), "\n");
    Print("lambda^2                  : ", Float(lam2), "\n");
    Print("gamma^2 = 4 mu^2/lambda^2 : ", Float(g2), "\n");
    Print("gamma                     : ", Sqrt(Float(g2)), "\n");
    Print("gamma_5                   : 1.4494568681328\n");
    Print("6 gamma^2 - 9             : ", Float(lhs), "\n");
    Print("(6 gamma^2 - 9)^2         : ", Float(lhs^2), " vs 13\n");
    is_correct:=g2=4*mu2/lam2 and lhs > 0 and lhs^2 < 13;
    if is_correct then
        Print("gamma < gamma_5, proved in exact rational arithmetic\n");
    else
        Print("gamma < gamma_5 is NOT proved\n");
    fi;
    return is_correct;
end;

test:=TestRecord();
CI_Decision_Reset();
if test then
    CI_Write_Ok();
fi;
