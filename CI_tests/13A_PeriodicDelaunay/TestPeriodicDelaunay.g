Read("../common.g");
Print("Beginning TestPeriodicDelaunay\n");

# Exercise of the periodic point set support of src_delaunay, through the
# two programs that expose it to the user, LATT_SerialPeriodicDelaunay and
# LATT_SerialPeriodic_IsoDelaunayDomain.
#
# Two kinds of check are made.
#
# --- Against values verified outside the C++ code. The orbit data of the
# A4 and HE6 families, and the single orbit of periodic iso-Delaunay
# domains of Z^2 + {(0,0), (1/3,0)}, were computed independently by the
# GAP Periodic_DelaunayDescriptionStandard and Periodic_EnumerationProcedureLtype
# and agreed exactly (August 2026).
#
# --- Against nothing at all, by invariance. A unimodular change of basis
# U sends the point set D Z^n + {c_i} with the form G to the point set
# with the form U^-1 G U^-T and the cosets c_i U taken modulo one, which
# is the same point set written in another basis: the orbit data has to
# come out identical, adjacency counts included. That check needs no
# reference value, and it is the one that exercises the machinery which
# used to have its own TEST binaries: the closest vector computation and
# the subgroup of the transformations preserving the point set are redone
# in the new basis, and the equivalence between domains has to match the
# orbits up again.

# GAP's Exec drops the exit status; recover it through the shell.
Exec_GetReturnValue:=function(TheCommand)
    local TmpDir, FileRC, list_lines;
    TmpDir:=DirectoryTemporary();
    FileRC:=Filename(TmpDir, "rc");
    Exec(Concatenation("( ", TheCommand, " ) ; echo $? > ", FileRC));
    list_lines:=ReadTextFile(FileRC);
    if Length(list_lines)=0 then
        return 1;
    fi;
    return Int(list_lines[1]);
end;

# The fractional part, GAP's Int truncating towards zero
FracPart:=function(x)
    local y;
    y:=x - Int(x);
    if y < 0 then
        y:=y + 1;
    fi;
    return y;
end;

# The point set in the basis given by the unimodular U
TransformCosets:=function(ListCoset, U)
    return List(ListCoset, c->List(c * U, FracPart));
end;

TransformGram:=function(GramMat, U)
    local V;
    V:=Inverse(U);
    return V * GramMat * TransposedMat(V);
end;

# A unimodular matrix of size n that is not a permutation: the identity
# with a one added above the diagonal.
ShearMatrix:=function(n)
    local U;
    U:=IdentityMat(n);
    U[1][2]:=1;
    return U;
end;

# One run of LATT_SerialPeriodic_IsoDelaunayDomain. Returns the record it
# wrote, or fail.
RunIsoDelaunay:=function(ListCoset, dim, OutFormat, tag)
    local eProg, TmpDir, FileCosets, FileNml, FileRes, FileErr, strNml, TheCommand, TheResult;
    eProg:=GetBinaryFilename("LATT_SerialPeriodic_IsoDelaunayDomain");
    TmpDir:=DirectoryTemporary();
    FileCosets:=Filename(TmpDir, Concatenation("Cosets_", tag));
    FileNml:=Filename(TmpDir, Concatenation("Iso_", tag, ".nml"));
    FileRes:=Filename(TmpDir, Concatenation("IsoRes_", tag));
    FileErr:=Filename(TmpDir, Concatenation("iso_", tag, ".err"));
    WriteMatrixFile(FileCosets, ListCoset);
    strNml:="&SYSTEM\n";
    strNml:=Concatenation(strNml, " max_runtime_second = 0\n");
    strNml:=Concatenation(strNml, " ApplyStdUnitbuf = T\n");
    strNml:=Concatenation(strNml, " Saving = F\n");
    strNml:=Concatenation(strNml, " Prefix = \"/irrelevant/\"\n");
    strNml:=Concatenation(strNml, " OutFile = \"", FileRes, "\"\n");
    strNml:=Concatenation(strNml, " OutFormat = \"", OutFormat, "\"\n");
    strNml:=Concatenation(strNml, "/\n\n");
    strNml:=Concatenation(strNml, "&DATA\n");
    strNml:=Concatenation(strNml, " arithmetic = \"gmp\"\n");
    strNml:=Concatenation(strNml, " FileDualDescription = \"unset\"\n");
    strNml:=Concatenation(strNml, " FileCosets = \"", FileCosets, "\"\n");
    strNml:=Concatenation(strNml, "/\n\n");
    strNml:=Concatenation(strNml, "&TSPACE\n");
    strNml:=Concatenation(strNml, " TypeTspace = \"Classic\"\n");
    strNml:=Concatenation(strNml, " ClassicDim = ", String(dim), "\n");
    strNml:=Concatenation(strNml, "/\n");
    WriteStringFile(FileNml, strNml);
    TheCommand:=Concatenation(eProg, " ", FileNml, " 2> ", FileErr);
    TheResult:=Exec_GetReturnValue(TheCommand);
    if TheResult<>0 then
        Print("LATT_SerialPeriodic_IsoDelaunayDomain failed on ", tag, ", its output:\n");
        Exec(Concatenation("cat ", FileErr));
        return fail;
    fi;
    return ReadAsFunction(FileRes)();
end;

# One run of LATT_SerialPeriodicDelaunay. Returns the record it wrote, or
# fail.
RunPeriodicDelaunay:=function(arg)
    local GramMat, ListCoset, tag, expect_failure, eProg, TmpDir, FileG, FileC, FileNml, FileRes, FileErr, strNml, TheCommand, TheResult;
    GramMat:=arg[1];
    ListCoset:=arg[2];
    tag:=arg[3];
    # A run that is meant to fail should not print as if something went
    # wrong, or the log reads as a failure while the test passes
    expect_failure:=false;
    if Length(arg) >= 4 then
        expect_failure:=arg[4];
    fi;
    eProg:=GetBinaryFilename("LATT_SerialPeriodicDelaunay");
    TmpDir:=DirectoryTemporary();
    FileG:=Filename(TmpDir, Concatenation("Gram_", tag));
    FileC:=Filename(TmpDir, Concatenation("Cosets_", tag));
    FileNml:=Filename(TmpDir, Concatenation("Run_", tag, ".nml"));
    FileRes:=Filename(TmpDir, Concatenation("Result_", tag));
    FileErr:=Filename(TmpDir, Concatenation("run_", tag, ".err"));
    WriteMatrixFile(FileG, GramMat);
    WriteMatrixFile(FileC, ListCoset);
    strNml:="&SYSTEM\n";
    strNml:=Concatenation(strNml, " max_runtime_second = 0\n");
    strNml:=Concatenation(strNml, " ApplyStdUnitbuf = T\n");
    strNml:=Concatenation(strNml, " Saving = F\n");
    strNml:=Concatenation(strNml, " Prefix = \"/irrelevant/\"\n");
    strNml:=Concatenation(strNml, " OutFile = \"", FileRes, "\"\n");
    strNml:=Concatenation(strNml, " OutFormat = \"SummaryGAP\"\n");
    strNml:=Concatenation(strNml, "/\n\n&DATA\n");
    strNml:=Concatenation(strNml, " arithmetic = \"gmp\"\n");
    strNml:=Concatenation(strNml, " GRAMfile = \"", FileG, "\"\n");
    strNml:=Concatenation(strNml, " FileCosets = \"", FileC, "\"\n");
    strNml:=Concatenation(strNml, "/\n");
    WriteStringFile(FileNml, strNml);
    TheCommand:=Concatenation(eProg, " ", FileNml, " 2> ", FileErr);
    TheResult:=Exec_GetReturnValue(TheCommand);
    if TheResult<>0 then
        if not expect_failure then
            Print("LATT_SerialPeriodicDelaunay failed on ", tag, "\n");
            Exec(Concatenation("cat ", FileErr));
        fi;
        return fail;
    fi;
    return ReadAsFunction(FileRes)();
end;

# The invariants of an iso-Delaunay run, in a shape that does not depend
# on the order the orbits came out in
IsoSummary:=function(TheRec)
    return SortedList(List(TheRec.ListEntry,
                           x->[x.GRPpermSize, x.n_ineq_red, x.det, x.n_shv]));
end;

DelaunaySummary:=function(TheRec)
    return SortedList(List(TheRec.ListRec, x->[x.nVert, x.ordStab, x.nAdj]));
end;

# The periodic iso-Delaunay domains of a point set in the T-space of all
# the symmetric matrices, and the same in a sheared basis.
RunIsoDelaunayCases:=function()
    local ListCase, n_error, eCase, TheRec, TheRecShear, U, ListCosetShear, summ;
    ListCase:=[
      rec(name:="Z2_third",   cosets:=[[0,0],[1/3,0]], dim:=2, nb:=1,
          summary:=[[1,3,2,4]]),
      rec(name:="Z2_quarter", cosets:=[[0,0],[1/4,0]], dim:=2, nb:=3,
          summary:=[[1,3,2,4],[2,2,1,4],[2,2,2,4]]),
      rec(name:="Z2_fifth",   cosets:=[[0,0],[1/5,0]], dim:=2, nb:=4,
          summary:=[[1,3,1,4],[1,3,2,4],[1,3,20,4],[1,3,35,4]])
    ];
    n_error:=0;
    for eCase in ListCase
    do
        TheRec:=RunIsoDelaunay(eCase.cosets, eCase.dim, "DetailedObjectGAP", eCase.name);
        if TheRec=fail then
            n_error:=n_error+1;
            continue;
        fi;
        summ:=IsoSummary(TheRec);
        Print(eCase.name, ": n_obj=", TheRec.n_obj, " summary=", summ, "\n");
        if TheRec.n_obj<>eCase.nb then
            Print("The number of orbits should be ", eCase.nb, "\n");
            n_error:=n_error+1;
        fi;
        if summ<>eCase.summary then
            Print("The orbit invariants differ from the expected ones\n");
            n_error:=n_error+1;
        fi;
        # The same point set in a sheared basis has to give the same thing
        U:=ShearMatrix(eCase.dim);
        ListCosetShear:=TransformCosets(eCase.cosets, U);
        TheRecShear:=RunIsoDelaunay(ListCosetShear, eCase.dim, "DetailedObjectGAP",
                                    Concatenation(eCase.name, "_shear"));
        if TheRecShear=fail then
            n_error:=n_error+1;
            continue;
        fi;
        if TheRecShear.n_obj<>TheRec.n_obj or IsoSummary(TheRecShear)<>summ then
            Print("The sheared basis of ", eCase.name, " gives ",
                  TheRecShear.n_obj, " orbits and ", IsoSummary(TheRecShear),
                  " instead of the same as the original\n");
            n_error:=n_error+1;
        else
            Print("  the sheared basis gives the same orbits\n");
        fi;
    od;
    return n_error;
end;

# The A4 + cosets family (c generating A4*/A4 = Z/5 in the A4 basis) and
# HE6 = E6 + {0, c6} (c6 of order 3 in E6*/E6), and the same in a sheared
# basis.
RunDelaunayCases:=function()
    local GramA4, c, c2, GramE6, c6, ListCase, n_error, eCase, TheRec, TheRecShear, ListSiz, U, GramShear, CosetShear;
    GramA4:=[[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,2]];
    c :=[4/5, 3/5, 2/5, 1/5];
    c2:=[3/5, 1/5, 4/5, 2/5];
    GramE6:=[[2,-1,0,0,0,0],[-1,2,-1,0,0,0],[0,-1,2,-1,0,-1],
             [0,0,-1,2,-1,0],[0,0,0,-1,2,0],[0,0,-1,0,0,2]];
    c6:=[1/3, 2/3, 0, 1/3, 2/3, 0];
    ListCase:=[
      rec(name:="A4+c",   gram:=GramA4, cosets:=[[0,0,0,0], c],     expected:=[ [5,120], [20,240] ]),
      rec(name:="A4+2c",  gram:=GramA4, cosets:=[[0,0,0,0], c2],    expected:=[ [5,120], [8,48], [10,240] ]),
      rec(name:="A4+c2c", gram:=GramA4, cosets:=[[0,0,0,0], c, c2], expected:=[ [5,8], [5,120], [8,48] ]),
      rec(name:="HE6",    gram:=GramE6, cosets:=[[0,0,0,0,0,0], c6], expected:=[ [54,103680] ])
    ];
    n_error:=0;
    for eCase in ListCase
    do
        TheRec:=RunPeriodicDelaunay(eCase.gram, eCase.cosets, eCase.name);
        if TheRec=fail then
            n_error:=n_error+1;
            continue;
        fi;
        ListSiz:=Set(List(TheRec.ListRec, x->[x.nVert, x.ordStab]));
        Print(eCase.name, ": nb=", TheRec.nb, " orbits=", ListSiz, "\n");
        if ListSiz<>Set(eCase.expected) then
            Print("The orbit data differs from the GAP-verified value\n");
            n_error:=n_error+1;
        fi;
        # The same point set in a sheared basis has to give the same thing
        U:=ShearMatrix(Length(eCase.gram));
        GramShear:=TransformGram(eCase.gram, U);
        CosetShear:=TransformCosets(eCase.cosets, U);
        TheRecShear:=RunPeriodicDelaunay(GramShear, CosetShear,
                                         Concatenation(eCase.name, "_shear"));
        if TheRecShear=fail then
            n_error:=n_error+1;
            continue;
        fi;
        if DelaunaySummary(TheRecShear)<>DelaunaySummary(TheRec) then
            Print("The sheared basis of ", eCase.name, " gives ",
                  DelaunaySummary(TheRecShear), " instead of ",
                  DelaunaySummary(TheRec), "\n");
            n_error:=n_error+1;
        else
            Print("  the sheared basis gives the same orbits, adjacencies included\n");
        fi;
    od;
    return n_error;
end;

# A set of cosets that is a group modulo Z^n describes a lattice, which
# the periodic programs refuse, pointing at LATT_SerialComputeDelaunay.
RunRejectLattice:=function()
    local GramMat, TheRec;
    GramMat:=[[2,-1],[-1,2]];
    Print("The cosets {0, 1/2} of Z, which form a group, have to be refused\n");
    TheRec:=RunPeriodicDelaunay(GramMat, [[0,0],[1/2,0],[0,1/2],[1/2,1/2]],
                                "reject", true);
    if TheRec<>fail then
        Print("The point set is a lattice and should have been refused\n");
        return 1;
    fi;
    Print("  refused, as it should be\n");
    return 0;
end;

FullTest:=function()
    local n_error;
    n_error:=0;
    n_error:=n_error + RunIsoDelaunayCases();
    n_error:=n_error + RunDelaunayCases();
    n_error:=n_error + RunRejectLattice();
    Print("FullTest: n_error=", n_error, "\n");
    return n_error;
end;

NestFunction:=function()
    local n_error;
    n_error:=FullTest();
    CI_Decision_Reset();
    if n_error > 0 then
        Print("Error case\n");
    else
        Print("Normal case\n");
        CI_Write_Ok();
    fi;
    CI_PrintExistConclusion();
end;

NestFunction();
