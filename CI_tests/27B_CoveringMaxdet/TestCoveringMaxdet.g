Read("../common.g");
Print("Beginning Test covering density optimization over iso-Delaunay domains\n");

prog_enum := GetBinaryFilename("LATT_SerialLattice_IsoDelaunayDomain");
prog_per := GetBinaryFilename("LATT_SerialPeriodic_IsoDelaunayDomain");
prog_ana := GetBinaryFilename("LATT_AnalysisIsoDelaunay");
prog_rec := GetBinaryFilename("PERIODIC_LookForRecordCovering");

TmpDir := DirectoryTemporary();
tmp := function(name)
    return Filename(TmpDir, name);
end;

# The namelist enumerating the iso-Delaunay domains of the classic GL_dim(Z)
# T-space, dumping each of them as a boost archive DomPrefix<i>.
write_enum_nml := function(FileNml, dim, OutFile, DomPrefix)
    local os;
    RemoveFileIfExist(FileNml);
    os := OutputTextFile(FileNml, true);
    # GAP's line wrapping would break the long temp-dir paths in two.
    SetPrintFormattingStatus(os, false);
    AppendTo(os, "&SYSTEM\n");
    AppendTo(os, " max_runtime_second = 0\n");
    AppendTo(os, " ApplyStdUnitbuf = T\n");
    AppendTo(os, " Saving = F\n");
    AppendTo(os, " Prefix = \"/irrelevant/\"\n");
    AppendTo(os, " OutFile = \"", OutFile, "\"\n");
    AppendTo(os, " OutFormat = \"NumberGAP\"\n");
    AppendTo(os, "/\n");
    AppendTo(os, "&DATA\n");
    AppendTo(os, " arithmetic = \"gmp\"\n");
    AppendTo(os, " FileDualDescription = \"unset\"\n");
    AppendTo(os, " CommonGramMat = \"unset\"\n");
    AppendTo(os, " PrefixIsoDelaunayDomains = \"", DomPrefix, "\"\n");
    AppendTo(os, " CVPmethod = \"SVexact\"\n");
    AppendTo(os, "/\n");
    AppendTo(os, "&TSPACE\n");
    AppendTo(os, " TypeTspace = \"Classic\"\n");
    AppendTo(os, " ClassicDim = ", dim, "\n");
    AppendTo(os, "/\n");
    CloseStream(os);
end;

# The same for a periodic point set: the cosets are read from FileCosets, and
# the T-space actually used is written out for the analysis to consume.
write_per_nml := function(FileNml, dim, OutFile, DomPrefix, FileCosets,
                          FileLinSpa)
    local os;
    RemoveFileIfExist(FileNml);
    os := OutputTextFile(FileNml, true);
    SetPrintFormattingStatus(os, false);
    AppendTo(os, "&SYSTEM\n");
    AppendTo(os, " max_runtime_second = 0\n");
    AppendTo(os, " ApplyStdUnitbuf = T\n");
    AppendTo(os, " Saving = F\n");
    AppendTo(os, " Prefix = \"/irrelevant/\"\n");
    AppendTo(os, " OutFile = \"", OutFile, "\"\n");
    AppendTo(os, " OutFormat = \"NumberGAP\"\n");
    AppendTo(os, "/\n");
    AppendTo(os, "&DATA\n");
    AppendTo(os, " arithmetic = \"gmp\"\n");
    AppendTo(os, " FileDualDescription = \"unset\"\n");
    AppendTo(os, " FileCosets = \"", FileCosets, "\"\n");
    AppendTo(os, " PrefixIsoDelaunayDomains = \"", DomPrefix, "\"\n");
    AppendTo(os, " FileLinSpaceOut = \"", FileLinSpa, "\"\n");
    AppendTo(os, "/\n");
    AppendTo(os, "&TSPACE\n");
    AppendTo(os, " TypeTspace = \"Classic\"\n");
    AppendTo(os, " ClassicDim = ", dim, "\n");
    AppendTo(os, "/\n");
    CloseStream(os);
end;

# The namelist optimizing the covering density over one domain. FileLinSpa and
# FileCosets are "null" for a lattice domain.
write_ana_nml := function(FileNml, FileDom, OutFile, FileCov, FileLinSpa,
                          FileCosets)
    local os;
    RemoveFileIfExist(FileNml);
    os := OutputTextFile(FileNml, true);
    SetPrintFormattingStatus(os, false);
    AppendTo(os, "&SYSTEM\n");
    AppendTo(os, " OutFile = \"", OutFile, "\"\n");
    AppendTo(os, "/\n");
    AppendTo(os, "&DATA\n");
    AppendTo(os, " arithmetic = \"gmp\"\n");
    AppendTo(os, " FileIsoDelaunay = \"", FileDom, "\"\n");
    AppendTo(os, " FileLinSpace = \"", FileLinSpa, "\"\n");
    AppendTo(os, " FileCosets = \"", FileCosets, "\"\n");
    AppendTo(os, "/\n");
    AppendTo(os, "&QUERIES\n");
    AppendTo(os, " FileCoveringOptimum = \"", FileCov, "\"\n");
    AppendTo(os, "/\n");
    CloseStream(os);
end;

read_nb := function(FileOut)
    if IsExistingFile(FileOut) = false then
        return fail;
    fi;
    return ReadAsFunction(FileOut)().nb;
end;

# The covering optima of every domain DomPrefix<i>, i in [0, nb-1].
optimize_all := function(nb, DomPrefix, FileLinSpa, FileCosets)
    local ListRec, i, FileNml, FileOut, FileCov, eRec;
    ListRec := [];
    for i in [0 .. nb - 1] do
        FileNml := tmp(Concatenation("ana", String(i), ".nml"));
        FileOut := tmp(Concatenation("ana", String(i), ".out"));
        FileCov := tmp(Concatenation("cov", String(i), ".g"));
        RemoveFileIfExist(FileCov);
        write_ana_nml(FileNml, Concatenation(DomPrefix, String(i)), FileOut,
                      FileCov, FileLinSpa, FileCosets);
        Exec(Concatenation(prog_ana, " ", FileNml));
        if IsExistingFile(FileCov) = false then
            Print("  FOUND ERROR: no covering optimum written for domain ", i,
                  "\n");
            return fail;
        fi;
        eRec := ReadAsFunction(FileCov)();
        if eRec.success <> true then
            Print("  FOUND ERROR: domain ", i, " did not converge: ",
                  eRec.message, "\n");
            return fail;
        fi;
        Add(ListRec, eRec);
    od;
    return ListRec;
end;

# The optimum is normalized to covering radius 1, and the density is the
# quantity the whole computation is about, so both are checked.
check_radius_one := function(ListRec)
    local eRec;
    for eRec in ListRec do
        if AbsoluteValue(eRec.covering_radius_sq - 1.0) > 1.0e-6 then
            Print("  FOUND ERROR: covering_radius_sq=",
                  eRec.covering_radius_sq, " is not 1 at the optimum\n");
            return false;
        fi;
    od;
    return true;
end;

# ---------------------------------------------------------------------------
# The lattice case, against the classical optima.
#
# dim 3: the single domain gives A_3^*, of density 5 pi sqrt(5) / 24.
# dim 4: the best of the 3 domains is A_4^*.
# dim 5: the best of the 222 domains is A_5^*.
# ---------------------------------------------------------------------------

tol := 1.0e-7;

test_lattice_dim := function(eCase)
    local FileNml, FileOut, DomPrefix, nb, ListRec, densities, best;
    Print("Lattice dimension ", eCase.dim, "\n");
    FileNml := tmp("enum.nml");
    FileOut := tmp("enum.out");
    DomPrefix := Filename(TmpDir, Concatenation("dom", String(eCase.dim), "_"));
    RemoveFileIfExist(FileOut);
    write_enum_nml(FileNml, eCase.dim, FileOut, DomPrefix);
    Exec(Concatenation(prog_enum, " ", FileNml));
    nb := read_nb(FileOut);
    Print("  nb_domains=", nb, " (expected ", eCase.nb_domains, ")\n");
    if nb <> eCase.nb_domains then
        Print("  FOUND ERROR: wrong number of iso-Delaunay domains\n");
        return false;
    fi;
    ListRec := optimize_all(nb, DomPrefix, "null", "null");
    if ListRec = fail then
        return false;
    fi;
    if check_radius_one(ListRec) = false then
        return false;
    fi;
    densities := List(ListRec, x -> x.covering_density);
    best := Minimum(densities);
    Print("  best_density=", best, " (expected ", eCase.best_density, ")\n");
    if AbsoluteValue(best - eCase.best_density) > tol then
        Print("  FOUND ERROR: wrong optimal covering density\n");
        return false;
    fi;
    return true;
end;

# ---------------------------------------------------------------------------
# The periodic case: Z^3 + {0, (1/3,1/3,1/3)}.
#
# The value of the best periodic covering here is not a published constant, so
# what is checked is what is known independently: the covering radius is 1 at
# every optimum, the point density factor is 2 / 3^3, and no domain beats
# A_3^*, the best lattice covering of dimension 3 (conjecturally the best
# covering of dimension 3 altogether).
# ---------------------------------------------------------------------------

test_periodic := function()
    local FileNml, FileOut, FileCosets, FileLinSpa, DomPrefix, nb, ListRec,
          eRec, best;
    Print("Periodic point set Z^3 + {0, (1/3,1/3,1/3)}\n");
    FileCosets := tmp("cosets.txt");
    RemoveFileIfExist(FileCosets);
    WriteMatrixFile(FileCosets, [[0, 0, 0], [1/3, 1/3, 1/3]]);
    FileNml := tmp("per.nml");
    FileOut := tmp("per.out");
    FileLinSpa := tmp("linspa.txt");
    DomPrefix := Filename(TmpDir, "perdom_");
    RemoveFileIfExist(FileOut);
    write_per_nml(FileNml, 3, FileOut, DomPrefix, FileCosets, FileLinSpa);
    Exec(Concatenation(prog_per, " ", FileNml));
    nb := read_nb(FileOut);
    Print("  nb_domains=", nb, "\n");
    if nb = fail or nb < 1 then
        Print("  FOUND ERROR: the periodic enumeration produced no domain\n");
        return false;
    fi;
    ListRec := optimize_all(nb, DomPrefix, FileLinSpa, FileCosets);
    if ListRec = fail then
        return false;
    fi;
    if check_radius_one(ListRec) = false then
        return false;
    fi;
    for eRec in ListRec do
        if AbsoluteValue(eRec.point_density - 2.0 / 27.0) > 1.0e-12 then
            Print("  FOUND ERROR: point_density=", eRec.point_density,
                  " instead of 2/27\n");
            return false;
        fi;
    od;
    best := Minimum(List(ListRec, x -> x.covering_density));
    Print("  best_density=", best, " (A_3^* is 1.4635030689668)\n");
    if best < 1.4635030689668 - tol then
        Print("  FOUND ERROR: a periodic covering of dimension 3 beating ",
              "A_3^*, which contradicts the lattice optimum\n");
        return false;
    fi;
    return true;
end;

# ---------------------------------------------------------------------------
# The random-walk record search, on the point set whose answer the full
# enumeration above gives.
#
# The walk descends on the per-domain covering optimum and jumps out of local
# minima; over the 6 domains of Z^3 + {0, (1/3,1/3,1/3)} it has to recover
# the same 1.856151 the enumeration found, and it must not claim a record,
# A_3^* being out of reach for that point set.
#
# It also checks the record the program computes for itself: "auto" takes
# Theta(A_3^*) from its closed form, which has to be the value the
# enumeration of dimension 3 produced.
# ---------------------------------------------------------------------------

write_rec_nml := function(FileNml, dim, OutFile, FileCosets, budget)
    local os;
    RemoveFileIfExist(FileNml);
    os := OutputTextFile(FileNml, true);
    SetPrintFormattingStatus(os, false);
    AppendTo(os, "&SYSTEM\n");
    AppendTo(os, " max_runtime_second = ", budget, "\n");
    AppendTo(os, " ApplyStdUnitbuf = T\n");
    AppendTo(os, " OutFile = \"", OutFile, "\"\n");
    AppendTo(os, "/\n");
    AppendTo(os, "&DATA\n");
    AppendTo(os, " arithmetic = \"gmp\"\n");
    AppendTo(os, " FileDualDescription = \"unset\"\n");
    AppendTo(os, " FileCosets = \"", FileCosets, "\"\n");
    AppendTo(os, "/\n");
    AppendTo(os, "&SEARCH\n");
    AppendTo(os, " RecordToBeat = \"auto\"\n");
    AppendTo(os, " n_walk_steps = 5\n");
    AppendTo(os, " max_gram_entry = 1000000\n");
    AppendTo(os, "/\n");
    AppendTo(os, "&TSPACE\n");
    AppendTo(os, " TypeTspace = \"Classic\"\n");
    AppendTo(os, " ClassicDim = ", dim, "\n");
    AppendTo(os, "/\n");
    CloseStream(os);
end;

test_record_search := function(expected_best, expected_record)
    local FileCosets, FileNml, FileOut, eRec;
    Print("Record search on Z^3 + {0, (1/3,1/3,1/3)}\n");
    FileCosets := tmp("rec_cosets.txt");
    RemoveFileIfExist(FileCosets);
    WriteMatrixFile(FileCosets, [[0, 0, 0], [1/3, 1/3, 1/3]]);
    FileNml := tmp("rec.nml");
    FileOut := tmp("rec.g");
    RemoveFileIfExist(FileOut);
    write_rec_nml(FileNml, 3, FileOut, FileCosets, 30);
    Exec(Concatenation(prog_rec, " ", FileNml));
    if IsExistingFile(FileOut) = false then
        Print("  FOUND ERROR: the record search wrote no result\n");
        return false;
    fi;
    eRec := ReadAsFunction(FileOut)();
    Print("  message=", eRec.message, "\n");
    # "auto" has to reproduce the A_3^* value the enumeration computed.
    Print("  record=", eRec.record, " (expected ", expected_record, ")\n");
    if AbsoluteValue(eRec.record - expected_record) > 1.0e-9 then
        Print("  FOUND ERROR: the automatic record is not Theta(A_3^*)\n");
        return false;
    fi;
    # A periodic point set of dimension 3 cannot beat the lattice optimum,
    # so a claimed record here is a bug, not a discovery.
    if eRec.found_record <> false then
        Print("  FOUND ERROR: a record was claimed in dimension 3\n");
        return false;
    fi;
    if eRec.has_best <> true then
        Print("  FOUND ERROR: the search optimized no domain at all\n");
        return false;
    fi;
    # A walk that drifts into skewed representatives loses most of its
    # evaluations; the restart guard is what keeps that at zero.
    Print("  n_domain_evaluated=", eRec.n_domain_evaluated,
          " n_domain_failed=", eRec.n_domain_failed,
          " n_restart=", eRec.n_restart, "\n");
    if eRec.n_domain_failed > eRec.n_domain_evaluated then
        Print("  FOUND ERROR: most domains could not be optimized, the walk ",
              "has drifted into skewed representatives\n");
        return false;
    fi;
    # The walk has to reach the minimum over the 6 domains, which the full
    # enumeration gives independently.
    Print("  best_density=", eRec.best_density, " (expected ", expected_best,
          ")\n");
    if AbsoluteValue(eRec.best_density - expected_best) > 1.0e-6 then
        Print("  FOUND ERROR: the walk did not reach the known optimum\n");
        return false;
    fi;
    return true;
end;

ListCases := [
    rec(dim := 3, nb_domains := 1, best_density := 1.4635030689668180),
    rec(dim := 4, nb_domains := 3, best_density := 1.7655285081493524),
    rec(dim := 5, nb_domains := 222, best_density := 2.1242859089916246),
];

n_error := 0;
for eCase in ListCases do
    if test_lattice_dim(eCase) = false then
        n_error := n_error + 1;
    fi;
od;
if test_periodic() = false then
    n_error := n_error + 1;
fi;
if test_record_search(1.856151125516228, 1.4635030689668180) = false then
    n_error := n_error + 1;
fi;
Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
