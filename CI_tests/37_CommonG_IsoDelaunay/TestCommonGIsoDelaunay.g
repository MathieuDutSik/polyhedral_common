Read("../common.g");
Print("Beginning Test CommonGramMat iso-Delaunay (rigid lattices and stars)\n");

prog_enum := GetBinaryFilename("LATT_SerialLattice_IsoDelaunayDomain");
prog_ana := GetBinaryFilename("LATT_AnalysisIsoDelaunay");
prog_canon := GetBinaryFilename("LATT_Canonicalize");

TmpDir := DirectoryTemporary();
tmp := function(name)
    return Filename(TmpDir, name);
end;


# Write the namelist that enumerates the iso-Delaunay domains of the classic
# GL_dim(Z) T-space. When DomPrefix <> "unset" every domain is dumped as a boost
# archive DomPrefix<i>. When CommonG <> "unset" the enumeration is restricted to
# the star of the Gram matrix stored in that file.
write_enum_nml := function(FileNml, dim, OutFile, DomPrefix, CommonG)
    local os;
    RemoveFileIfExist(FileNml);
    os := OutputTextFile(FileNml, true);
    # Disable GAP's automatic line wrapping: the temp-dir paths are long and a
    # break in the middle of a "..." field corrupts the namelist.
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
    AppendTo(os, " CommonGramMat = \"", CommonG, "\"\n");
    AppendTo(os, " PrefixIsoDelaunayDomains = \"", DomPrefix, "\"\n");
    AppendTo(os, " CVPmethod = \"SVexact\"\n");
    AppendTo(os, "/\n");
    AppendTo(os, "&TSPACE\n");
    AppendTo(os, " TypeTspace = \"Classic\"\n");
    AppendTo(os, " ClassicDim = ", dim, "\n");
    AppendTo(os, "/\n");
    CloseStream(os);
end;

# Read the "return rec(nb:=K);" number written by the enumeration.
read_nb := function(FileOut)
    if IsExistingFile(FileOut) = false then
        return fail;
    fi;
    return ReadAsFunction(FileOut)().nb;
end;

# Enumerate the iso-Delaunay domains of dimension dim, dumping each of them.
# Returns rec(nb, DomPrefix).
enumerate_and_dump := function(dim)
    local FileNml, FileOut, DomPrefix, nb;
    FileNml := tmp("enum.nml");
    FileOut := tmp("enum.out");
    DomPrefix := Filename(TmpDir, "dom_");
    RemoveFileIfExist(FileOut);
    write_enum_nml(FileNml, dim, FileOut, DomPrefix, "unset");
    Exec(Concatenation(prog_enum, " ", FileNml));
    nb := read_nb(FileOut);
    return rec(nb := nb, DomPrefix := DomPrefix);
end;

# Canonical form (up to GL_dim(Z)) of a Gram matrix, used as the dedup key.
canonical_form := function(GramMat)
    local FileI, FileO;
    FileI := tmp("canon.gram");
    FileO := tmp("canon.out");
    RemoveFileIfExist(FileI);
    RemoveFileIfExist(FileO);
    WriteMatrixFile(FileI, GramMat);
    Exec(Concatenation(prog_canon, " gmp ", FileI, " GAP ", FileO));
    return ReadAsFunction(FileO)();
end;

# For every dumped domain, extract its full-rank extreme rays and reduce them up
# to GL_dim(Z)-equivalence. Returns the list of distinct canonical rigid forms
# sorted by determinant.
extract_rigid := function(DomPrefix, nbDomains)
    local distinct, i, FileNml, FileRays, DomFile, rays, M, C, os;
    distinct := [];
    for i in [0..nbDomains-1] do
        DomFile := Concatenation(DomPrefix, String(i));
        FileNml := tmp("ana.nml");
        FileRays := tmp("rays.g");
        RemoveFileIfExist(FileNml);
        RemoveFileIfExist(FileRays);
        os := OutputTextFile(FileNml, true);
        SetPrintFormattingStatus(os, false);
        AppendTo(os, "&DATA\n");
        AppendTo(os, " arithmetic = \"gmp\"\n");
        AppendTo(os, " FileIsoDelaunay = \"", DomFile, "\"\n");
        AppendTo(os, "/\n");
        AppendTo(os, "&SYSTEM\n");
        AppendTo(os, " OutFile = \"stderr\"\n");
        AppendTo(os, "/\n");
        AppendTo(os, "&QUERIES\n");
        AppendTo(os, " FileFullRankRays = \"", FileRays, "\"\n");
        AppendTo(os, "/\n");
        CloseStream(os);
        Exec(Concatenation(prog_ana, " ", FileNml));
        rays := ReadAsFunction(FileRays)();
        for M in rays do
            C := canonical_form(M);
            if Position(distinct, C) = fail then
                Add(distinct, C);
            fi;
        od;
    od;
    SortBy(distinct, DeterminantMat);
    return distinct;
end;

# Count the iso-Delaunay domains of the star of the Gram matrix GramMat.
compute_star := function(dim, GramMat)
    local FileNml, FileOut, FileG;
    FileG := tmp("common.gram");
    FileNml := tmp("star.nml");
    FileOut := tmp("star.out");
    RemoveFileIfExist(FileG);
    RemoveFileIfExist(FileOut);
    WriteMatrixFile(FileG, GramMat);
    write_enum_nml(FileNml, dim, FileOut, "unset", FileG);
    Exec(Concatenation(prog_enum, " ", FileNml));
    return read_nb(FileOut);
end;

# Run the full pipeline for one dimension and compare against the reference.
TestDim := function(eRec)
    local res, rigid, dets, stars;
    Print("=== dimension ", eRec.dim, " ===\n");
    res := enumerate_and_dump(eRec.dim);
    Print("  nb_domains=", res.nb, " (expected ", eRec.nb_domains, ")\n");
    if res.nb <> eRec.nb_domains then
        Print("  FOUND ERROR: wrong number of iso-Delaunay domains\n");
        return false;
    fi;
    rigid := extract_rigid(res.DomPrefix, res.nb);
    dets := List(rigid, DeterminantMat);
    Print("  n_rigid=", Length(rigid), " (expected ", eRec.n_rigid,
          ") dets=", dets, " (expected ", eRec.dets, ")\n");
    if Length(rigid) <> eRec.n_rigid then
        Print("  FOUND ERROR: wrong number of rigid lattices\n");
        return false;
    fi;
    if AsSortedList(dets) <> AsSortedList(eRec.dets) then
        Print("  FOUND ERROR: wrong rigid-lattice determinants\n");
        return false;
    fi;
    stars := List(rigid, x -> compute_star(eRec.dim, x));
    Print("  star_counts=", stars, " (expected ", eRec.star_counts, ")\n");
    if AsSortedList(stars) <> AsSortedList(eRec.star_counts) then
        Print("  FOUND ERROR: wrong star counts\n");
        return false;
    fi;
    return true;
end;

# Reference values.
# dim 4: 3 iso-Delaunay domains, 1 rigid lattice (D4, det 4), star 2.
# dim 5: 222 iso-Delaunay domains, 7 rigid lattices, star counts as below.
ListCases := [
    rec(dim := 4, nb_domains := 3, n_rigid := 1, dets := [4],
        star_counts := [2]),
    rec(dim := 5, nb_domains := 222, n_rigid := 7,
        dets := [4, 16, 48, 96, 432, 648, 1024],
        star_counts := [258, 172, 21, 237, 87, 24, 48]),
];

n_error := 0;
for eCase in ListCases do
    if TestDim(eCase) = false then
        n_error := n_error + 1;
    fi;
od;
Print("n_error=", n_error, "\n");

CI_Decision_Reset();
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
