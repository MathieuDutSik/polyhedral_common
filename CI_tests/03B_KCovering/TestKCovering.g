Read("../common.g");
Read("../access_points.g");
Print("Beginning TestKCovering\n");

# Exercise of the order-k Delaunay programs of src_k_coverings:
#   LATT_SerialComputeKDelaunay          -> order-k tiling of one lattice
#   LATT_SerialLattice_IsoKDelaunayDomain -> enumeration of the (L,k)-types
# Phase 1 checks, per lattice and k, the number of orbits of tiles, the
# squared k-covering radius and the volume identity of the tiling.
# Phase 2 checks, per T-space and k, the number of (L,k)-types and the minimum
# over the types of the normalized k-covering density obtained by the
# determinant maximization.

# --------------------------------------------------------------------------- #
# Phase 1: the tiling of one lattice. The Gram matrix and the namelist are
# written in TmpDir; the output file, the volume file and the covering file
# are returned for the caller to read.
RunKDelaunay:=function(TmpDir, eG, k, OutFormat, FileVol)
    local FileG, FileN, FileO, FileE, strOut, eProg, TheCommand;
    FileG:=Filename(TmpDir, "Gram.in");
    FileN:=Filename(TmpDir, "KDel.nml");
    FileO:=Filename(TmpDir, "KDel.out");
    FileE:=Filename(TmpDir, "KDel.err");
    WriteMatrixFile(FileG, eG);
    #
    strOut:="&SYSTEM\n";
    strOut:=Concatenation(strOut, " OutFormat = \"", OutFormat, "\"\n");
    strOut:=Concatenation(strOut, " OutFile = \"", FileO, "\"\n");
    strOut:=Concatenation(strOut, " max_runtime_second = 0\n");
    strOut:=Concatenation(strOut, "/\n\n");
    strOut:=Concatenation(strOut, "&DATA\n");
    strOut:=Concatenation(strOut, " arithmetic = \"gmp\"\n");
    strOut:=Concatenation(strOut, " GRAMfile = \"", FileG, "\"\n");
    strOut:=Concatenation(strOut, " k = ", String(k), "\n");
    strOut:=Concatenation(strOut, " CacheFile = \"none\"\n");
    strOut:=Concatenation(strOut, "/\n\n");
    strOut:=Concatenation(strOut, "&QUERIES\n");
    strOut:=Concatenation(strOut, " FileTilingVolume = \"", FileVol, "\"\n");
    strOut:=Concatenation(strOut, "/\n");
    WriteStringFile(FileN, strOut);
    #
    eProg:=GetBinaryFilename("LATT_SerialComputeKDelaunay");
    TheCommand:=Concatenation(eProg, " ", FileN, " 2> ", FileE);
    Exec(TheCommand);
    RemoveFile(FileN);
    RemoveFile(FileG);
    return FileO;
end;

TestTiling:=function(eRec)
    local TmpDir, FileVol, FileO, U, V, W, is_correct;
    TmpDir:=DirectoryTemporary();
    FileVol:=Filename(TmpDir, "KDel.vol");
    # The tiles and the volume identity.
    FileO:=RunKDelaunay(TmpDir, eRec.eG, eRec.k, "GAP", FileVol);
    if IsExistingFile(FileO)=false or IsExistingFile(FileVol)=false then
        Print("The output files are not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileO)();
    V:=ReadAsFunction(FileVol)();
    RemoveFile(FileO);
    RemoveFile(FileVol);
    is_correct:=Length(U.ListTile)=eRec.n_tile and U.k=eRec.k and V.correct=true;
    if is_correct=false then
        Print("name=", eRec.name, " k=", eRec.k, "\n");
        Print("  |ListTile|=", Length(U.ListTile), " expected ", eRec.n_tile, "\n");
        Print("  volume check total=", V.total, " expected=", V.expected, " correct=", V.correct, "\n");
        return false;
    fi;
    # The squared k-covering radius.
    FileO:=RunKDelaunay(TmpDir, eRec.eG, eRec.k, "GAP_Covering", FileVol);
    if IsExistingFile(FileO)=false then
        Print("The covering output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    W:=ReadAsFunction(FileO)();
    RemoveFile(FileO);
    RemoveFile(FileVol);
    if W.TheCov<>eRec.cov_sq then
        Print("name=", eRec.name, " k=", eRec.k, " TheCov=", W.TheCov, " expected ", eRec.cov_sq, "\n");
        return false;
    fi;
    return true;
end;

# The stored values: the number of orbits of tiles under the affine
# isometries and the squared k-covering radius. For A2 with k = 2 the tiling
# is the kagome tiling (hexagons around the lattice points and triangles) and
# the 2-covering radius is the minimum norm.
ListRecTiling:=function()
    local eA2, eA3, eD4;
    eA2:=ClassicalSporadicLattices("A2");
    eA3:=ClassicalSporadicLattices("A3");
    eD4:=ClassicalSporadicLattices("D4");
    return [rec(name:="A2_k2", eG:=eA2, k:=2, n_tile:=2, cov_sq:=2),
            rec(name:="A2_k3", eG:=eA2, k:=3, n_tile:=1, cov_sq:=2),
            rec(name:="A2_k4", eG:=eA2, k:=4, n_tile:=2, cov_sq:=8/3),
            rec(name:="Z2_k2", eG:=IdentityMat(2), k:=2, n_tile:=2, cov_sq:=1),
            rec(name:="Z3_k2", eG:=IdentityMat(3), k:=2, n_tile:=2, cov_sq:=1),
            rec(name:="Z3_k3", eG:=IdentityMat(3), k:=3, n_tile:=3, cov_sq:=5/4),
            rec(name:="A3_k2", eG:=eA3, k:=2, n_tile:=3, cov_sq:=2),
            rec(name:="D4_k2", eG:=eD4, k:=2, n_tile:=2, cov_sq:=2)];
end;

# --------------------------------------------------------------------------- #
# Phase 2: the (L,k)-types of the space of all forms of dimension n.
RunIsoKDelaunay:=function(TmpDir, n, k, FileCov)
    local FileN, FileO, FileE, strOut, eProg, TheCommand;
    FileN:=Filename(TmpDir, "IsoKDel.nml");
    FileO:=Filename(TmpDir, "IsoKDel.out");
    FileE:=Filename(TmpDir, "IsoKDel.err");
    #
    strOut:="&SYSTEM\n";
    strOut:=Concatenation(strOut, " OutFormat = \"NumberGAP\"\n");
    strOut:=Concatenation(strOut, " OutFile = \"", FileO, "\"\n");
    strOut:=Concatenation(strOut, " max_runtime_second = 0\n");
    strOut:=Concatenation(strOut, "/\n\n");
    strOut:=Concatenation(strOut, "&DATA\n");
    strOut:=Concatenation(strOut, " arithmetic = \"gmp\"\n");
    strOut:=Concatenation(strOut, " k = ", String(k), "\n");
    strOut:=Concatenation(strOut, " FileCoveringOptimum = \"", FileCov, "\"\n");
    strOut:=Concatenation(strOut, "/\n\n");
    strOut:=Concatenation(strOut, "&TSPACE\n");
    strOut:=Concatenation(strOut, " TypeTspace = \"Classic\"\n");
    strOut:=Concatenation(strOut, " ClassicDim = ", String(n), "\n");
    strOut:=Concatenation(strOut, "/\n");
    WriteStringFile(FileN, strOut);
    #
    eProg:=GetBinaryFilename("LATT_SerialLattice_IsoKDelaunayDomain");
    TheCommand:=Concatenation(eProg, " ", FileN, " 2> ", FileE);
    Exec(TheCommand);
    RemoveFile(FileN);
    return FileO;
end;

TestEnumeration:=function(eRec)
    local TmpDir, FileCov, FileO, U, ListCov, eCov, best, tol;
    TmpDir:=DirectoryTemporary();
    FileCov:=Filename(TmpDir, "IsoKDel.cov");
    FileO:=RunIsoKDelaunay(TmpDir, eRec.n, eRec.k, FileCov);
    if IsExistingFile(FileO)=false or IsExistingFile(FileCov)=false then
        Print("The output files are not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileO)();
    ListCov:=ReadAsFunction(FileCov)();
    RemoveFile(FileO);
    RemoveFile(FileCov);
    if U.nb<>eRec.n_dom then
        Print("name=", eRec.name, " nb=", U.nb, " expected ", eRec.n_dom, "\n");
        return false;
    fi;
    # The minimum over the domains of the normalized k-covering density. The
    # optimization is numerical, hence the tolerance.
    best:=fail;
    for eCov in ListCov
    do
        if eCov.success=false then
            Print("name=", eRec.name, " one optimization did not converge: ", eCov.message, "\n");
            return false;
        fi;
        if best=fail or eCov.k_covering_density_normalized < best then
            best:=eCov.k_covering_density_normalized;
        fi;
    od;
    tol:=1.0e-6;
    if best=fail or AbsoluteValue(best - eRec.best_density) > tol then
        Print("name=", eRec.name, " best=", best, " expected ", eRec.best_density, "\n");
        return false;
    fi;
    return true;
end;

# The stored values. In the plane the optima are the thinnest k-fold lattice
# coverings of Blundon: 4 pi / sqrt(27) = 2 x (2 pi / sqrt(27)) for k = 2 and
# 25 pi / 18 (the square lattice) for k = 4. The dimension 3 case is the
# longest one (a few minutes).
ListRecEnumeration:=function()
    return [rec(name:="dim2_k2", n:=2, k:=2, n_dom:=1, best_density:=1.2091995762740997),
            rec(name:="dim2_k3", n:=2, k:=3, n_dom:=3, best_density:=1.1452093073591294),
            rec(name:="dim2_k4", n:=2, k:=4, n_dom:=9, best_density:=1.0908307829719477),
            rec(name:="dim3_k2", n:=3, k:=2, n_dom:=8, best_density:=1.3921133441560509)];
end;

# --------------------------------------------------------------------------- #
AllCategories:=[
rec(label:="Tiling",      tester:=TestTiling,      cases:=ListRecTiling()),
rec(label:="Enumeration", tester:=TestEnumeration, cases:=ListRecEnumeration())
];

FullTest:=function()
    local n_error, eCat, iRec, eRec;
    n_error:=0;
    for eCat in AllCategories
    do
        Print("=== Category ", eCat.label, " : ", Length(eCat.cases),
              " case(s) ===\n");
        iRec:=0;
        for eRec in eCat.cases
        do
            Print("iRec=", iRec, " / ", Length(eCat.cases), " name=", eRec.name,
                  "\n");
            if eCat.tester(eRec)=false then
                n_error:=n_error+1;
                return n_error;
            fi;
            iRec:=iRec + 1;
        od;
    od;
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
