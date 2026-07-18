Read("../common.g");
Read("../access_points.g");
Print("Beginning TestSerialDelaunayQueries\n");

# Unified exercise of LATT_SerialComputeDelaunay. Most categories below drive the
# binary with one &QUERIES option and validate the emitted output file:
#   FileDeformation    -> deformation of the quantizer along Q + t H
#   FileFreeVectors    -> free vectors of a lattice
#   FileHessian        -> Hessian of the normalized quantizer constant
#   FileRigidityDegree -> rigidity degree of a lattice
#   FileIsotropy       -> isotropy ("white" quantizer) test
# The last category instead exercises the GAP_Covering output format (covering
# density of a lattice), reached through the get_lattice_covering access point.
# The categories previously lived in separate CI entries (33_QuantizationDeformation,
# 34_FreeVectors, 35_QuantizationHessian, 35_RigidityDegree, 36_QuantizationIsotropy,
# 29_CoveringPerfectDim9); they are merged here since they all exercise
# LATT_SerialComputeDelaunay.

# Shared driver: writes the Gram matrix, assembles the namelist with the caller's
# &QUERIES body, and runs LATT_SerialComputeDelaunay. Returns nothing; the caller
# reads whatever output file its query was told to produce.
RunSerialDelaunay:=function(TmpDir, FileG, eG, queriesBody)
    local FileN, FileE, strOut, eProg, TheCommand;
    FileN:=Filename(TmpDir, "Query.nml");
    FileE:=Filename(TmpDir, "Query.err");
    WriteMatrixFile(FileG, eG);
    #
    strOut:="&SYSTEM\n";
    strOut:=Concatenation(strOut, " OutFormat = \"nothing\"\n");
    strOut:=Concatenation(strOut, " OutFile = \"unset.out\"\n");
    strOut:=Concatenation(strOut, " max_runtime_second = 0\n");
    strOut:=Concatenation(strOut, "/\n\n");
    strOut:=Concatenation(strOut, "&DATA\n");
    strOut:=Concatenation(strOut, " arithmetic = \"gmp\"\n");
    strOut:=Concatenation(strOut, " GRAMfile = \"", FileG, "\"\n");
    strOut:=Concatenation(strOut, " SVRfile = \"unset.svr\"\n");
    strOut:=Concatenation(strOut, " CacheFile = \"none\"\n");
    strOut:=Concatenation(strOut, "/\n\n");
    strOut:=Concatenation(strOut, "&QUERIES\n");
    strOut:=Concatenation(strOut, queriesBody);
    strOut:=Concatenation(strOut, "/\n");
    WriteStringFile(FileN, strOut);
    #
    eProg:=GetBinaryFilename("LATT_SerialComputeDelaunay");
    TheCommand:=Concatenation(eProg, " ", FileN, " 2> ", FileE);
    Exec(TheCommand);
    RemoveFile(FileN);
end;

# --------------------------------------------------------------------------- #
# Deformation of the quantizer along Q + t H (FileDeformation). The query value
# is the file holding the symmetric direction H; the result is written to that
# file with the suffix ".output". We test the E6 root lattice with H = e1 e1^T,
# for which the exact answer is known:
#   SecMoment(t) = (15/28 + 5/6 t + 1/9 t^2) / (1 + 4/3 t)
#   SecMoment(0)=15/28, SecMoment'(0)=5/42, SecMoment''(0)=-2/21
TestDeformation:=function(eRec)
    local TmpDir, FileQ, FileH, FileO, queriesBody, U, is_correct;
    TmpDir:=DirectoryTemporary();
    FileQ:=Filename(TmpDir, "Q.in");
    FileH:=Filename(TmpDir, "H.in");
    FileO:=Concatenation(FileH, ".output");
    WriteMatrixFile(FileH, eRec.H);
    queriesBody:=Concatenation(" FileDeformation = \"", FileH, "\"\n");
    RunSerialDelaunay(TmpDir, FileQ, eRec.Q, queriesBody);
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileQ);
    RemoveFile(FileH);
    RemoveFile(FileO);
    is_correct:=U.SecMoment0=eRec.SecMoment0
                and U.SecMoment1=eRec.SecMoment1
                and U.SecMoment2=eRec.SecMoment2
                and U.numerator=eRec.numerator
                and U.denominator=eRec.denominator;
    Print("name=", eRec.name, " is_correct=", is_correct, "\n");
    Print("  SecMoment0=", U.SecMoment0, " SecMoment1=", U.SecMoment1,
          " SecMoment2=", U.SecMoment2, "\n");
    return is_correct;
end;

ListRecDeformation:=[
rec(name:="E6_e1",
    Q:=[ [ 2, 1, 1, 0, 1, 1 ], [ 1, 4, 1, 1, 1, 3 ], [ 1, 1, 2, 1, 1, 1 ],
         [ 0, 1, 1, 2, 1, 2 ], [ 1, 1, 1, 1, 2, 2 ], [ 1, 3, 1, 2, 2, 4 ] ],
    H:=[ [ 1, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0 ],
         [ 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0 ], [ 0, 0, 0, 0, 0, 0 ] ],
    SecMoment0:=15/28, SecMoment1:=5/42, SecMoment2:=-2/21,
    numerator:=[ 15/28, 5/6, 1/9 ],
    denominator:=[ 1, 4/3, 0, 0, 0, 0, 0 ])
];

# --------------------------------------------------------------------------- #
# Free vectors of a lattice (FileFreeVectors). We test a few lattices for which
# the answer is known/stable.
TestFreeVectors:=function(eRec)
    local TmpDir, FileG, FileO, queriesBody, U, obtained, is_correct;
    TmpDir:=DirectoryTemporary();
    FileG:=Filename(TmpDir, "Gram.in");
    FileO:=Filename(TmpDir, "FreeVect.out");
    queriesBody:=Concatenation(" FileFreeVectors = \"", FileO, "\"\n");
    RunSerialDelaunay(TmpDir, FileG, eRec.eG, queriesBody);
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileG);
    RemoveFile(FileO);
    # Stable signature: [nbRelevantVector, sorted [SubspaceDim,OrbitSize,nMatched]].
    obtained:=rec(nbRelevantVector:=U.nbRelevantVector,
                  nbFreeVectorOrbit:=U.nbFreeVectorOrbit,
                  Signature:=Set(List(U.ListFreeVectorOrbit,
                      x->[x.SubspaceDim, x.OrbitSize, x.nMatched])));
    is_correct:=obtained.nbRelevantVector=eRec.nbRelevantVector
                and obtained.nbFreeVectorOrbit=eRec.nbFreeVectorOrbit
                and obtained.Signature=eRec.Signature;
    Print("name=", eRec.name, " is_correct=", is_correct, "\n");
    Print("  obtained=", obtained, "\n");
    return is_correct;
end;

ListRecFreeVectors:=[
rec(name:="A2", eG:=[ [ 2, 1 ], [ 1, 2 ] ],
    nbRelevantVector:=3, nbFreeVectorOrbit:=1,
    Signature:=Set([ [ 1, 3, 1 ] ])),
rec(name:="A4", eG:=[ [ 2, -1, 0, 0 ], [ -1, 2, -1, 0 ], [ 0, -1, 2, -1 ],
                      [ 0, 0, -1, 2 ] ],
    nbRelevantVector:=10, nbFreeVectorOrbit:=2,
    Signature:=Set([ [ 1, 5, 6 ], [ 1, 10, 4 ] ])),
rec(name:="E6",
    eG:=[ [ 2, 1, 1, 0, 1, 1 ], [ 1, 4, 1, 1, 1, 3 ], [ 1, 1, 2, 1, 1, 1 ],
          [ 0, 1, 1, 2, 1, 2 ], [ 1, 1, 1, 1, 2, 2 ], [ 1, 3, 1, 2, 2, 4 ] ],
    nbRelevantVector:=36, nbFreeVectorOrbit:=1,
    Signature:=Set([ [ 1, 27, 20 ] ]))
];

# --------------------------------------------------------------------------- #
# Hessian of the normalized quantizer constant G at a lattice Q (FileHessian).
# The method uses the derivative of the second-moment matrix M(Q) along rank-one
# directions v v^T (which span Sym^n by linearity of the moment derivative),
# assembled via the Aut(Q) orbit equivariance. We check the exact signature of
# the m=n(n+1)/2-1 dimensional Hessian, together with the two internal
# consistency residuals (radial: beta Q = 0 at a critical point; cross-check: an
# independent scalar deformation reproduces the Hessian prediction on a held-out
# direction).
TestHessian:=function(eRec)
    local TmpDir, FileG, FileO, queriesBody, U, obtained, is_correct;
    TmpDir:=DirectoryTemporary();
    FileG:=Filename(TmpDir, "Gram.in");
    FileO:=Filename(TmpDir, "Hessian.out");
    queriesBody:=Concatenation(" FileHessian = \"", FileO, "\"\n");
    RunSerialDelaunay(TmpDir, FileG, eRec.eG, queriesBody);
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileG);
    RemoveFile(FileO);
    obtained:=rec(nbPlus:=U.signature.nbPlus, nbMinus:=U.signature.nbMinus,
                  nbZero:=U.signature.nbZero);
    is_correct:=U.solved=true
                and U.radialResidual=0
                and U.checkResidual=0
                and obtained.nbPlus=eRec.nbPlus
                and obtained.nbMinus=eRec.nbMinus
                and obtained.nbZero=eRec.nbZero;
    Print("name=", eRec.name, " is_correct=", is_correct, "\n");
    Print("  signature=", obtained, " radialResidual=", U.radialResidual,
          " checkResidual=", U.checkResidual, " nbEval=", U.nbEval, "\n");
    return is_correct;
end;

ListRecHessian:=[
# A3 (FCC) is a degenerate critical point: 2 flat directions on the m=5 space.
rec(name:="A3", eG:=[ [ 2, -1, 0 ], [ -1, 2, -1 ], [ 0, -1, 2 ] ],
    nbPlus:=3, nbMinus:=0, nbZero:=2),
# E6 is a local minimum: positive definite Hessian on the m=20 space.
rec(name:="E6",
    eG:=[ [ 2, 1, 1, 0, 1, 1 ], [ 1, 4, 1, 1, 1, 3 ], [ 1, 1, 2, 1, 1, 1 ],
          [ 0, 1, 1, 2, 1, 2 ], [ 1, 1, 1, 1, 2, 2 ], [ 1, 3, 1, 2, 2, 4 ] ],
    nbPlus:=20, nbMinus:=0, nbZero:=0)
];

# --------------------------------------------------------------------------- #
# Rigidity degree of a lattice (FileRigidityDegree). Mirrors GAP's
# MyPolyhedral/lib/LatticeDelaunays.g::GetRigidityDegree.
TestRigidity:=function(eRec)
    local TmpDir, FileG, FileO, queriesBody, obtained, is_correct;
    TmpDir:=DirectoryTemporary();
    FileG:=Filename(TmpDir, "Gram.in");
    FileO:=Filename(TmpDir, "Rigid.out");
    queriesBody:=Concatenation(" FileRigidityDegree = \"", FileO, "\"\n");
    RunSerialDelaunay(TmpDir, FileG, eRec.eG, queriesBody);
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    obtained:=ReadAsFunction(FileO)();
    RemoveFile(FileG);
    RemoveFile(FileO);
    is_correct:=obtained=eRec.rigidity;
    Print("name=", eRec.name, " obtained=", obtained,
          " expected=", eRec.rigidity, " is_correct=", is_correct, "\n");
    return is_correct;
end;

ListRecRigidity:=[
rec(name:="Z2", eG:=[ [ 2, 0 ], [ 0, 2 ] ], rigidity:=2),
rec(name:="A2", eG:=[ [ 2, 1 ], [ 1, 2 ] ], rigidity:=3),
rec(name:="Z3", eG:=[ [ 2, 0, 0 ], [ 0, 2, 0 ], [ 0, 0, 2 ] ], rigidity:=3),
rec(name:="D4",
    eG:=[ [ 2, -1, 0, 0 ], [ -1, 2, -1, -1 ], [ 0, -1, 2, 0 ],
          [ 0, -1, 0, 2 ] ],
    rigidity:=1),
rec(name:="D5",
    eG:=[ [ 2, -1, 0, 0, 0 ], [ -1, 2, -1, 0, 0 ], [ 0, -1, 2, -1, -1 ],
          [ 0, 0, -1, 2, 0 ], [ 0, 0, -1, 0, 2 ] ],
    rigidity:=1),
rec(name:="E6",
    eG:=[ [ 2, -1, 0, 0, 0, 0 ], [ -1, 2, -1, 0, 0, 0 ],
          [ 0, -1, 2, -1, 0, -1 ], [ 0, 0, -1, 2, -1, 0 ],
          [ 0, 0, 0, -1, 2, 0 ], [ 0, 0, -1, 0, 0, 2 ] ],
    rigidity:=1)
];

# --------------------------------------------------------------------------- #
# Isotropy (extremal / "white" quantizer) test (FileIsotropy). A lattice is
# isotropic when the second-moment matrix M of its Voronoi cell is proportional
# to the inverse Gram matrix, equivalently GramMat * M = Lambda * I for a scalar
# Lambda. The test is exact: the defect DefectMat = GramMat*M - Lambda*I, with
# Lambda the average of the diagonal of GramMat*M, must be the zero matrix. We
# check the boolean and its internal consistency with the defect matrix.
TestIsotropy:=function(eRec)
    local TmpDir, FileG, FileO, queriesBody, U, n, defect_is_zero, is_correct;
    TmpDir:=DirectoryTemporary();
    FileG:=Filename(TmpDir, "Gram.in");
    FileO:=Filename(TmpDir, "Isotropy.out");
    queriesBody:=Concatenation(" FileIsotropy = \"", FileO, "\"\n");
    RunSerialDelaunay(TmpDir, FileG, eRec.eG, queriesBody);
    if IsExistingFile(FileO)=false then
        Print("The output file is not existing. That qualifies as a fail\n");
        return false;
    fi;
    U:=ReadAsFunction(FileO)();
    RemoveFile(FileG);
    RemoveFile(FileO);
    n:=Length(eRec.eG);
    # The boolean must agree with the defect being exactly zero, and with the
    # expected isotropy of the lattice.
    defect_is_zero:=U.DefectMat=NullMat(n, n);
    is_correct:=U.IsIsotropic=defect_is_zero
                and U.IsIsotropic=eRec.isIsotropic;
    Print("name=", eRec.name, " is_correct=", is_correct, "\n");
    Print("  IsIsotropic=", U.IsIsotropic, " Lambda=", U.Lambda, "\n");
    return is_correct;
end;

ListRecIsotropy:=[
# A3 (FCC / body-centered cubic dual) is an isotropic quantizer.
rec(name:="A3", eG:=[ [ 2, -1, 0 ], [ -1, 2, -1 ], [ 0, -1, 2 ] ],
    isIsotropic:=true),
# E6 is an isotropic quantizer.
rec(name:="E6",
    eG:=[ [ 2, 1, 1, 0, 1, 1 ], [ 1, 4, 1, 1, 1, 3 ], [ 1, 1, 2, 1, 1, 1 ],
          [ 0, 1, 1, 2, 1, 2 ], [ 1, 1, 1, 1, 2, 2 ], [ 1, 3, 1, 2, 2, 4 ] ],
    isIsotropic:=true),
# A rectangular lattice with unequal sides is not isotropic.
rec(name:="rect12", eG:=[ [ 1, 0 ], [ 0, 2 ] ],
    isIsotropic:=false)
];

# --------------------------------------------------------------------------- #
# Covering density of a lattice (GAP_Covering output format), reached through the
# get_lattice_covering access point. We run over a family of dimension-9 perfect
# forms and only check that each computation completes without error (the exact
# covering density is not pinned down here).
TestCovering:=function(eRec)
    local result;
    result:=get_lattice_covering(eRec.eMat);
    if is_error(result) then
        Print("name=", eRec.name, " is_correct=false (program error)\n");
        return false;
    fi;
    Print("name=", eRec.name, " TheCov=", result.TheCov, " det=", result.TheDet,
          " CovDensity=", result.CovDensity, "\n");
    return true;
end;

# The dimension-9 perfect forms are read from the PerfectMatrices directory.
ListRecCovering:=function()
    local eDir, ListFiles, L, eFile, eMat;
    eDir:="PerfectMatrices";
    ListFiles:=ListFileDirectory(eDir);
    L:=[];
    for eFile in ListFiles
    do
        eMat:=ReadMatrixFile(Concatenation(eDir, "/", eFile));
        Add(L, rec(name:=eFile, eMat:=eMat));
    od;
    return L;
end;

# --------------------------------------------------------------------------- #
AllCategories:=[
rec(label:="Deformation", tester:=TestDeformation, cases:=ListRecDeformation),
rec(label:="FreeVectors",  tester:=TestFreeVectors,  cases:=ListRecFreeVectors),
rec(label:="Hessian",      tester:=TestHessian,      cases:=ListRecHessian),
rec(label:="Rigidity",     tester:=TestRigidity,     cases:=ListRecRigidity),
rec(label:="Isotropy",     tester:=TestIsotropy,     cases:=ListRecIsotropy),
rec(label:="Covering",     tester:=TestCovering,     cases:=ListRecCovering())
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
