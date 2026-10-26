Read("../common.g");
Print("Beginning TestPolyNorm\n");

# Runs the program with the arguments (a string) and the standard output
# written to FileOut, the standard error being discarded. Returns true when
# the exit status is 0.
RunProgram:=function(eProg, TheArg, FileOut)
    local FileRC, list_lines;
    FileRC:=Filename(DirectoryTemporary(), "rc");
    Exec(Concatenation("( ", GetBinaryFilename(eProg), " ", TheArg, " > ", FileOut,
                       " 2> /dev/null ) ; echo $? > ", FileRC));
    list_lines:=ReadTextFile(FileRC);
    if Length(list_lines)=0 or list_lines[1]<>"0" then
        Print("Failure of ", eProg, " ", TheArg, "\n");
        return false;
    fi;
    return true;
end;

# Runs a program with the GAP output on stdout and returns the record, or
# fail when the run fails.
RunProgramGAP:=function(eProg, TheArg)
    local FileOut, answer;
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    if RunProgram(eProg, TheArg, FileOut)=false then
        return fail;
    fi;
    answer:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    return answer;
end;

# Runs a program with the text output on stdout and returns that text, or
# fail when the run fails.
RunProgramText:=function(eProg, TheArg)
    local FileOut, text;
    FileOut:=Filename(DirectoryTemporary(), "Test.out");
    if RunProgram(eProg, TheArg, FileOut)=false then
        return fail;
    fi;
    text:=StringFile(FileOut);
    RemoveFile(FileOut);
    return text;
end;

# Checks that the text contains each of the expected lines.
CheckLines:=function(name, text, list_line)
    local eLine, test;
    test:=true;
    for eLine in list_line
    do
        if PositionSublist(text, eLine)=fail then
            Print(name, ": the line \"", eLine, "\" is missing\n");
            test:=false;
        fi;
    od;
    if test=false then
        Print(text, "\n");
    fi;
    return test;
end;

# Packing and covering scalars of polytopes with respect to Z^n, against
# their known values. When brute_force is true, the covering radius is also
# computed by the brute force enumeration of Cslovjecsek, Malikiosis,
# Naszodi and Schymura.
TestPackingCovering:=function(eRec)
    local res_pack, res_cov, res_bf, test;
    Print("Packing and covering of ", eRec.file, " alpha=", eRec.alpha, " mu=", eRec.mu, "\n");
    test:=true;
    res_pack:=RunProgramGAP("POLYNORM_Packing", Concatenation("gmp ", eRec.file, " GAP stdout"));
    if res_pack=fail or res_pack.alpha<>eRec.alpha then
        Print("Failure of the packing of ", eRec.file, "\n");
        test:=false;
    fi;
    res_cov:=RunProgramGAP("POLYNORM_Covering", Concatenation("gmp ", eRec.file, " GAP stdout"));
    if res_cov=fail or res_cov.mu<>eRec.mu then
        Print("Failure of the covering of ", eRec.file, "\n");
        test:=false;
    fi;
    if eRec.brute_force then
        res_bf:=RunProgramGAP("POLYNORM_Covering", Concatenation("gmp ", eRec.file, " GAP stdout brute_force"));
        if res_bf=fail or res_bf.mu<>eRec.mu then
            Print("Failure of the brute force covering of ", eRec.file, "\n");
            test:=false;
        fi;
    fi;
    return test;
end;

# The brute force is run on the cases of dimension 2 and on the cube (4
# million systems).
ListPackingCovering:=[
  rec(file:="square.ext", alpha:=1, mu:=1, brute_force:=true),
  rec(file:="diamond.ext", alpha:=1/2, mu:=1, brute_force:=true),
  rec(file:="triangle.ext", alpha:=1, mu:=2, brute_force:=true),
  rec(file:="hexagon.ext", alpha:=1/2, mu:=2/3, brute_force:=true),
  rec(file:="pentagon.ext", alpha:=1/4, mu:=1/2, brute_force:=true),
  rec(file:="asym2d.ext", alpha:=1/2, mu:=5/6, brute_force:=true),
  rec(file:="cube3.ext", alpha:=1, mu:=1, brute_force:=true),
  rec(file:="octa3.ext", alpha:=1/2, mu:=3/2, brute_force:=false),
  rec(file:="simplex3.ext", alpha:=1, mu:=3, brute_force:=false),
  rec(file:="asym3d.ext", alpha:=1/2, mu:=4/3, brute_force:=false),
  rec(file:="frac3d.ext", alpha:=2/3, mu:=12/5, brute_force:=false),
  rec(file:="cube4.ext", alpha:=1, mu:=1, brute_force:=false),
  rec(file:="cross4.ext", alpha:=1/2, mu:=2, brute_force:=false),
  rec(file:="simplex4.ext", alpha:=1, mu:=4, brute_force:=false),
  rec(file:="cell24.ext", alpha:=1/2, mu:=1, brute_force:=false)
];

# The densest lattice packings of the tetrahedron (18/49, Hoylman) and of
# the octahedron (18/19, Minkowski), as given by Betke and Henk, and the
# thinnest lattice covering of the octahedron (9/8 with the body-centered
# cubic lattice, Lian and Xue) are evaluated directly, and must be
# certified strict local optima.
TestCertify:=function(eRec)
    local text;
    Print("Certification of ", eRec.file, " with ", eRec.lat, "\n");
    text:=RunProgramText("POLYNORM_LatticeCertify", Concatenation("gmp ", eRec.file, " ", eRec.lat, " stdout"));
    if text=fail then
        return false;
    fi;
    return CheckLines(eRec.file, text, eRec.lines);
end;

ListCertify:=[
  rec(file:="tetra_betke_henk.ext", lat:="tetra_betke_henk.lat",
      lines:=["Packing: alpha=1 density=18/49",
              "certified_local_optimum=1 (strict local optimum on every branch)"]),
  rec(file:="octa3.ext", lat:="octa3_betke_henk.lat",
      lines:=["Packing: alpha=1 density=18/19",
              "certified_local_optimum=1 (strict local optimum on every branch)"]),
  rec(file:="octa3.ext", lat:="octa3_bcc.lat",
      lines:=["Covering: mu=3/2 density=9/8",
              "certified_local_optimum=1 (strict local optimum)"])
];

# The search from a few random starts (the generator is portable, so the
# seed fixes the run) must find the known optimum exactly, polished and
# certified. Only the part of the output from "Best lattice" on is checked.
TestOptimalLattice:=function(eRec)
    local TheArg, text, pos;
    Print("Optimal lattice of ", eRec.file, " for ", eRec.problem, " density=", eRec.density, "\n");
    TheArg:=Concatenation("gmp ", eRec.file, " ", eRec.problem, " ", String(eRec.nstart), " 1 none stdout");
    text:=RunProgramText("POLYNORM_OptimalLattice", TheArg);
    if text=fail then
        return false;
    fi;
    pos:=PositionSublist(text, "Best lattice");
    if pos=fail then
        Print("No best lattice found for ", eRec.file, "\n");
        return false;
    fi;
    return CheckLines(eRec.file, text{[pos..Length(text)]},
                      [Concatenation("exact=", String(eRec.density)), "certified_local_optimum=1"]);
end;

# Dimension 2: Fary (triangle), Sriamorn (quadrilateral). Dimension 3:
# Betke and Henk.
ListOptimalLattice:=[
  rec(file:="triangle.ext", problem:="packing", nstart:=2, density:=2/3),
  rec(file:="asym2d.ext", problem:="packing", nstart:=2, density:=3/4),
  rec(file:="triangle.ext", problem:="covering", nstart:=2, density:=3/2),
  rec(file:="asym2d.ext", problem:="covering", nstart:=2, density:=4/3),
  rec(file:="simplex3.ext", problem:="packing", nstart:=3, density:=18/49),
  rec(file:="octa3.ext", problem:="packing", nstart:=3, density:=18/19)
];

# The algorithm of Betke and Henk: densest lattice packings of 3-polytopes,
# over the rationals and the fields Q(sqrt 2), Q(sqrt 5), against the
# values of their paper (12 digits).
TestBetkeHenk:=function(eRec)
    local text;
    Print("Betke-Henk of ", eRec.file, " (", eRec.arith, ") density ", eRec.density, "\n");
    text:=RunProgramText("POLYNORM_BetkeHenk", Concatenation(eRec.arith, " ", eRec.file, " stdout 4"));
    if text=fail then
        return false;
    fi;
    return CheckLines(eRec.file, text,
                      [Concatenation("Densest lattice packing density ~ ", eRec.density),
                       "Every candidate of smaller determinant is resolved"]);
end;

# Octahedron 18/19, truncated cube 9/(5 + 3 sqrt 2), dodecahedron (2 + tau)/4.
ListBetkeHenk:=[
  rec(arith:="gmp", file:="octa3.ext", density:="0.947368421053"),
  rec(arith:="Qsqrt2", file:="truncated_cube_qsqrt2.ext", density:="0.973747687989"),
  rec(arith:="Qsqrt5", file:="dodecahedron_qsqrt5.ext", density:="0.904508497187")
];

# The maximal kissing number: the number of touching translates of the
# densest-contact lattice packing, with no undecided configuration left.
TestMaxKissing:=function(eRec)
    local res;
    Print("Maximal kissing number of ", eRec.file, " kiss=", eRec.kiss, "\n");
    res:=RunProgramGAP("POLYNORM_MaxKissing", Concatenation("gmp ", eRec.file, " GAP stdout"));
    if res=fail then
        return false;
    fi;
    if res.kiss<>eRec.kiss or res.n_undecided<>0 then
        Print("kiss=", res.kiss, " n_undecided=", res.n_undecided, "\n");
        return false;
    fi;
    return true;
end;

# The square and the cube are parallelotopes and reach 3^n - 1; the
# triangle and the hexagon reach 6, the most a plane convex body other
# than a parallelogram can.
ListMaxKissing:=[
  rec(file:="square.ext", kiss:=8),
  rec(file:="triangle.ext", kiss:=6),
  rec(file:="hexagon.ext", kiss:=6),
  rec(file:="cube3.ext", kiss:=26)
];

# All the cases are run, so that every failure is reported.
ListResult:=Concatenation(List(ListPackingCovering, TestPackingCovering),
                          List(ListCertify, TestCertify),
                          List(ListOptimalLattice, TestOptimalLattice),
                          List(ListBetkeHenk, TestBetkeHenk),
                          List(ListMaxKissing, TestMaxKissing));
test:=ForAll(ListResult, x->x);
Print("test=", test, "\n");
CI_Decision_Reset();
if test then
    CI_Write_Ok();
fi;
