Read("../common.g");

# The regular N-gon has vertices (cos(2*pi*k/N), sin(2*pi*k/N)), which for
# these N are real algebraic numbers, so the polytope lives over the field
# Q(x) with x = sin(2*pi/N). RegularNgons/FileDesc<N> describes that field
# and RegularNgons/Regular<N>gon holds the vertices over it. Both the
# automorphism group and the dual description therefore exercise the real
# algebraic arithmetic rather than the rational one.
#
# Two properties are checked, both known in advance:
#   * the automorphism group of the N-gon is dihedral of order 2N,
#   * the N-gon has exactly N facets, and every dual description method
#     must find that same number.

ListNgon:=[5, 7, 9];
ListMethod:=["bb", "cdd", "lrs", "normaliz"];

GetFieldFile:=function(n)
    return Concatenation("RealAlgebraic=RegularNgons/FileDesc", String(n));
end;

GetPolytopeFile:=function(n)
    return Concatenation("RegularNgons/Regular", String(n), "gon");
end;

# The automorphism group of the regular N-gon is the dihedral group of order
# 2N acting on the vertices: a rotation of order N and a reflection.
TestAutomorphismNgon:=function(n)
    local eProg, FileOut, TheCommand, GRP, ord;
    eProg:=GetBinaryFilename("GRP_LinPolytope_Automorphism");
    FileOut:=Filename(DirectoryTemporary(), "Ngon.grp");
    TheCommand:=Concatenation(eProg, " ", GetFieldFile(n), " ",
                              GetPolytopeFile(n), " GAP ", FileOut);
    Exec(TheCommand);
    if IsExistingFile(FileOut)=false then
        Print("N=", n, " the automorphism program created no output\n");
        return false;
    fi;
    GRP:=ReadAsFunction(FileOut)();
    RemoveFile(FileOut);
    ord:=Order(GRP);
    Print("N=", n, " |Aut|=", ord, " expected=", 2*n, "\n");
    if ord <> 2*n then
        Print("N=", n, " the automorphism group is not of order 2N\n");
        return false;
    fi;
    if LargestMovedPoint(GRP) > n then
        Print("N=", n, " the group moves a point beyond the vertices\n");
        return false;
    fi;
    if IsTransitive(GRP, [1..n])=false then
        Print("N=", n, " the group is not transitive on the vertices\n");
        return false;
    fi;
    # A dihedral group of order 2N has a cyclic subgroup of index 2.
    if IsCyclic(DerivedSubgroup(GRP))=false then
        Print("N=", n, " the group does not look dihedral\n");
        return false;
    fi;
    return true;
end;

# Read back the "|FAC|=<k>" that the Number output format writes.
ParseNbFacet:=function(FileName)
    local ListLines, eLine, eSplit;
    ListLines:=ReadTextFile(FileName);
    for eLine in ListLines
    do
        eSplit:=SplitString(eLine, "=");
        if Length(eSplit)=2 then
            return Int(eSplit[2]);
        fi;
    od;
    return fail;
end;

# The N-gon has N facets. Every method must agree on that.
TestDualDescNgon:=function(n)
    local eProg, eMethod, FileOut, TheCommand, nbFacet;
    eProg:=GetBinaryFilename("POLY_dual_description");
    for eMethod in ListMethod
    do
        FileOut:=Filename(DirectoryTemporary(), "Ngon.fac");
        RemoveFileIfExist(FileOut);
        TheCommand:=Concatenation(eProg, " ", GetFieldFile(n), " ", eMethod,
                                  " Number ", GetPolytopeFile(n), " ", FileOut);
        Exec(TheCommand);
        if IsExistingFile(FileOut)=false then
            Print("N=", n, " method=", eMethod, " created no output\n");
            return false;
        fi;
        nbFacet:=ParseNbFacet(FileOut);
        RemoveFile(FileOut);
        Print("N=", n, " method=", eMethod, " |FAC|=", nbFacet, " expected=", n, "\n");
        if nbFacet <> n then
            Print("N=", n, " method=", eMethod, " found ", nbFacet,
                  " facets instead of ", n, "\n");
            return false;
        fi;
    od;
    return true;
end;

#

CI_Decision_Reset();
n_error:=0;
for eN in ListNgon
do
    Print("---------------------------------------- N=", eN,
          " ----------------------------------------\n");
    if TestAutomorphismNgon(eN)=false then
        n_error:=n_error + 1;
    fi;
    if TestDualDescNgon(eN)=false then
        n_error:=n_error + 1;
    fi;
od;
Print("-------------------------------------------------------\n");
Print("n_error=", n_error, "\n");
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
