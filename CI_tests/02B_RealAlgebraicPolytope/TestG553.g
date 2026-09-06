Read("../common.g");

# G553 is a 7-dimensional cone on 150 generators over the real algebraic field
# Q(x) with x = 2*sin(2*pi/5), described by G553/FileDesc5. It carries a
# symmetry group (G553.grp) and is the real-size test of the dual description
# over a real algebraic field: the computation runs over the underlying ring
# Z[x] of that field, not over the field itself.
#
# Two independent computations of the same object are checked against each
# other:
#
#   * POLY_SerialDualDesc, the recursive dual description using the symmetry,
#     which returns the orbits of facets. Its output is compared with the
#     orbits_G553 reference shipped with the test: the number of orbits and
#     the sorted list of their incidences. The incidences are invariants of
#     the orbits, so the comparison does not depend on which representative
#     the algorithm happens to pick.
#
#   * POLY_dual_description, the direct dual description ignoring the
#     symmetry, run by each of its methods, which returns the number of
#     facets. It must agree with the total that the orbits account for,
#     computed here by expanding them under the group. The methods run the
#     computation over the underlying ring of the field, each with its own
#     kernel, so they are worth exercising separately.
#
# The two agreeing pins down the facet count without any hardcoded number:
# the reference file provides the orbits, the group file provides their sizes.

prefix:="G553";
FileOrbitRef:=Concatenation(prefix, "/orbits_G553");
FileOrbitOut:=Concatenation(prefix, "/orbits");

# The sorted incidences of a list of orbits. Two runs that find the same
# orbits agree on this, whichever representatives they return.
OrbitIncidences:=function(LOrb)
    return SortedList(List(LOrb, Length));
end;

# The recursive dual description, checked against the reference orbits.
TestSerialDualDesc:=function()
    local eProg, TheCommand, LOrbRef, LOrbOut;
    LOrbRef:=ReadAsFunction(FileOrbitRef)();
    eProg:=GetBinaryFilename("POLY_SerialDualDesc");
    RemoveFileIfExist(FileOrbitOut);
    TheCommand:=Concatenation("(cd ", prefix, " && ", eProg, " input_G553.nml)");
    Exec(TheCommand);
    if IsExistingFile(FileOrbitOut)=false then
        Print("POLY_SerialDualDesc created no output\n");
        return false;
    fi;
    LOrbOut:=ReadAsFunction(FileOrbitOut)();
    RemoveFile(FileOrbitOut);
    Print("|orbits|=", Length(LOrbOut), " expected=", Length(LOrbRef), "\n");
    if Length(LOrbOut)<>Length(LOrbRef) then
        Print("The number of orbits of facets is not the expected one\n");
        return false;
    fi;
    if OrbitIncidences(LOrbOut)<>OrbitIncidences(LOrbRef) then
        Print("The incidences of the orbits are not the expected ones\n");
        Print("obtained=", OrbitIncidences(LOrbOut), "\n");
        Print("expected=", OrbitIncidences(LOrbRef), "\n");
        return false;
    fi;
    return true;
end;

# The number of facets that the reference orbits account for: the orbits are
# expanded under the group, and their union must have that many elements,
# which is the statement that the orbits are pairwise distinct.
GetNbFacetFromOrbits:=function()
    local GRP, LOrb, ListOrbit, nb_sum, TheUnion;
    GRP:=ReadGroupFile(Concatenation(prefix, "/G553.grp"));
    LOrb:=ReadAsFunction(FileOrbitRef)();
    ListOrbit:=List(LOrb, x->Orbit(GRP, Set(x), OnSets));
    nb_sum:=Sum(List(ListOrbit, Length));
    TheUnion:=Union(ListOrbit);
    if Length(TheUnion)<>nb_sum then
        Print("The reference orbits are not pairwise distinct: |union|=",
              Length(TheUnion), " sum=", nb_sum, "\n");
        return fail;
    fi;
    return nb_sum;
end;

# The direct dual description, without the symmetry, checked against that
# total. lrs is left out: it is correct here but takes minutes on this cone,
# where the other methods take about a second.
ListMethod:=["bb", "cdd", "normaliz"];

TestDirectDualDesc:=function(nb_facet_expected)
    local eProg, eMethod, FileOut, TheCommand, nbFacet, result;
    eProg:=GetBinaryFilename("POLY_dual_description");
    result:=true;
    for eMethod in ListMethod
    do
        FileOut:=Filename(DirectoryTemporary(), "G553.fac");
        RemoveFileIfExist(FileOut);
        TheCommand:=Concatenation("(cd ", prefix, " && ", eProg,
                                  " RealAlgebraic=FileDesc5 ", eMethod,
                                  " Number G553.ext ", FileOut, ")");
        Exec(TheCommand);
        if IsExistingFile(FileOut)=false then
            Print("method=", eMethod, " created no output\n");
            result:=false;
            continue;
        fi;
        nbFacet:=ParseNbFacet(FileOut);
        RemoveFile(FileOut);
        Print("method=", eMethod, " |FAC|=", nbFacet, " expected=",
              nb_facet_expected, "\n");
        if nbFacet<>nb_facet_expected then
            Print("method=", eMethod,
                  " disagrees with the orbits\n");
            result:=false;
        fi;
    od;
    return result;
end;

#

CI_Decision_Reset();
n_error:=0;
Print("---------------------------------------- G553 ",
      "----------------------------------------\n");
if TestSerialDualDesc()=false then
    n_error:=n_error + 1;
fi;
nb_facet:=GetNbFacetFromOrbits();
if nb_facet=fail then
    n_error:=n_error + 1;
else
    Print("The ", Length(ReadAsFunction(FileOrbitRef)()),
          " orbits account for ", nb_facet, " facets\n");
    if TestDirectDualDesc(nb_facet)=false then
        n_error:=n_error + 1;
    fi;
fi;
Print("-------------------------------------------------------\n");
Print("n_error=", n_error, "\n");
if n_error > 0 then
    Print("Error case\n");
else
    Print("Normal case\n");
    CI_Write_Ok();
fi;
