

GeneratePartialHypermetric:=function(n, n_orbit)
    local ListOrbit, GRP, FAC, eRepr, eOrb, ListIneq, FileSave, output;
    ListOrbit := GetOrbitInequalityHypermetricCone(n);
    Print("|ListOrbit|=", Length(ListOrbit), "\n");
    GRP:=SymmetricGroup(n);
    FAC:=[];
    for eRepr in ListOrbit{[1..n_orbit]}
    do
        eOrb:=Orbit(GRP, eRepr, Permuted);
        ListIneq:=List(eOrb, GetVectorFromHypermetric);
        Print("eRepr=", eRepr, " |ListIneq|=", Length(ListIneq), "\n");
        Append(FAC, ListIneq);
    od;
    FileSave:=Concatenation("HYPpartial_", String(n), "_", String(n_orbit));
    output:=OutputTextFile(FileSave, true);
    CPP_WriteMatrix(output, FAC);
    CloseStream(output);
end;



GeneratePartialHypermetric(5,2);
GeneratePartialHypermetric(6,2);
GeneratePartialHypermetric(7,2);
GeneratePartialHypermetric(8,2);
