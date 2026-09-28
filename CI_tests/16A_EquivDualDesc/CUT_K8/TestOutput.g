


LOrb:=ReadAsFunction("CutK8.orbitFacet")();


if Length(LOrb) <> 147 then
    Print("Error case\n");
    GAP_EXIT_CODE(1);
else
    Print("Normal case\n");
    GAP_EXIT_CODE(0);
fi;
