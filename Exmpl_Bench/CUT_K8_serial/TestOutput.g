


LOrb:=ReadAsFunction("CutK8.orbitFacet")();


if Length(LOrb) <> 147 then
    GAP_EXIT_CODE(1);
else
    GAP_EXIT_CODE(0);
fi;
