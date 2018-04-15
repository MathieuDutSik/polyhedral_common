U:=ReadVectorFile("EXMPL_ITER48");

eSHV:=U{[2..13]};
eRecTest2:=SHORT_TestRealizabilityShortestFamily_General(eSHV, "perfect", "glpk_secure");
#eRecTest2:=SHORT_TestRealizabilityShortestFamily_General(eSHV, "perfect", "cdd");
