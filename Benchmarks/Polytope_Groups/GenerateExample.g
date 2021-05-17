n:=8;

GRA:=CompleteGraph(Group(()), n);
eRec:=CMC_GetCutPolytope(GRA);

eRecCUT:=eRec.GetCUT_info();

FileOut:=Concatenation("CUT_", String(n));
output:=OutputTextFile(FileOut, true);
CPP_WriteMatrix(output, eRecCUT.EXT);
CloseStream(output);


