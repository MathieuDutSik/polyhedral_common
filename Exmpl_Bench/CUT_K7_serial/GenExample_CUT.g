n:=7;

GRA:=CompleteGraph(Group(()), n);
eRec:=CMC_GetCutPolytope(GRA);

eRecCUT:=eRec.GetCUT_info();

FileOutExt:=Concatenation("CUT_", String(n), ".ext");
RemoveFileIfExist(FileOutExt);
SYMPOL_PrintMatrix(FileOutExt, eRecCUT.EXT);


FileOutGrp:=Concatenation("CUT_", String(n), ".grp");
RemoveFileIfExist(FileOutGrp);
SYMPOL_PrintGroup(FileOutGrp, Length(eRecCUT.EXT), eRecCUT.GRP);




