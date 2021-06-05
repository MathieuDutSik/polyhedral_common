n:=8;

ListMat:=GetAllPerfectForm(n);
nbMat:=Length(ListMat);

TheFile:=Concatenation("ListFileForm_", String(n));
RemoveFileIfExist(TheFile);


output:=OutputTextFile(TheFile, true);
AppendTo(output, nbMat, "\n");
#
for iMat in [1..nbMat]
do
  CPP_WriteMatrix(output, ListMat[iMat]);
od;

CloseStream(output);
