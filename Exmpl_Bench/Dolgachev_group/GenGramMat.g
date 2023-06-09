LoadPackage("grape");


ListEdges:=[];
Nvert:=1;
AppendLine:=function(k)
    local i;
    Add(ListEdges, [1, Nvert + 1]);
    for i in [1..k-1]
    do
        Add(ListEdges, [Nvert + i, Nvert + i + 1]);
    od;
    Nvert:=Nvert+k;
end;

AppendLine(2);
AppendLine(4);
AppendLine(4);


TheMat:=NullMat(Nvert, Nvert);
for i in [1..Nvert]
do
    TheMat[i][i]:=2;
od;
for eEdge in ListEdges
do
    i:=eEdge[1];
    j:=eEdge[2];
    TheMat[i][j]:=-1;
    TheMat[j][i]:=-1;
od;


TheFile:="GramMat";
output:=OutputTextFile(TheFile, true);
AppendTo(output, Nvert, " ", Nvert, "\n");
for i in [1..Nvert]
do
    for j in [1..Nvert]
    do
        AppendTo(output, " ", TheMat[i][j]);
    od;
    AppendTo(output, "\n");
od;
CloseStream(output);
