


for i in [1..7]
do
    FileIn :=Concatenation("code", String(i), ".txt");
    FileOut:=Concatenation("codeword.", String(i), ".out");
    Print("FileIn=", FileIn, "\n");
    eProg:="../../src_spherical_code/CODE_Analysis";
    TheCommand:=Concatenation(eProg, " rational ", FileIn, " GramA3.txt 2> ", FileOut);
    Exec(TheCommand);
od;




for i in [8..37]
do
    FileIn :=Concatenation("codeword.", String(i), ".txt");
    FileOut:=Concatenation("codeword.", String(i), ".out");
    Print("FileIn=", FileIn, "\n");
    eProg:="../../src_spherical_code/CODE_Analysis";
    TheCommand:=Concatenation(eProg, " Qsqrt2 ", FileIn, " identity 2> ", FileOut);
    Exec(TheCommand);
    #
od;
