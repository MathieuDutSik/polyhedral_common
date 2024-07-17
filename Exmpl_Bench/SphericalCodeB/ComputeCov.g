
for i in [1..37]
do
    FileIn :=Concatenation("codeword.", String(i), ".txt");
    FileOut:=Concatenation("codeword.", String(i), ".out");
    Print("FileIn=", FileIn, "\n");
    eProg:="../../src_spherical_code/CODE_Analysis";
    TheCommand:=Concatenation(eProg, " Qsqrt2 ", FileIn, " 2> ", FileOut);
    Exec(TheCommand);
    #
od;
