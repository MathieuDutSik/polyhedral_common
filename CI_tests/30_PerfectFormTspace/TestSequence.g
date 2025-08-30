Read("../common.g");
Print("Beginning Test enumeration of perfect forms in T-spaces\n");

GetFundamentalInfo:=function(d)
  local res, IsCorrect, eSum, eProd, Dval, eQuot;
  res:=d mod 4;
  IsCorrect:=false;
  eSum:=0;
  eProd:=0;
  if res=0 then
    eSum:=0;
    eProd:=-d/4;
    Dval:=-eProd;
    IsCorrect:=true;
    if IsInt((Dval-1)/4)=true or IsInt(Dval/4)=true then
      IsCorrect:=false;
    fi;
  fi;
  if res=1 then
    eQuot:=(1-d)/4;
    eSum:=1;
    eProd:=eQuot;
    IsCorrect:=true;
  fi;
  return rec(eSum:=eSum, eProd:=eProd, IsCorrect:=IsCorrect);
end;

GetRecInfo:=function(d, n)
    local FileNml, FileResult, output, eProg, TheCommand, U, is_correct, info;
    #
    FileNml:="PerfectForm.nml";
    FileResult:="Result";
    RemoveFileIfExist(FileResult);
    #
    info:=GetFundamentalInfo(d);
    if info.IsCorrect=false then
        Print("Discriminant d=", d, " is not valid, skipping\n");
        return fail;
    fi;
    #
    # Create the namelist file
    output:=OutputTextFile(FileNml, true);
    AppendTo(output, "&DATA\n");
    AppendTo(output, " arithmetic_T=\"gmp_rational\"\n");
    AppendTo(output, " arithmetic_Tint=\"gmp_integer\"\n");
    AppendTo(output, " OutFormat=\"ObjectGAP\"\n");
    AppendTo(output, " OutFile=\"", FileResult, "\"\n");
    AppendTo(output, " FileDualDescription=\"../../ExternalLib/DualDescriptionStandard.txt\"\n");
    AppendTo(output, " max_runtime_second=0\n");
    AppendTo(output, " ApplyStdUnitbuf=F\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&STORAGE\n");
    AppendTo(output, " Saving=F\n");
    AppendTo(output, " Prefix=\"/tmp/PerfectForm/\"\n");
    AppendTo(output, "/\n");
    AppendTo(output, "\n");
    AppendTo(output, "&TSPACE\n");
    AppendTo(output, " TypeTspace=\"RealImaginary\"\n");
    AppendTo(output, " FileLinSpa=\"unset.linspa\"\n");
    AppendTo(output, " SuperMatMethod=\"NotNeeded\"\n");
    AppendTo(output, " ListComm=\"Trivial\"\n");
    AppendTo(output, " PtGroupMethod=\"Trivial\"\n");
    AppendTo(output, " FileListSubspaces=\"unset\"\n");
    AppendTo(output, " RealImagDim=", n, "\n");
    AppendTo(output, " RealImagSum=", info.eSum, "\n");
    AppendTo(output, " RealImagProd=", info.eProd, "\n");
    AppendTo(output, "/\n");
    CloseStream(output);
    #
    eProg:="../../src_perfect/PERF_SerialEnumeratePerfectCones";
    TheCommand:=Concatenation(eProg, " ", FileNml);
    Exec(TheCommand);
    #
    if IsExistingFile(FileResult)=false then
        Error("The output file is not existing. That qualifies as a fail");
    fi;
    U:=ReadAsFunction(FileResult)();
    RemoveFile(FileNml);
    RemoveFile(FileResult);
    return U;
end;

# Test cases for n=3 with various discriminants
ListDiscriminants_n3:=[-3, -4, -7, -8, -11, -15, -19, -20, -23, -24];
# Test cases for n=4 with discriminants -3, -4
ListDiscriminants_n4:=[-3, -4];

Print("Testing perfect forms for n=3\n");
ListRec_n3:=[];
for d in ListDiscriminants_n3
do
    Print("  ----------\n");
    Print("Testing d=", d, " n=3\n");
    U:=GetRecInfo(d, 3);
    if U<>fail then
        nb:=Length(U);
        Print("d=", d, " n=3: Found ", nb, " perfect forms\n");
        eRec:=rec(d:=d, n:=3, nb:=nb);
        Add(ListRec_n3, eRec);
    fi;
od;

Print("Testing perfect forms for n=4\n");
ListRec_n4:=[];
for d in ListDiscriminants_n4
do
    Print("  ----------\n");
    Print("Testing d=", d, " n=4\n");
    U:=GetRecInfo(d, 4);
    if U<>fail then
        nb:=Length(U);
        Print("d=", d, " n=4: Found ", nb, " perfect forms\n");
        eRec:=rec(d:=d, n:=4, nb:=nb);
        Add(ListRec_n4, eRec);
    fi;
od;

Print("Summary for n=3:\n");
for eRec in ListRec_n3
do
    Print("d=", eRec.d, " nb=", eRec.nb, "\n");
od;

Print("Summary for n=4:\n");
for eRec in ListRec_n4
do
    Print("d=", eRec.d, " nb=", eRec.nb, "\n");
od;

Print("Test enumeration of perfect forms in T-spaces completed\n");

# Create the CI_CONCLUSION file to indicate successful completion
output:=OutputTextFile("CI_CONCLUSION", true);
AppendTo(output, "CI test 30_PerfectFormTspace completed successfully\n");
AppendTo(output, "Results summary:\n");
AppendTo(output, "n=3 cases tested: ", Length(ListRec_n3), "\n");
AppendTo(output, "n=4 cases tested: ", Length(ListRec_n4), "\n");
AppendTo(output, "Total perfect forms enumerated:\n");
for eRec in ListRec_n3
do
    AppendTo(output, "d=", eRec.d, " n=", eRec.n, " nb=", eRec.nb, "\n");
od;
for eRec in ListRec_n4
do
    AppendTo(output, "d=", eRec.d, " n=", eRec.n, " nb=", eRec.nb, "\n");
od;
CloseStream(output);