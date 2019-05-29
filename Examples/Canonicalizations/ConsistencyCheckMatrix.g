


RandomCheckProducts:=function(eG, nbIter)
  local eGcan, n, iIter, Umat, eProd, eProdCan;
  eGcan:=PositiveDefiniteCanonicalization(eG);
  n:=Length(eG);
  for iIter in [1..nbIter]
  do
    Print("iIter=", iIter, "\n");
    Umat:=RandomIntegralGLnZmatrix(n);
    eProd:=Umat * eG * TransposedMat(Umat);
    eProdCan:=PositiveDefiniteCanonicalization(eProd);
    if eProdCan <> eGcan then
      Print("The canonicalization procedure gives incoherent results\n");
      Error("Please debug");
    fi;
  od;
end;




GramMat:=ClassicalSporadicLattices("E8");
nbIter:=1000;
RandomCheckProducts(GramMat, nbIter);
