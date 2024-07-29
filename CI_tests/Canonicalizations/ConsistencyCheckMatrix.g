


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
            return false;
        fi;
    od;
    return true;
end;


GramMat:=ClassicalSporadicLattices("E8");
nbIter:=1000;
result:=RandomCheckProducts(GramMat, nbIter);
if result = false then
    # Error case
    Print("Error case\n")
    GAP_EXIT_CODE(1);
else
    # No error case
    Print("Normal case\n")
    GAP_EXIT_CODE(0);
fi;
