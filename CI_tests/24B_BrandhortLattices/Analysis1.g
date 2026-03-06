Read("/Users/mathieudutoursikiric/GITall/GITmathieu/polyhedral_common/CI_tests/common.g");
Read("/Users/mathieudutoursikiric/GITall/GITmathieu/polyhedral_common/CI_tests/access_points.g");



LMat:=ReadAsFunction("finite_symmetry.gap")();


for iMat in [1..Length(LMat)]
do
    eMat:=LMat[iMat];
    Print("iMat=", iMat, " det(eMat)=", DeterminantMat(eMat), "\n");
    LorMat:=-eMat;
    RecGRP:=GetReflectivityInformation(LorMat);
    GRPmatr:=RecGRP.GrpIsomCoxMatr;
    Print("  |GRPmatr|=", Order(GRPmatr), "\n");
od;
