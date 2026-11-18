# The oracle of the tests of GENUS_GlueIsometry factors for G = Z/2 x Z/2:
# for a genus given by its classes Lambda and generators of the O(Lambda),
# the sum over the classes and over the ordered pairs (g1, g2) of commuting
# elements of O(Lambda) with g1^2 = g2^2 = 1 of 1 / |O(Lambda)|, by the
# dimensions of the joint eigenspaces: the total mass of the genera of pairs
# (L, <g1, g2>) of that type. It is
#   sum over the classes [g1] and the classes [g2] of C(g1) of 1 / |C(g1, g2)|.
#
# Input (InFile) as for CharPolyOracle.g. Output (OutFile): lines
# "m++,m+-,m-+,m-- mass".
R := ReadAsFunction(InFile)();;
masses := rec();;
for i in [1 .. Length(R.ListGram)] do
  A := R.ListGram[i];
  gens := R.ListListGen[i];
  if ForAll(gens, g -> TransposedMat(g) * A * g = A) then
    gens := List(gens, TransposedMat);
  fi;
  if not ForAll(gens, g -> g * A * TransposedMat(g) = A) then
    Error("the generators are not isometries");
  fi;
  n := Length(A);
  M := Group(gens);
  sv := ShortestVectors(A, Maximum(DiagonalOfMat(A))).vectors;
  sv := Concatenation(sv, -sv);
  Sort(sv);
  hom := ActionHomomorphism(M, sv, OnRight);
  P := Image(hom);
  Id := IdentityMat(n);
  for c1 in ConjugacyClasses(P) do
    p1 := Representative(c1);
    if Order(p1) > 2 then
      continue;
    fi;
    C1 := Centralizer(P, p1);
    g1 := PreImagesRepresentative(hom, p1);
    for c2 in ConjugacyClasses(C1) do
      p2 := Representative(c2);
      if Order(p2) > 2 then
        continue;
      fi;
      g2 := PreImagesRepresentative(hom, p2);
      dims := [];
      for s1 in [1, -1] do
        for s2 in [1, -1] do
          # x (g1 - s1) = 0 and x (g2 - s2) = 0
          K := NullspaceMat(TransposedMat(Concatenation(
                   TransposedMat(g1 - s1 * Id), TransposedMat(g2 - s2 * Id))));
          Add(dims, Length(K));
        od;
      od;
      if Sum(dims) <> n then
        Error("the joint eigenspaces do not add up");
      fi;
      s := JoinStringsWithSeparator(List(dims, String), ",");
      if not IsBound(masses.(s)) then
        masses.(s) := 0;
      fi;
      masses.(s) := masses.(s) + 1 / Size(Centralizer(C1, p2));
    od;
  od;
od;
out := OutputTextFile(OutFile, false);;
SetPrintFormattingStatus(out, false);;
for s in Set(RecNames(masses)) do
  PrintTo(out, s, " ", masses.(s), "\n");
od;
CloseStream(out);
QUIT;
