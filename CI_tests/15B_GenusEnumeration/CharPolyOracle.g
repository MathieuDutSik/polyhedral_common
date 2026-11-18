# The oracle of the tests of GENUS_GlueIsometry factors: for a genus given
# by its classes Lambda and generators of the O(Lambda),
# - by default, the sum over the classes and over the sigma in O(Lambda) of
#   characteristic polynomial prod_d Phi_d^{b_d} of 1 / |O(Lambda)|, the
#   total mass of the genera of pairs (L, sigma) of that characteristic
#   polynomial; output lines "sig mass", sig = d1^b1,d2^b2,... by
#   increasing d;
# - with Klein := true, for G = Z/2 x Z/2, the sum over the classes and
#   over the ordered pairs (g1, g2) of commuting elements of O(Lambda) with
#   g1^2 = g2^2 = 1 of 1 / |O(Lambda)|, that is the sum over the classes
#   [g1] and the classes [g2] of C(g1) of 1 / |C(g1, g2)|, by the dimensions
#   of the joint eigenspaces; output lines "m++,m+-,m-+,m-- mass".
#
# Input (InFile): return rec(ListGram := [...], ListListGen := [...]), the
# generators acting on the rows, g A g^T = A, or on the columns.
#
# Usage: gap -q with InFile, OutFile (and Klein) bound first, e.g.
#   echo 'InFile:="in.g";; OutFile:="out.txt";; Read("CharPolyOracle.g");' | gap -q
x := Indeterminate(Rationals, "x");;
if not IsBound(Klein) then
  Klein := false;
fi;
CharPolySignature := function(g)
  local P, n, sig, d, Q, b;
  P := CharacteristicPolynomial(g);
  n := Length(g);
  sig := [];
  for d in [1 .. 4 * n * n] do
    if Phi(d) <= n then
      Q := CyclotomicPolynomial(Rationals, d);
      b := 0;
      while IsZero(P mod Q) do
        P := P / Q;
        b := b + 1;
      od;
      if b > 0 then
        Add(sig, Concatenation(String(d), "^", String(b)));
      fi;
    fi;
  od;
  if not IsOne(P) then
    Error("not a product of cyclotomic polynomials");
  fi;
  return JoinStringsWithSeparator(sig, ",");
end;;
# the dimensions of the joint eigenspaces of g1 and g2 (involutions)
JointSignature := function(g1, g2)
  local n, Id, dims, s1, s2, K;
  n := Length(g1);
  Id := IdentityMat(n);
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
  return JoinStringsWithSeparator(List(dims, String), ",");
end;;
R := ReadAsFunction(InFile)();;
masses := rec();;
AddMass := function(s, m)
  if not IsBound(masses.(s)) then
    masses.(s) := 0;
  fi;
  masses.(s) := masses.(s) + m;
end;;
for i in [1 .. Length(R.ListGram)] do
  A := R.ListGram[i];
  gens := R.ListListGen[i];
  if ForAll(gens, g -> TransposedMat(g) * A * g = A) then
    gens := List(gens, TransposedMat);
  fi;
  if not ForAll(gens, g -> g * A * TransposedMat(g) = A) then
    Error("the generators are not isometries");
  fi;
  M := Group(gens);
  # a faithful permutation action: the vectors of norm at most the largest
  # diagonal entry, which span
  sv := ShortestVectors(A, Maximum(DiagonalOfMat(A))).vectors;
  sv := Concatenation(sv, -sv);
  Sort(sv);
  hom := ActionHomomorphism(M, sv, OnRight);
  P := Image(hom);
  for c1 in ConjugacyClasses(P) do
    p1 := Representative(c1);
    g1 := PreImagesRepresentative(hom, p1);
    if not Klein then
      AddMass(CharPolySignature(g1), Size(c1) / Size(P));
    elif Order(p1) <= 2 then
      C1 := Centralizer(P, p1);
      for c2 in ConjugacyClasses(C1) do
        p2 := Representative(c2);
        if Order(p2) <= 2 then
          g2 := PreImagesRepresentative(hom, p2);
          AddMass(JointSignature(g1, g2), 1 / Size(Centralizer(C1, p2)));
        fi;
      od;
    fi;
  od;
od;
out := OutputTextFile(OutFile, false);;
SetPrintFormattingStatus(out, false);;
for s in Set(RecNames(masses)) do
  PrintTo(out, s, " ", masses.(s), "\n");
od;
CloseStream(out);
QUIT;
