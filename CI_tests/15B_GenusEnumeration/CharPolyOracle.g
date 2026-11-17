# The oracle of the tests of GENUS_GlueIsometry factors: for a genus given
# by its classes Lambda and generators of the O(Lambda), the sum over the
# classes and over the sigma in O(Lambda) of characteristic polynomial
# prod_d Phi_d^{b_d} of 1 / |O(Lambda)|, that is the total mass of the
# genera of pairs (L, sigma) of that characteristic polynomial.
#
# Input (InFile): return rec(ListGram := [...], ListListGen := [...]), the
# generators acting on the rows, g A g^T = A, or on the columns.
# Output (OutFile): one line "sig mass" per characteristic polynomial,
# sig = d1^b1,d2^b2,... by increasing d.
#
# Usage: gap -q CharPolyOracle.g with InFile and OutFile bound first, e.g.
#   echo 'InFile:="in.g";; OutFile:="out.txt";; Read("CharPolyOracle.g");' | gap -q
x := Indeterminate(Rationals, "x");;
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
  M := Group(gens);
  # a faithful permutation action: the vectors of norm at most the largest
  # diagonal entry, which span
  sv := ShortestVectors(A, Maximum(DiagonalOfMat(A))).vectors;
  sv := Concatenation(sv, -sv);
  Sort(sv);
  hom := ActionHomomorphism(M, sv, OnRight);
  P := Image(hom);
  for c in ConjugacyClasses(P) do
    g := PreImagesRepresentative(hom, Representative(c));
    s := CharPolySignature(g);
    if not IsBound(masses.(s)) then
      masses.(s) := 0;
    fi;
    masses.(s) := masses.(s) + Size(c) / Size(P);
  od;
od;
out := OutputTextFile(OutFile, false);;
SetPrintFormattingStatus(out, false);;
for s in Set(RecNames(masses)) do
  PrintTo(out, s, " ", masses.(s), "\n");
od;
CloseStream(out);
QUIT;
