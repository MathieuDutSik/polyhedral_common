The edgewalk of Allcock, "An alternative to Vinberg's algorithm", the
programs being built with flint and run on it.

`TestReflectivity.g`: LORENTZ_FundDomain_AllcockEdgewalk on the reflective
Lorentzian lattices of ListReflect, the number of simple roots being
compared with the known one.

`TestIsomorphism.g`: LORENTZ_FundDomain_AllcockEdgewalk_Isomorphism on
lattices of ListReflect in dimension 3 to 8. A lattice and a scrambled copy
must be found isometric, the isometry being checked; two lattices of
different determinants, or of equal determinant and different parity, must
not.

`TestTwoDim.g`: the algorithms of Section 8 of the paper for the
2-dimensional Lorentzian forms a x^2 + 2 b x y + c y^2, by
LORENTZ_TwoDimAnisotropic_Allcock, each result being checked by enumerating
the lattice vectors of a bounded region.
- TwoDimAnisotropicCases (a b c M): the generator g of SO^+(L) and the list
  R of Algorithm 8.7 from the start vector r0. g preserves the form and is
  not +-1, and R is exactly the set of the primitive vectors of norm <= M
  after r0 and before g(R[1]), in ascending order. The forms with a <= 0 and
  c <= 0, where the start vector comes from the positive eigenvector, are
  about half of the cases.
- TwoDimFirstNextCases (a b c x y N): the first vector of norm N after
  r0 = (x,y) in its sector, for isotropic and anisotropic forms, multiples
  of primitive vectors included.
The cases were drawn at random among the forms whose region to enumerate
fits in a box of side 60.

Note on the paper: step 2 of Algorithm 8.2 (Promised) reads "m^2 <= M or
m.l < 0" and calls Isotropic; the code tests m.r < 0 and calls Promised.
Compared with an enumeration on random anisotropic forms, the test of the
paper gives wrong vectors and the one of the code does not.
