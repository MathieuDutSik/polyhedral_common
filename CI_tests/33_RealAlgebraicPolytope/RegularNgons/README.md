Regular N-gons over real algebraic fields
=========================================

The regular N-gon is the polytope whose vertices are `(cos(2*pi*k/N), sin(2*pi*k/N))`
for `k = 0, ..., N-1`. For most N these coordinates are real algebraic numbers, so
the polytope is given over a real algebraic field `Q(x)` with generator
`x = sin(2*pi/N)`. The vertices are stored (in homogeneous coordinates
`[1, cos, sin]`) as polynomials in `x`; each entry is a token of the form used by
`NumberTheoryRealField.h` (e.g. `-3/2+2*x^2`).

Files
-----

For each N there are three files:
  * `generate_file_desc_<N>.sage` : Sage script producing the field description.
  * `FileDesc<N>`                 : the real algebraic field description (output of
                                    the script), used in the type input
                                    `RealAlgebraic=FileDesc<N>`.
  * `Regular<N>gon`               : the N-gon vertex matrix over that field.

The `FileDesc<N>` file contains:
  * the degree of the minimal polynomial of `x`,
  * the minimal polynomial (coefficients in ascending degree order),
  * a `double` approximation of `x`,
  * a sequence of lower/upper continued-fraction approximants bracketing `x`.

To regenerate a description file: `sage generate_file_desc_<N>.sage`.


N = 5  (x = sin(2*pi/5))
------------------------

Minimal polynomial: `64 x^4 - 80 x^2 + 20 = 0`.

  cos(0)      = 1
  cos(2*pi/5) = 2*x^2 - 3/2
  cos(4*pi/5) = 1 - 2*x^2
  cos(6*pi/5) = cos(4*pi/5) = 1 - 2*x^2
  cos(8*pi/5) = cos(2*pi/5) = 2*x^2 - 3/2

  sin(2*pi/5) = x
  sin(4*pi/5) = 4*x^3 - 3*x
  (sin(6*pi/5) = -sin(4*pi/5), sin(8*pi/5) = -sin(2*pi/5))


N = 7  (x = sin(2*pi/7))
------------------------

Minimal polynomial: `64 x^6 - 112 x^4 + 56 x^2 - 7 = 0`.

  cos(0)       = 1
  cos(2*pi/7)  = -8*x^4 + 10*x^2 - 5/2
  cos(4*pi/7)  = 1 - 2*x^2
  cos(6*pi/7)  = 8*x^4 - 8*x^2 + 1
  cos(8*pi/7)  = cos(6*pi/7)  = 8*x^4 - 8*x^2 + 1
  cos(10*pi/7) = cos(4*pi/7)  = 1 - 2*x^2
  cos(12*pi/7) = cos(2*pi/7)  = -8*x^4 + 10*x^2 - 5/2

  sin(2*pi/7)  = x
  sin(4*pi/7)  = -16*x^5 + 20*x^3 - 5*x
  sin(6*pi/7)  = 3*x - 4*x^3
  (sin(8*pi/7) = -sin(6*pi/7), sin(10*pi/7) = -sin(4*pi/7), sin(12*pi/7) = -sin(2*pi/7))


N = 9  (x = sin(2*pi/9))
------------------------

Minimal polynomial: `64 x^6 - 96 x^4 + 36 x^2 - 3 = 0`.

  cos(0)       = 1
  cos(2*pi/9)  = -8*x^4 + 10*x^2 - 2
  cos(4*pi/9)  = 1 - 2*x^2
  cos(6*pi/9)  = -1/2                       (= cos(2*pi/3))
  cos(8*pi/9)  = 8*x^4 - 8*x^2 + 1
  cos(10*pi/9) = cos(8*pi/9) = 8*x^4 - 8*x^2 + 1
  cos(12*pi/9) = cos(6*pi/9) = -1/2
  cos(14*pi/9) = cos(4*pi/9) = 1 - 2*x^2
  cos(16*pi/9) = cos(2*pi/9) = -8*x^4 + 10*x^2 - 2

  sin(2*pi/9)  = x
  sin(4*pi/9)  = -16*x^5 + 20*x^3 - 4*x
  sin(6*pi/9)  = 3*x - 4*x^3                (= sqrt(3)/2)
  sin(8*pi/9)  = -16*x^5 + 20*x^3 - 5*x
  (sin(2*pi*(9-k)/9) = -sin(2*pi*k/9))


How to use it
-------------

The field is selected by the arithmetic type `RealAlgebraic=FileDesc<N>`, where
`FileDesc<N>` is read from the current directory. Examples:

Automorphism group of the polytope (the dihedral group of order 2N):

  ./GRP_LinPolytope_Automorphism RealAlgebraic=FileDesc5 Regular5gon
  # -> |GRP|=10   (dihedral D_5)

  ./GRP_LinPolytope_Automorphism RealAlgebraic=FileDesc7 Regular7gon
  # -> |GRP|=14   (dihedral D_7)

  ./GRP_LinPolytope_Automorphism RealAlgebraic=FileDesc9 Regular9gon
  # -> |GRP|=18   (dihedral D_9)

Dual description (facets); a regular N-gon has N facets:

  ./POLY_dual_description RealAlgebraic=FileDesc5 cdd CPP Regular5gon out5
  ./POLY_dual_description RealAlgebraic=FileDesc7 cdd CPP Regular7gon out7
  ./POLY_dual_description RealAlgebraic=FileDesc9 cdd CPP Regular9gon out9

Note: for a real algebraic field use the output format `CPP` (or `control`); the
`GAP` output format does not support algebraic field entries.
