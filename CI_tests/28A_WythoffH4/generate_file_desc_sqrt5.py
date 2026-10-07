# Generates FileDescSqrt5, the real algebraic field description of Q(x) with
# x = sqrt(5), whose minimal polynomial is x^2 - 5 = 0, used by the type
# input "RealAlgebraic=FileDescSqrt5". The WythoffH4 polytope (WythoffH4.ext)
# is written with that x (entries such as 3/4+1/4*x).
#
# The format is the one of CI_tests/02B_RealAlgebraicPolytope/RegularNgons
# (generate_file_desc_5.sage):
#   * the degree of the minimal polynomial
#   * the minimal polynomial (ascending coefficients)
#   * a double approximation of the value
#   * a sequence of lower/upper continued-fraction approximants bracketing x
# The continued fraction of sqrt(5) is [2; 4, 4, 4, ...], so that the
# convergents are computed exactly, without Sage.
from fractions import Fraction
import math


def convergent(n_terms):
    # the value of [2; 4, ..., 4] with n_terms partial quotients
    quotients = [2] + [4] * (n_terms - 1)
    value = Fraction(quotients[-1])
    for q in reversed(quotients[:-1]):
        value = q + 1 / value
    return value


def is_below_sqrt5(r):
    return r * r < 5


l_coeff = [-5, 0, 1]
n_expo = 100
with open("FileDescSqrt5", "w") as f:
    f.write(str(len(l_coeff) - 1) + "\n")
    f.write(" ".join(str(c) for c in l_coeff) + "\n")
    f.write(str(math.sqrt(5)) + "\n")
    f.write(str(n_expo) + "\n")
    for i in range(1, n_expo + 1):
        a1 = convergent(5 * i)
        a2 = convergent(5 * i + 1)
        low, upp = min(a1, a2), max(a1, a2)
        if not (is_below_sqrt5(low) and not is_below_sqrt5(upp)):
            raise RuntimeError("the approximants do not bracket sqrt(5)")
        f.write(str(low) + " " + str(upp) + "\n")
