import os

# Generate the real algebraic field description file used for the type input
# "RealAlgebraic=FileDesc9" in polyhedral_common. The generator field element is
#   x = sin(2*pi/9)
# whose minimal polynomial is
#   64 x^6 - 96 x^4 + 36 x^2 - 3 = 0
# (l_coeff below lists the coefficients in ascending degree order).
#
# The regular 9-gon (Regular9gon) is written in terms of this x using:
#   sin(2pi/9) = x
#   cos(2pi/9) = -8x^4 + 10x^2 - 2
#   sin(4pi/9) = -16x^5 + 20x^3 - 4x
#   cos(4pi/9) = 1 - 2x^2
#   sin(6pi/9) = 3x - 4x^3            (= sqrt(3)/2)
#   cos(6pi/9) = -1/2
#   sin(8pi/9) = -16x^5 + 20x^3 - 5x
#   cos(8pi/9) = 8x^4 - 8x^2 + 1
#
# The description file contains:
#   * the degree of the minimal polynomial
#   * the minimal polynomial (ascending coefficients)
#   * a double approximation of the value
#   * a sequence of lower/upper continued-fraction approximants bracketing x


def create_real_algebraic_input(val, l_coeff, FileName):

    def get_error(the_rat):
        thediff = abs(val - the_rat)
        thepow = 1
        while(True):
            thepow_new = thepow * 2
            if thediff > 1/thepow_new:
                return 1/thepow
            thepow = thepow_new


    def get_estimate_error(the_val):
        thepow = 1
        expo = 0
        while(True):
            thepow_new = thepow * 2
            if the_val > 1/thepow_new:
                return expo
            thepow = thepow_new
            expo += 1


    def get_approx_order(n):
        elist = continued_fraction_list(val, nterms=n)
        cf = continued_fraction(elist)
        the_rat = cf.value()
        return the_rat

    def get_approximant(n):
        approx1 = get_approx_order(n)
        approx2 = get_approx_order(n+1)
        if approx1 < approx2:
            the_diff = approx2 - approx1
            expo = get_estimate_error(the_diff)
            return [approx1, approx2, expo]
        else:
            the_diff = approx1 - approx2
            expo = get_estimate_error(the_diff)
            return [approx2, approx1, expo]

    f = open(FileName, "w")
    # the degree
    the_deg = len(l_coeff) - 1
    f.write(str(the_deg) + "\n")
    # the minimal polynomial
    for i in range(len(l_coeff)):
        if i>0:
            f.write(" ")
        f.write(str(l_coeff[i]))
    f.write("\n");
    # The double approximation of the value
    f.write(str(float(val)) + "\n");
    # the rational approximations
    n_expo = 100
    f.write(str(n_expo) + "\n")
    for i in range(1,n_expo+1):
        print("i=", i, " / ", n_expo)
        [approx1, approx2, expo] = get_approximant(5*i)
        if approx1 > val or approx2 < val:
            print("BIG ERROR")
            os.sys.exit(1)
        print("approx1=", approx1, " approx2=", approx2, " expo=", expo)
        f.write(str(approx1) + " " + str(approx2) + "\n")

    f.close()


# x = sin(2*pi/9), minimal polynomial 64 x^6 - 96 x^4 + 36 x^2 - 3 = 0.
val = sin(2*pi/9)
l_coeff = [-3, 0, 36, 0, -96, 0, 64]
FileName = "FileDesc9"
create_real_algebraic_input(val, l_coeff, FileName)
