import os

# Generate the real algebraic field description file used for the type input
# "RealAlgebraic=FileDesc7" in polyhedral_common. The generator field element is
#   x = 2*sin(2*pi/7)
# whose minimal polynomial is
#   x^6 - 7 x^4 + 14 x^2 - 7 = 0
# (l_coeff below lists the coefficients in ascending degree order). The
# generator is 2*sin(2*pi/7) and not sin(2*pi/7) so that the minimal polynomial
# is monic, that is so that x is an algebraic integer.
#
# The regular 7-gon (Regular7gon) is written in terms of this x using:
#   sin(2pi/7) = x/2
#   cos(2pi/7) = -x^4/2 + 5x^2/2 - 5/2
#   sin(4pi/7) = -x^5/2 + 5x^3/2 - 5x/2
#   cos(4pi/7) = 1 - x^2/2
#   sin(6pi/7) = 3x/2 - x^3/2
#   cos(6pi/7) = x^4/2 - 2x^2 + 1
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


# x = 2*sin(2*pi/7), minimal polynomial x^6 - 7 x^4 + 14 x^2 - 7 = 0.
val = 2*sin(2*pi/7)
l_coeff = [-7, 0, 14, 0, -7, 0, 1]
FileName = "FileDesc7"
create_real_algebraic_input(val, l_coeff, FileName)
