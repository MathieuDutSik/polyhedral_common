import os

# Generate the real algebraic field description file used for the type input
# "RealAlgebraic=FileDesc5" in polyhedral_common. The generator field element is
#   x = sin(2*pi/5)
# whose minimal polynomial is
#   64 x^4 - 80 x^2 + 20 = 0
# (l_coeff below lists the coefficients in ascending degree order).
#
# The regular 5-gon (Regular5gon) is written in terms of this x using:
#   sin(2pi/5) = x
#   cos(2pi/5) = 2x^2 - 3/2
#   sin(4pi/5) = 4x^3 - 3x
#   cos(4pi/5) = 1 - 2x^2
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


# x = sin(2*pi/5), minimal polynomial 64 x^4 - 80 x^2 + 20 = 0.
val = sin(2*pi/5)
l_coeff = [20, 0, -80, 0, 64]
FileName = "FileDesc5"
create_real_algebraic_input(val, l_coeff, FileName)
