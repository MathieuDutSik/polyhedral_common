# Vendored Python used only to generate reference values

`siegel_series/` is the Sage package of S. Takemori
(https://github.com/stakemori/siegel_series, GPL-3, see its LICENSE), which
implements Katsurada's recursion for Siegel series (H. Katsurada, An explicit
formula for Siegel series, Amer. J. Math. 121 (1999) 415-452). Only the
modules needed for `siegel_series_polynomial` are kept. Two Python-3 fixes,
both in `rec_coeff.py`: three divisions of integers by 2 that became float
divisions, and `(-1) ** e` with a negative integer `e`, which Python 3
evaluates to the float -1.0; with them the polynomials are exact again.

`degree2/utils.py` replaces the one function (`list_group_by`) the package
imports from Takemori's `degree2` package.

The C++ port is `src_genera/SiegelSeries.h`; these files are used by
`../make_density_references.py` to produce the reference values the GAP test
`TestLocalDensity.g` compares it against.
