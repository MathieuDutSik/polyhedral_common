#!/bin/bash
# The programs of the genera of pairs with every arithmetic compiled in
# (arithmetic=gmp, flint, gmp_boost, multi_boost): the masses must be those
# of gmp, and the enumerations must close with the same number of classes.
# An arithmetic that was not compiled is skipped.
set -e
. ./MassPrograms.sh
ENUM=../../src_genera/GENUS_EnumerateIsometry

# lattice, isometry or group file
MASS_CASES="
E8 refl
E8 E8_inv3
E6 E6_g24
det30_L557 det30_g19
A6A6 A6A6_F21
A4A4A4 A4A4A4_D5
qD4x2 qD4x2_Q8
"

# the masses of MassAbelianHermitian.h and MassAbelianQuadratic.h, generic
# in the rational type (GENUS_ComponentMass): the field (abelian_cm or
# abelian_real), the files tag_gram.txt and tag_omega.txt, n and the
# generators of H
ABMASS_CASES="
abelian_cm hab_q7r3 7 1,2,4
abelian_cm hab_s17m2 136 9,19
abelian_real quadm_z7pm2 7 1
abelian_real quad_s2m4 8 7
abelian_real quad_s17s41m3 697 2,21
"

abmass_of() {
    arith=$1; field=$2; tag=$3; n=$4; H=$5
    output=mass
    [ "$field" = abelian_real ] && output=local_factors
    component_mass field=$field output=$output GramFile=${tag}_gram.txt \
          OmegaFile=${tag}_omega.txt n=$n H=$H arithmetic=$arith \
          2>/dev/null | sed -n 's/^return rec(mass:=\([^,)]*\).*/\1/p'
}
# lattice, isometries, expected classes
ENUM_CASES="
d25 refl8 4
E8E8 c6 3
E6x2 A2diagx2 2
A8 F21A8 2
qD4W2glue qD4W2glueQS 3
"

compiled() {
    # the dispatch refuses a name that was not compiled
    component_mass field=cyclotomic n=3 rank=1 arithmetic=$1 > /dev/null 2>&1
}

mass_of() {
    arith=$1; lat=$2; isom=$3
    isometry_mass GramFile=isom_${lat}.txt GeneratorsFile=isom_${isom}.txt \
          arithmetic=$arith 2>/dev/null | sed -n 's/^  mass:=\(.*\));$/\1/p'
}

n_arith=0
for arith in gmp flint gmp_boost multi_boost; do
    if ! compiled $arith; then
        echo "--- arithmetic=$arith: not compiled, skipped ---"
        continue
    fi
    n_arith=$((n_arith + 1))
    echo "--- arithmetic=$arith ---"
    echo "$MASS_CASES" | while read lat isom; do
        [ -z "$lat" ] && continue
        ref=$(mass_of gmp $lat $isom)
        val=$(mass_of $arith $lat $isom)
        echo "mass ($lat, $isom): $val"
        if [ -z "$ref" ] || [ "$val" != "$ref" ]; then
            echo "ERROR: ($lat, $isom) with $arith gave '$val', gmp '$ref'" >&2
            exit 1
        fi
    done
    echo "$ABMASS_CASES" | while read field tag n H; do
        [ -z "$field" ] && continue
        ref=$(abmass_of gmp $field $tag $n $H)
        val=$(abmass_of $arith $field $tag $n $H)
        echo "mass ($tag): $val"
        if [ -z "$ref" ] || [ "$val" != "$ref" ]; then
            echo "ERROR: $tag with $arith gave '$val', gmp '$ref'" >&2
            exit 1
        fi
    done
    echo "$ENUM_CASES" | while read lat isom expected; do
        [ -z "$lat" ] && continue
        out=$($ENUM eiso_${lat}.txt eiso_${isom}.txt eiso_mass_${isom}.txt \
                  Summary stdout arithmetic=$arith 2>/dev/null)
        nb=$(echo "$out" | sed -n 's/^classes of pairs = //p')
        complete=$(echo "$out" | sed -n 's/^complete = //p')
        echo "enumeration ($lat, $isom): $nb classes, complete = $complete"
        if [ "$nb" != "$expected" ] || [ "$complete" != "true" ]; then
            echo "ERROR: ($lat, $isom) with $arith gave $nb classes" >&2
            exit 1
        fi
    done
done
if [ "$n_arith" = "0" ]; then
    echo "ERROR: no arithmetic available" >&2
    exit 1
fi

echo "All arithmetic tests passed"
touch CI_CONCLUSION
