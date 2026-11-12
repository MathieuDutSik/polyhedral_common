#!/bin/bash
# The programs of the genera of pairs with every arithmetic compiled in
# (arithmetic=gmp, flint, gmp_boost, multi_boost): the masses must be those
# of gmp, and the enumerations must close with the same number of classes.
# An arithmetic that was not compiled is skipped.
set -e
MASS=../../src_genera/GENUS_IsometryMass
ENUM=../../src_genera/GENUS_EnumerateIsometry

# lattice, isometry or group file, mode (mass or group)
MASS_CASES="
E8 refl mass
E8 E8_inv3 mass
E6 E6_g24 mass
det30_L557 det30_g19 mass
A6A6 A6A6_F21 group
A4A4A4 A4A4A4_D5 group
qD4x2 qD4x2_Q8 group
"
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
    $MASS hermitian 3 1 arithmetic=$1 > /dev/null 2>&1
}

mass_of() {
    arith=$1; lat=$2; isom=$3; mode=$4
    if [ "$mode" = "mass" ]; then
        $MASS mass isom_${lat}.txt isom_${isom}.txt none stdout \
              arithmetic=$arith 2>/dev/null
    else
        $MASS group isom_${lat}.txt isom_${isom}.txt stdout \
              arithmetic=$arith 2>/dev/null
    fi | sed -n 's/^  mass:=\(.*\));$/\1/p'
}

n_arith=0
for arith in gmp flint gmp_boost multi_boost; do
    if ! compiled $arith; then
        echo "--- arithmetic=$arith: not compiled, skipped ---"
        continue
    fi
    n_arith=$((n_arith + 1))
    echo "--- arithmetic=$arith ---"
    echo "$MASS_CASES" | while read lat isom mode; do
        [ -z "$lat" ] && continue
        ref=$(mass_of gmp $lat $isom $mode)
        val=$(mass_of $arith $lat $isom $mode)
        echo "mass ($lat, $isom): $val"
        if [ -z "$ref" ] || [ "$val" != "$ref" ]; then
            echo "ERROR: ($lat, $isom) with $arith gave '$val', gmp '$ref'" >&2
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
