#!/bin/bash
# The seeds of the enumeration of pairs by gluing (GENUS_GlueIsometry): the
# sigma-stable primitive extensions of (L_1, g_1) (+) (L_2, g_2), each
# genus of pairs met enumerated by its mass.
#
# - A1 with -1 and E7 with 1, determinant 1: E8 with a root reflection,
#   one class, centraliser 2 |W(E7)|;
# - D4 with 1 and D4 with -1, any determinant: the glue orders 1, 2 and 4,
#   three genera of one class each, the last one E8 with an involution of
#   type D4 (centraliser 6 |W(D4)|^2);
# - A1 with -1 and E7+E8 with 1, determinant 1: the genus of II_16 with a
#   root reflection, whose second class (on D16+) the enumeration finds;
# - A2 with an element of order 3 and E6+E8 with 1, determinant 1: the
#   genus of II_16 with an element of order 3 of an A2, two classes;
# - D4+D4 with 1 and with -1, determinant 1: glue (Z/2)^4, one class.
# In every case all the genera must close their mass: complete = true.
set -e
BIN=../../src_genera/GENUS_GlueIsometry

glue_case() {
    l1=$1
    i1=$2
    l2=$3
    i2=$4
    det=$5
    expected_genera=$6
    expected_classes=$7
    echo "--- ($l1, $i1) + ($l2, $i2), det $det, expecting $expected_genera genera, $expected_classes classes ---"
    out=$($BIN enumerate glue_${l1}.txt glue_${i1}.txt glue_${l2}.txt \
              glue_${i2}.txt $det Summary stdout 2>/dev/null)
    echo "$out"
    nb_genera=$(echo "$out" | sed -n 's/^genera of pairs = //p')
    nb_classes=$(echo "$out" | sed -n 's/^classes of pairs = //p')
    complete=$(echo "$out" | sed -n 's/^complete = //p')
    if [ "$nb_genera" != "$expected_genera" ] ||
           [ "$nb_classes" != "$expected_classes" ] ||
           [ "$complete" != "true" ]; then
        echo "ERROR: ($l1, $i1) + ($l2, $i2) gave $nb_genera genera, $nb_classes classes, complete=$complete" >&2
        exit 1
    fi
}

glue_case A1 A1_mI E7 E7_id 1 1 1
glue_case D4 D4_id D4 D4_mI any 3 3
glue_case A1 A1_mI E7E8 E7E8_id 1 1 2
glue_case A2 A2_w E6E8 E6E8_id 1 1 2
glue_case D4D4 D4D4_id D4D4 D4D4_mI 1 1 1

echo "All gluing tests passed"
touch CI_CONCLUSION
