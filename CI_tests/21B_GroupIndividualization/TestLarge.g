# The cases above the thresholds, with the binary built as in the other
# tests.
Read("IndividualizationFunctions.g");
Print("Beginning TestLarge\n");
ListRec:=[CaseWythoffH4(),
          CaseOrbit("B6 regular", 6, [1, 2, 3, 5, 8, 13] / 7),
          CaseOrbit("B6 stab2", 6, [1, 1, 2, 3, 4, 5]),
          CaseSmallSubset(6, [1, 2, 3, 5, 8, 13] / 7, 7)];
RunIndividualizationTests(ListRec);
