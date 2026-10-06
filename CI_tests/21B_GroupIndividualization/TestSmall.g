# The small cases, with the binary built with TEST_GROUP_THRESHOLD (the
# individualization above 100 points) and the sanity checks.
Read("IndividualizationFunctions.g");
Print("Beginning TestSmall\n");
ListRec:=[CaseOrbit("B4 regular", 4, [1, 2, 3, 5] / 7),
          CaseOrbit("B4 stab2", 4, [1, 1, 2, 3]),
          CaseOrbit("B5 regular", 5, [1, 2, 3, 5, 8] / 7),
          CaseSmallSubset(4, [1, 2, 3, 5] / 7, 7)];
RunIndividualizationTests(ListRec);
