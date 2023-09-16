&HEURISTIC_PRIOR
  DefaultPrior = "state_ppl_ext"
  ListConclusion = "state_small_polytopes", "state_lrs_ring"
  ListFullCond = "delta < 2", "delta < 8"
/


&IO
  ProcessExistingDataIfExist = F
  WriteLog = F
  LogFileToProcess = "irrelevant"
  name = "unset"
/


&KEY_COMPRESSION
  ListDescription = "superfine"
  ListKey = "delta"
/


&PROBABILITY_DISTRIBUTIONS
  ListNmax = 100
  ListNstart = 100
  ListDescription = "145.3"
  ListName = "distri1"
  ListNature = "dirac"
/


&THOMPSON_PRIOR
  ListAnswer = "small_polytopes", "lrs_ring", "ppl_ext"
  ListDescription = "small_polytopes:distri1", "lrs_ring:distri1", "ppl_ext:distri1"
  ListName = "state_small_polytopes", "state_lrs_ring", "state_ppl_ext"
/
