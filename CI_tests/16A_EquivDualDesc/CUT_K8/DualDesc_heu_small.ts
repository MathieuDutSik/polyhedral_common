&HEURISTIC_PRIOR
  DefaultPrior = "state_cdd"
  ListConclusion = "state_small_polytopes", "state_lrs"
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
  ListAnswer = "small_polytopes", "lrs", "cdd"
  ListDescription = "small_polytopes:distri1", "lrs:distri1", "cdd:distri1"
  ListName = "state_small_polytopes", "state_lrs", "state_cdd"
/
