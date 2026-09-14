&HEURISTIC_PRIOR
  DefaultPrior = "state_cdd"
  ListConclusion = "state_small", "state_lrs", "state_cdd"
  ListFullCond = "delta <= 1", "incidence < 35", "incidence < 45"
/


&IO
  ProcessExistingDataIfExist = F
  WriteLog = F
  LogFileToProcess = "irrelevant"
  name = "unset"
/


&KEY_COMPRESSION
  ListDescription = "superfine", "0-1,2-infinity"
  ListKey = "incidence", "delta"
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
  ListName = "state_small", "state_lrs", "state_cdd"
/
