&HEURISTIC_PRIOR
  DefaultPrior = "state_lrs"
  ListConclusion = "state_small", "state_opts"
  ListFullCond = "delta <= 1", "delta < 16"
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
  ListNmax = 25, 1
  ListNstart = 2, 1
  ListDescription = "0.0", "0.0"
  ListName = "distri1", "distriTriv"
  ListNature = "dirac", "dirac"
/


&THOMPSON_PRIOR
  ListAnswer = "small_polytopes", "lrs_cdd", "lrs"
  ListDescription = "small_polytopes:distriTriv", "lrs:distri1 cdd:distri1", "lrs:distriTriv"
  ListName = "state_small", "state_opts", "state_lrs"
/
