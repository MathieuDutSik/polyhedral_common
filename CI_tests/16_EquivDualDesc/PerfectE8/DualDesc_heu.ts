&HEURISTIC_PRIOR
  DefaultPrior = "state_normaliz"
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
  ListDescription = "superfine", "0-1,2-infinity"
  ListKey = "incidence", "delta"
/


&PROBABILITY_DISTRIBUTIONS
  ListNmax = 25, 1
  ListNstart = 2, 1
  ListDescription = "0.0", "0.0"
  ListName = "distri1", "distriTriv"
  ListNature = "dirac", "dirac"
/


&THOMPSON_PRIOR
  ListAnswer = "small_polytopes", "rs_cdd", "normaliz"
  ListDescription = "small_polytopes:distriTriv", "rs:distri1 cdd:distri1", "normaliz:distriTriv"
  ListName = "state_small", "state_opts", "state_normaliz"
/
