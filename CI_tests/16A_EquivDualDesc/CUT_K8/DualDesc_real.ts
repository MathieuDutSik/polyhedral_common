&HEURISTIC_PRIOR
  DefaultPrior = "state_cdd"
  ListConclusion = "state_opts"
  ListFullCond = "delta < 17"
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
  ListAnswer = "lrsring_cdd", "cdd"
  ListDescription = "lrs:distri1 cdd:distri1", "cdd:distriTriv"
  ListName = "state_opts", "state_cdd"
/
