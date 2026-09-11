&HEURISTIC_PRIOR
  DefaultPrior = "state_normaliz"
  ListConclusion = "state_normaliz", "state_lrs"
  ListFullCond = "incidence > 44", "delta < 16"
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
  ListNmax = 1
  ListNstart = 1
  ListDescription = "0.0"
  ListName = "distriTriv"
  ListNature = "dirac"
/


&THOMPSON_PRIOR
  ListAnswer = "lrs", "normaliz"
  ListDescription = "lrs:distriTriv", "normaliz:distriTriv"
  ListName = "state_lrs", "state_normaliz"
/
