&HEURISTIC_PRIOR
  DefaultPrior = "state_glrs"
  ListConclusion = "state_lrs_ring"
  ListFullCond = "delta <= 6"
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
  ListAnswer = "glrs", "lrs_ring"
  ListDescription = "glrs:distri1", "lrs_ring:distri1"
  ListName = "state_glrs", "state_lrs_ring"
/
