&HEURISTIC_PRIOR
  DefaultPrior = "state"
  ListConclusion = "state"
  ListFullCond = "incidence < 30"
/


&IO
  ProcessExistingDataIfExist = F
  WriteLog = F
  LogFileToProcess = "irrelevant"
  name = "unset"
/


&KEY_COMPRESSION
  ListDescription = "superfine"
  ListKey = "incidence"
/


&PROBABILITY_DISTRIBUTIONS
  ListNmax = 100
  ListNstart = 100
  ListDescription = "145.3"
  ListName = "distri1"
  ListNature = "dirac"
/


&THOMPSON_PRIOR
  ListAnswer = "lrsring_pplext_cdd"
  ListDescription = "lrs_ring:distri1 ppl_ext:distri1 cdd:distri1"
  ListName = "state"
/
