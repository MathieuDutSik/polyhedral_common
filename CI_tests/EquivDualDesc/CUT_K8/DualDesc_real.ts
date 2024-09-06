&HEURISTIC_PRIOR
  DefaultPrior = "state_ppl_ext"
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
  ListAnswer = "lrsring_pplext_cdd", "ppl_ext"
  ListDescription = "lrs_ring:distri1 ppl_ext:distri1 cdd:distri1", "ppl_ext:distriTriv"
  ListName = "state_opts", "state_ppl_ext"
/
