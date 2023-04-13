&HEURISTIC_PRIOR
  DefaultPrior = "state_ppl_ext"
  ListConclusion = "state_lrs_ring", "state_lrs_ring", "state_ppl_ext"
  ListFullCond = "incidence < 30", "incidence < 35", "incidence < 40"
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
  ListAnswer = "lrs_ring", "ppl_ext"
  ListDescription = "lrs_ring:distri1", "ppl_ext:distri1"
  ListName = "state_lrs_ring", "state_ppl_ext"
/
