#ifndef BALINSKI_BASIC_H
#define BALINSKI_BASIC_H


template <typename Tint> struct UndoneOrbitInfo {
  size_t nbOrbitDone;
  Tint nbUndone;
  Face eSetUndone;
};

template <typename Tint>
UndoneOrbitInfo<Tint> get_default_undoneinfo(int n_rows) {
  Face f(n_rows);
  return {0, 0, f};
}




template <typename Tint>
UndoneOrbitInfo<Tint>
CombineUndoneOrbitInfo(const std::vector<UndoneOrbitInfo<Tint>> &LComb) {
  size_t nbOrbitDone = LComb[0].nbOrbitDone;
  Tint nbUndone = LComb[0].nbUndone;
  Face f = LComb[0].eSetUndone;
  for (size_t i = 1; i < LComb.size(); i++) {
    nbOrbitDone += LComb[i].nbOrbitDone;
    nbUndone += LComb[i].nbUndone;
    f &= LComb[i].eSetUndone;
  }
  return {nbOrbitDone, nbUndone, f};
}

template <typename Tint>
bool ComputeStatusUndone(const UndoneOrbitInfo<Tint> &eComb,
                         const Tint &CritSiz) {
  if (eComb.nbOrbitDone > 0)
    if (eComb.nbUndone <= CritSiz || eComb.eSetUndone.count() > 0)
      return true;
  return false;
}


// The condition on nbOrbitDone make the check more complex.
// For parallel, we use this monotonic partial check as heuristic
// about whether to do the major checks or not.
template <typename Tint>
bool MonotonicCheckStatusUndone(const UndoneOrbitInfo<Tint> &eComb,
                                const Tint &CritSiz) {
  if (eComb.nbUndone <= CritSiz || eComb.eSetUndone.count() > 0)
    return true;
  return false;
}



template <typename Tint> struct StatusUndoneOrbitInfo {
  bool status;
  UndoneOrbitInfo<Tint> erec;
};


namespace boost::serialization {

template <class Archive, typename Tint>
inline void serialize(Archive &ar, StatusUndoneOrbitInfo<Tint> &mesg,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("status", mesg.status);
  ar &make_nvp("nborbitdone", mesg.erec.nbOrbitDone);
  ar &make_nvp("nbundone", mesg.erec.nbUndone);
  ar &make_nvp("setundone", mesg.erec.eSetUndone);
}

} // namespace boost::serialization





#endif
