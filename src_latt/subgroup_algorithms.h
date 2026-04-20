// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_LATT_SUBGROUP_ALGORITHMS_H_
#define SRC_LATT_SUBGROUP_ALGORITHMS_H_

#ifdef DEBUG
#define DEBUG_SUBGROUP_ALGORITHM
#endif


template<typename Tgroup, typename Fcorrect>
std::pair<std::vector<typename Tgroup::Telt>, Tgroup> get_intermediate_group(typename Tgroup::Telt::Tidx const& n_act,
                                                                             std::vector<typename Tgroup::Telt> const& LGenSma,
                                                                             std::vector<typename Tgroup::Telt> const& LGenBig,
                                                                             Fcorrect f_correct,
                                                                             std::ostream& os) {
  using Telt = typename Tgroup::Telt;
  using LeftCosets = typename Tgroup::LeftCosets;
  std::vector<Telt> LGenSma_work = LGenSma;
  Tgroup GRPbig(LGenBig, n_act);
  Tgroup GRPsma(LGenSma_work, n_act);
#ifdef DEBUG_SUBGROUP_ALGORITHM
  os << "SUBA: get_intermediate_group |GRPbig|=" << GRPbig.size() << " |GRPsma|=" << GRPsma.size() << "\n";
#endif
  auto try_upgrade = [&]() -> std::optional<Telt> {
    // We can use either the left or right cosets.
    // This is because the right thing to use is the double cosets.
    // However, we do not have the formalism for having iterator
    // over the double cosets. We build all of them.
    //
    // For the left/right cosets we have efficient iterators
    // and that is why we use them here. We cannot afford at all
    // to enumerate all the double cosets because we will have
    // some scenario where FullGRP is big, GRPsub very small and
    // that would mean enumerating all the elements of the group.
    //
    // Left  transversals are g H
    // Right transversals are H g
    LeftCosets rc = GRPbig.left_cosets(GRPsma);
    for (auto &eCosReprPerm : rc) {
      bool test = f_correct(eCosReprPerm);
      if (test) {
        // We have this problem that the first cosets is not necessarily the one
        // of GRPsub and that the coset of the GRPsub is also not necessarily
        // the identity.
        if (!GRPsma.isin(eCosReprPerm)) {
#ifdef DEBUG_SUBGROUP_ALGORITHM
          os << "SUBA: get_intermediate_group Finding a new eCosReprPerm\n";
#endif
          return eCosReprPerm;
        }
      }
    }
    return {};
  };
  while (true) {
    std::optional<Telt> opt = try_upgrade();
    if (opt) {
      // Found another stabilizing element, upgrading the group and retry.
      LGenSma_work.push_back(*opt);
      GRPsma = Tgroup(LGenSma_work, n_act);
#ifdef DEBUG_SUBGROUP_ALGORITHM
      os << "SUBA: get_intermediate_group Now |GRPsub|=" << GRPsma.size() << "\n";
#endif
    } else {
      break;
    }
  }
  std::pair<std::vector<Telt>, Tgroup> pair{std::move(LGenSma_work), std::move(GRPsma)};
  return pair;
}



// clang-format off
#endif  // SRC_LATT_SUBGROUP_ALGORITHMS_H_
// clang-format on
