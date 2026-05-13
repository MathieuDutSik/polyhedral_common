// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPNEST_H_
#define SRC_GROUP_MATRIXGROUPNEST_H_

// clang-format off
#include "MatrixGroupBasic.h"
#include "Group.h"
#include "Permutation.h"
#include <optional>
#include <utility>
#include <vector>
// clang-format on

/*
  Thin wrappers around permutalib primitives that take TeltMatr =
  MyMatrixContainer<T>.  Wrapping MyMatrix<T> inside MyMatrixContainer<T>
  (which lives in the global namespace) keeps ADL happy when permutalib
  templates call Inverse / IsIdentity / operator< / operator== on the
  matrix type, regardless of the namespace in which T itself lives
  (e.g. boost::multiprecision::mpq_rational).
 */

template <typename T, typename Tgroup>
std::vector<MyMatrix<T>>
PreImageSubgroupContainer(std::vector<MyMatrix<T>> const &ListMatrGens,
                          std::vector<typename Tgroup::Telt> const &ListPermGens,
                          MyMatrix<T> const &id_matr, Tgroup const &stab,
                          [[maybe_unused]] std::ostream &os) {
  std::vector<MyMatrixContainer<T>> ListMatrGens_ct =
      get_vector_mmc(ListMatrGens);
  MyMatrixContainer<T> id_matr_ct(id_matr);
  std::vector<MyMatrixContainer<T>> l_ret_ct =
      permutalib::PreImageSubgroup<Tgroup, MyMatrixContainer<T>>(
          ListMatrGens_ct, ListPermGens, id_matr_ct, stab);
  std::vector<MyMatrix<T>> l_ret;
  for (auto &eMatr : l_ret_ct) {
    l_ret.push_back(eMatr.get_m());
  }
  return l_ret;
}

template <typename T, typename Tgroup>
std::optional<MyMatrix<T>>
RepresentativeActionMatrixPermSubsetContainer(
    std::vector<MyMatrix<T>> const &ListMatr,
    std::vector<typename Tgroup::Telt> const &ListPerm,
    MyMatrix<T> const &id_matr, Face const &eFace1, Face const &eFace2) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  std::vector<MyMatrixContainer<T>> ListMatr_ct = get_vector_mmc(ListMatr);
  MyMatrixContainer<T> id_matr_ct(id_matr);
  std::optional<MyMatrixContainer<T>> opt =
      permutalib::RepresentativeActionMatrixPermSubset<
          Telt, MyMatrixContainer<T>, TintGroup>(ListMatr_ct, ListPerm,
                                                  id_matr_ct, eFace1, eFace2);
  if (opt) {
    return opt->get_const_m();
  }
  return {};
}

template <typename T, typename Tgroup>
std::pair<std::vector<MyMatrix<T>>, std::vector<MyMatrix<T>>>
StabilizerRightCosetMatrixPermSubsetContainer(
    std::vector<MyMatrix<T>> const &ListMatr,
    std::vector<typename Tgroup::Telt> const &ListPerm,
    MyMatrix<T> const &id_matr, Face const &eFace) {
  using Telt = typename Tgroup::Telt;
  using TintGroup = typename Tgroup::Tint;
  std::vector<MyMatrixContainer<T>> ListMatr_ct = get_vector_mmc(ListMatr);
  MyMatrixContainer<T> id_matr_ct(id_matr);
  std::pair<std::vector<MyMatrixContainer<T>>,
            std::vector<MyMatrixContainer<T>>>
      pair = permutalib::StabilizerRightCosetMatrixPermSubset<
          Telt, MyMatrixContainer<T>, TintGroup>(ListMatr_ct, ListPerm,
                                                 id_matr_ct, eFace);
  std::vector<MyMatrix<T>> list1;
  for (auto &M : pair.first) {
    list1.push_back(M.get_const_m());
  }
  std::vector<MyMatrix<T>> list2;
  for (auto &M : pair.second) {
    list2.push_back(M.get_const_m());
  }
  return {std::move(list1), std::move(list2)};
}

/*
  Wrapper around permutalib::PreImagerElement<Telt, MyMatrixContainer<T>,
  TintGroup>.  Exposes a MyMatrix<T>-flavored interface so callers do not
  have to know about MyMatrixContainer<T> or include the permutalib
  headers directly.
 */
template <typename T, typename Telt, typename TintGroup>
class PreImagerElementContainer {
private:
  permutalib::PreImagerElement<Telt, MyMatrixContainer<T>, TintGroup> inner;

  static std::vector<MyMatrixContainer<T>>
  to_container(std::vector<MyMatrix<T>> const &l_matr) {
    return get_vector_mmc(l_matr);
  }

public:
  PreImagerElementContainer(std::vector<MyMatrix<T>> const &l_matr,
                            std::vector<Telt> const &l_perm,
                            MyMatrix<T> const &id_matr)
      : inner(to_container(l_matr), l_perm, MyMatrixContainer<T>(id_matr)) {}

  std::optional<MyMatrix<T>> get_preimage(Telt const &elt) const {
    std::optional<MyMatrixContainer<T>> opt = inner.get_preimage(elt);
    if (!opt) {
      return {};
    }
    return opt->get_const_m();
  }

  std::vector<MyMatrix<T>> get_list_matr_gens() const {
    std::vector<MyMatrixContainer<T>> const &l_ct = inner.get_list_matr_gens();
    std::vector<MyMatrix<T>> l_ret;
    l_ret.reserve(l_ct.size());
    for (auto const &mc : l_ct) {
      l_ret.push_back(mc.get_const_m());
    }
    return l_ret;
  }
};

// clang-format off
#endif  //  SRC_GROUP_MATRIXGROUPNEST_H_
// clang-format on
