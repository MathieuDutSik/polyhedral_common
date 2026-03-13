// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_K_COVERING_CIRCUMSCRIBING_SPHERE_H_
#define SRC_K_COVERING_CIRCUMSCRIBING_SPHERE_H_

// clang-format off
#include "MAT_Matrix.h"
#include <algorithm>
#include <numeric>
#include <optional>
#include <random>
#include <vector>
// clang-format on

template <typename T> struct CircumscribingSphere {
  MyVector<T> Center;
  T SquareRadius;
  std::vector<int> Support;
};

namespace circumscribing_sphere_impl {

template <typename T> T MetricScalarProduct(MyMatrix<T> const &GramMat,
                                            MyVector<T> const &V1,
                                            MyVector<T> const &V2) {
  T eSum = 0;
  int dim = GramMat.rows();
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      eSum += V1(i) * GramMat(i, j) * V2(j);
  return eSum;
}

template <typename T> T SquaredNorm(MyMatrix<T> const &GramMat,
                                    MyVector<T> const &V) {
  return MetricScalarProduct(GramMat, V, V);
}

template <typename T> T SquaredDistance(MyMatrix<T> const &GramMat,
                                        MyVector<T> const &V1,
                                        MyVector<T> const &V2) {
  MyVector<T> eDiff(V1.size());
  for (int i = 0; i < V1.size(); i++) {
    eDiff(i) = V1(i) - V2(i);
  }
  return SquaredNorm(GramMat, eDiff);
}

template <typename T> MyVector<T> GetPoint(MyMatrix<T> const &PointSet, int iRow) {
  return PointSet.row(iRow).transpose();
}

template <typename T>
bool Contains(MyMatrix<T> const &GramMat, CircumscribingSphere<T> const &Sphere,
              MyVector<T> const &Point) {
  T eDist = SquaredDistance(GramMat, Sphere.Center, Point);
  return eDist <= Sphere.SquareRadius;
}

template <typename T>
std::optional<CircumscribingSphere<T>>
SphereFromIndependentSupport(MyMatrix<T> const &GramMat,
                             MyMatrix<T> const &PointSet,
                             std::vector<int> const &Support) {
  int dim = PointSet.cols();
  if (Support.empty()) {
    return CircumscribingSphere<T>{ZeroVector<T>(dim), 0, {}};
  }
  if (Support.size() == 1) {
    MyVector<T> ePoint = GetPoint(PointSet, Support[0]);
    return CircumscribingSphere<T>{std::move(ePoint), 0, Support};
  }
  int nbVect = static_cast<int>(Support.size()) - 1;
  MyVector<T> p0 = GetPoint(PointSet, Support[0]);
  MyMatrix<T> U(nbVect, dim);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    MyVector<T> p = GetPoint(PointSet, Support[iVect + 1]);
    for (int i = 0; i < dim; i++)
      U(iVect, i) = p(i) - p0(i);
  }
  if (RankMat(U) != nbVect)
    return {};
  MyMatrix<T> Gram = U * GramMat * U.transpose();
  MyVector<T> B(nbVect);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    MyVector<T> eRow = U.row(iVect).transpose();
    B(iVect) = SquaredNorm(GramMat, eRow) / 2;
  }
  std::optional<MyVector<T>> optLambda = SolutionMat(Gram, B);
  if (!optLambda)
    return {};
  MyVector<T> Center = p0;
  for (int iVect = 0; iVect < nbVect; iVect++)
    for (int i = 0; i < dim; i++)
      Center(i) += (*optLambda)(iVect) * U(iVect, i);
  T SquareRadius = SquaredDistance(GramMat, Center, p0);
  return CircumscribingSphere<T>{std::move(Center), std::move(SquareRadius),
                                 Support};
}

template <typename T>
CircumscribingSphere<T>
SphereFromBoundary(MyMatrix<T> const &GramMat, MyMatrix<T> const &PointSet,
                   std::vector<int> const &Boundary) {
  int nbBound = static_cast<int>(Boundary.size());
  std::optional<CircumscribingSphere<T>> Best;
  int nbCase = 1 << nbBound;
  for (int iCase = 0; iCase < nbCase; iCase++) {
    std::vector<int> Support;
    for (int i = 0; i < nbBound; i++)
      if ((iCase >> i) & 1)
        Support.push_back(Boundary[i]);
    std::optional<CircumscribingSphere<T>> optSphere =
        SphereFromIndependentSupport(GramMat, PointSet, Support);
    if (!optSphere)
      continue;
    bool test = true;
    for (int idx : Boundary)
      if (!Contains(GramMat, *optSphere, GetPoint(PointSet, idx))) {
        test = false;
        break;
      }
    if (!test)
      continue;
    if (!Best || optSphere->SquareRadius < Best->SquareRadius ||
        (optSphere->SquareRadius == Best->SquareRadius &&
         optSphere->Support.size() < Best->Support.size()))
      Best = std::move(optSphere);
  }
  if (Best)
    return *Best;
  return CircumscribingSphere<T>{ZeroVector<T>(PointSet.cols()), 0, {}};
}

template <typename T>
CircumscribingSphere<T> WelzlRecursive(MyMatrix<T> const &GramMat,
                                       MyMatrix<T> const &PointSet,
                                       std::vector<int> const &Order,
                                       int nbPoint,
                                       std::vector<int> &Boundary) {
  if (nbPoint == 0 || Boundary.size() == static_cast<size_t>(PointSet.cols() + 1))
    return SphereFromBoundary(GramMat, PointSet, Boundary);
  int idx = Order[nbPoint - 1];
  CircumscribingSphere<T> Sphere =
      WelzlRecursive(GramMat, PointSet, Order, nbPoint - 1, Boundary);
  if (Contains(GramMat, Sphere, GetPoint(PointSet, idx)))
    return Sphere;
  Boundary.push_back(idx);
  CircumscribingSphere<T> NewSphere =
      WelzlRecursive(GramMat, PointSet, Order, nbPoint - 1, Boundary);
  Boundary.pop_back();
  return NewSphere;
}

template <typename T> bool HasHomogeneousCoordinate(MyMatrix<T> const &PointSet) {
  if (PointSet.cols() == 0)
    return false;
  if (PointSet.rows() == 0)
    return false;
  for (int iRow = 0; iRow < PointSet.rows(); iRow++)
    if (PointSet(iRow, 0) != 1)
      return false;
  return true;
}

template <typename T>
MyMatrix<T> ExtractAffinePointSet(MyMatrix<T> const &PointSet,
                                  bool FirstColumnIsHomogeneous) {
  // In the common EXT convention the first column is the homogenizing 1.
  if (!FirstColumnIsHomogeneous)
    return PointSet;
  if (PointSet.cols() == 0)
    return PointSet;
  int nbPoint = PointSet.rows();
  int dim = PointSet.cols() - 1;
  MyMatrix<T> RetMat(nbPoint, dim);
  for (int iRow = 0; iRow < nbPoint; iRow++)
    for (int i = 0; i < dim; i++)
      RetMat(iRow, i) = PointSet(iRow, i + 1);
  return RetMat;
}

template <typename Tfield>
CircumscribingSphere<Tfield>
ComputeCircumscribingSphereField(MyMatrix<Tfield> const &GramMat,
                                 MyMatrix<Tfield> const &PointSet) {
  int nbPoint = PointSet.rows();
  std::vector<int> Order(nbPoint);
  std::iota(Order.begin(), Order.end(), 0);
  std::mt19937 rng(static_cast<uint32_t>(PointSet.rows() * 257 + PointSet.cols()));
  std::shuffle(Order.begin(), Order.end(), rng);
  std::vector<int> Boundary;
  return WelzlRecursive(GramMat, PointSet, Order, nbPoint, Boundary);
}

}  // namespace circumscribing_sphere_impl

template <typename T>
CircumscribingSphere<typename overlying_field<T>::field_type>
ComputeCircumscribingSphere(MyMatrix<T> const &GramMat,
                            MyMatrix<T> const &PointSet,
                            bool FirstColumnIsHomogeneous = false) {
  using Tfield = typename overlying_field<T>::field_type;
  // This follows Welzl's randomized incremental scheme, with a deterministic
  // shuffle so the output is reproducible.
  bool test_homogeneous =
      FirstColumnIsHomogeneous ||
      circumscribing_sphere_impl::HasHomogeneousCoordinate(PointSet);
  MyMatrix<T> PointSetAff =
      circumscribing_sphere_impl::ExtractAffinePointSet(PointSet,
                                                        test_homogeneous);
  MyMatrix<Tfield> GramMatField =
      UniversalMatrixConversion<Tfield, T>(GramMat);
  MyMatrix<Tfield> PointSetField =
      UniversalMatrixConversion<Tfield, T>(PointSetAff);
  return circumscribing_sphere_impl::ComputeCircumscribingSphereField(
      GramMatField, PointSetField);
}

template <typename T>
CircumscribingSphere<typename overlying_field<T>::field_type>
ComputeCircumscribingSphereEXT(MyMatrix<T> const &GramMat,
                               MyMatrix<T> const &EXT) {
  return ComputeCircumscribingSphere(GramMat, EXT, true);
}

#endif  // SRC_K_COVERING_CIRCUMSCRIBING_SPHERE_H_
