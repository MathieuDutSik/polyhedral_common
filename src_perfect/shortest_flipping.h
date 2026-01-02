// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_PERFECT_SHORTEST_FLIPPING_H_
#define SRC_PERFECT_SHORTEST_FLIPPING_H_

// clang-format off>
#include "EquiStabMemoization.h"
#include "Positivity.h"
#include "fractions.h"
#include <map>
#include <string>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_FLIP
#endif

#ifdef DISABLE_DEBUG_FLIP
#undef DEBUG_FLIP
#endif

#ifdef SANITY_CHECK
#define SANITY_CHECK_FLIP
#endif

#ifdef TIMINGS
#define TIMINGS_FLIP
#endif

// Use a set, better asymptotics, but copy needed
template <typename Tint>
bool TestInclusionSHV_set(MyMatrix<Tint> const &TheSHVbig,
                          MyMatrix<Tint> const &TheSHVsma) {
  std::unordered_set<MyVector<Tint>> set_big;
  int nbRowBig = TheSHVbig.rows();
  for (int iRowBig = 0; iRowBig < nbRowBig; iRowBig++) {
    MyVector<Tint> V = GetMatrixRow(TheSHVbig, iRowBig);
    set_big.insert(V);
  }
  int nbRowSma = TheSHVsma.rows();
  for (int iRowSma = 0; iRowSma < nbRowSma; iRowSma++) {
    MyVector<Tint> V = GetMatrixRow(TheSHVsma, iRowSma);
    if (set_big.count(V) == 0) {
      return false;
    }
  }
  return true;
}

// Use a quadratic algorithm. But not copy needed
template <typename Tint>
bool TestInclusionSHV_quad(MyMatrix<Tint> const &TheSHVbig,
                           MyMatrix<Tint> const &TheSHVsma) {
  int nbRowSma = TheSHVsma.rows();
  int nbRowBig = TheSHVbig.rows();
  int n = TheSHVsma.cols();
  auto is_equal_vector = [&](int iRowSma, int iRowBig) -> bool {
    for (int i = 0; i < n; i++) {
      Tint const &eVal = TheSHVbig(iRowBig, i);
      Tint const &fVal = TheSHVsma(iRowSma, i);
      if (eVal != fVal) {
        return false;
      }
    }
    return true;
  };
  auto is_small_vect_in_big = [&](int iRowSma) -> bool {
    for (int iRowBig = 0; iRowBig < nbRowBig; iRowBig++) {
      if (is_equal_vector(iRowSma, iRowBig)) {
        return true;
      }
    }
    return false;
  };
  for (int iRowSma = 0; iRowSma < nbRowSma; iRowSma++) {
    if (!is_small_vect_in_big(iRowSma)) {
      return false;
    }
  }
  return true;
}

template <typename Tint>
bool TestInclusionSHV(MyMatrix<Tint> const &TheSHVbig,
                      MyMatrix<Tint> const &TheSHVsma) {
  if (TheSHVbig.rows() > 100) {
    return TestInclusionSHV_set(TheSHVbig, TheSHVsma);
  } else {
    return TestInclusionSHV_quad(TheSHVbig, TheSHVsma);
  }
}



/*
  This is the code for given a form eMatIn, a direction of change eMatDir,
  to compute a coefficient lambda>0 such that eMatIn + lambda eMatDir has
  more shortest vector that eMatIn + epsilon eMatDir for epsilon infinitely
  small.
  ---
  Two functions need to be provided on input:
  + The f_admissible function that tells you whether you are in the admissible
    domain or not.
  + The f_shortest that returns the shortest vectors.
  ---
  There are two cases in which this routine is used:
  + The perfect forms enumeration
  + The copositive simplex algorithm
  See "A simplex algorithm for rational cp-factorization" arXiv:1807.01382
  ---
  Those flipping algorithms tend to be more intricate than one would expect.
  See the following for some example description:
  + Algorithm 2 (page 8) in "A simplex algorithm for rational cp-factorization"
  arXiv:1807.01382
  + Algorithm 2 (page 10) in "Enumerating perfect forms" arXiv:0901.1587
 */
template <typename T, typename Tint, typename Fadmissible, typename Fshortest>
std::pair<MyMatrix<T>, Tshortest<T, Tint>>
Kernel_Flipping_Perfect(Fadmissible f_admissible, Fshortest f_shortest,
                        MyMatrix<T> const &eMatIn, MyMatrix<T> const &eMatDir,
                        [[maybe_unused]] std::ostream &os) {
#ifdef SANITY_CHECK_FLIP
  if (f_admissible(eMatDir)) {
    std::cerr << "PERF: If the direction of change belong to the space, then "
                 "we are not going to find new vectors\n";
    throw TerminalException{1};
  }
#endif
  Memoization<Fadmissible, MyMatrix<T>, bool> memo_admissible(f_admissible);
  Memoization<Fshortest, MyMatrix<T>, Tshortest<T, Tint>> memo_shortest(
      f_shortest);
#ifdef DEBUG_FLIP
  os << "PERF: Kernel_Flipping_Perfect, case 1 (eMatIn)\n";
#endif
  Tshortest<T, Tint> const &rec_shv_in = memo_shortest.comp(eMatIn);
  //
  // Initial checks of the input
  //
#ifdef DEBUG_FLIP_DISABLE
  os << "PERF: Kernel_Flipping_Perfect : SHVinformation=\n";
  int nbSHV = rec_shv_in.SHV.rows();
  int n = rec_shv_in.SHV.cols();
  int nbZero_sumMat = 0;
  std::vector<int> ListScal(nbSHV);
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    MyVector<Tint> V = GetMatrixRow(rec_shv_in.SHV, iSHV);
    T sumMatIn = EvaluationQuadForm(eMatIn, V);
    T sumMatDir = EvaluationQuadForm(eMatDir, V);
    int eVal = 1;
    if (sumMatDir == 0) {
      nbZero_sumMat++;
      eVal = 0;
    }
    ListScal[iSHV] = eVal;
    os << "PERF: iSHV=" << iSHV << " V=";
    for (int i = 0; i < n; i++)
      os << V(i) << " ";
    os << "sumMatIn=" << sumMatIn << " sumMatDir=" << sumMatDir << "\n";
  }
  MyMatrix<Tint> SHVface(nbZero_sumMat, n);
  int idx = 0;
  for (int iSHV = 0; iSHV < nbSHV; iSHV++) {
    if (ListScal[iSHV] == 0) {
      for (int i = 0; i < n; i++)
        SHVface(idx, i) = rec_shv_in.SHV(iSHV, i);
      idx++;
    }
  }
  os << "PERF: SHVface=\n";
  WriteMatrix(os, SHVface);
  os << "PERF: nbZero_sumMat=" << nbZero_sumMat << "\n";
  MyMatrix<T> ConeClassicalInput =
      GetNakedPerfectConeClassical<T>(rec_shv_in.SHV);
  os << "PERF: RankMat(ConeClassicalInput)=" << RankMat(ConeClassicalInput)
     << "\n";
  MyMatrix<T> ConeClassicalFace = GetNakedPerfectConeClassical<T>(SHVface);
  os << "PERF: RankMat(ConeClassicalFace)=" << RankMat(ConeClassicalFace)
     << "\n";
#endif
  //
  // Running the algorithm
  // First loop where we find an interval where the solution may exist.
  //
  T bound_upp(1);
  T bound_low(0);
#ifdef DEBUG_FLIP
  int iterLoop_first = 0;
  int n_shortest = 0;
#endif
  while (true) {
    MyMatrix<T> Q_upp = eMatIn + bound_upp * eMatDir;
#ifdef DEBUG_FLIP
    iterLoop_first += 1;
    os << "PERF: iterLoop_first=" << iterLoop_first << " bound=" << bound_low
       << " | " << bound_upp << "\n";
#endif
#ifdef TIMINGS_PERFECT_FORM
    MicrosecondTime time_admissible;
#endif
    bool test = memo_admissible.comp(Q_upp);
#ifdef TIMINGS_PERFECT_FORM
    os << "|PERF: memo_admissible|=" << time_admissible << "\n";
#endif
    if (!test) {
      // We are outside, so reduce the range.
      bound_upp = get_mid_val(bound_low, bound_upp);
    } else {
#ifdef TIMINGS_PERFECT_FORM
      MicrosecondTime time_shortest;
#endif
      Tshortest<T, Tint> const &rec_shv_upp = memo_shortest.comp(Q_upp);
#ifdef TIMINGS_PERFECT_FORM
      os << "|ITER: memo_shortest|=" << time_shortest << "\n";
#endif
#ifdef DEBUG_FLIP
      n_shortest += 1;
#endif
#ifdef DEBUG_FLIP_DISABLE
      os << "ITER: rec_shv_upp.min=" << rec_shv_upp.min << "\n";
      os << "ITER: rec_shv_upp.SHV=\n";
      WriteMatrix(os, rec_shv_upp.SHV);
#endif
      if (rec_shv_upp.min == rec_shv_in.min) {
        // That kind of scheme is rather primitive.
        // However, the alternative is to write the equation det(xA + yB) = 0
        // and to look for rational solutions. A little risky.
        // For the Lorentzian stuff, things are easier as we have a quadratic
        // form to consider. But there is no reason to have a general solution.
        T nLow = bound_upp;
        T nUpp = 2 * bound_upp;
        bound_low = nLow;
        bound_upp = nUpp;
      } else {
        break;
      }
    }
  }
#ifdef SANITY_CHECK_FLIP
  {
    MyMatrix<T> Q_low = eMatIn + bound_low * eMatDir;
    MyMatrix<T> Q_upp = eMatIn + bound_upp * eMatDir;
    if (!memo_admissible.comp(Q_low)) {
      std::cerr << "PERF: We should have Qlow being admissible\n";
      throw TerminalException{1};
    }
    if (!memo_admissible.comp(Q_upp)) {
      std::cerr << "PERF: We should have Qlow being admissible\n";
      throw TerminalException{1};
    }
    Tshortest<T, Tint> const &rec_shv_upp = memo_shortest.comp(Q_upp);
    if (rec_shv_upp.min >= rec_shv_in.min) {
      std::cerr << "PERF: We should have rec_shv_upp.min < rec_shv_in.min\n";
      throw TerminalException{1};
    }
  }
#endif
  //
  // Second while loop
  //
  // Now we must have Q_low / Q_upp admissible
  // and rec_shv_upp.min < rec_shv_in.min
#ifdef DEBUG_FLIP
  os << "PERF: FIRST LOOP FINISHED bound_upp=" << bound_upp
     << " bound_low=" << bound_low << " iterLoop_first=" << iterLoop_first
     << " n_shortest=" << n_shortest << "\n";
  int iterLoop_second = 0;
#endif


  while (true) {
#ifdef DEBUG_FLIP
    os << "PERF: -------------- iterLoop_second=" << iterLoop_second
       << " bound=" << bound_low << " | " << bound_upp << " --------------\n";
    iterLoop_second += 1;
#endif
    MyMatrix<T> Q_low = eMatIn + bound_low * eMatDir;
    MyMatrix<T> Q_upp = eMatIn + bound_upp * eMatDir;
    Tshortest<T, Tint> const &rec_shv_low = memo_shortest.comp(Q_low);
    Tshortest<T, Tint> const &rec_shv_upp = memo_shortest.comp(Q_upp);
#ifdef DEBUG_FLIP_DISABLE
    os << "PERF: Kernel_Flipping_Perfect, case 4 SHV_low.min=" << rec_shv_low.min
       << " SHV_upp.min=" << rec_shv_upp.min << "\n";
    MyMatrix<T> ConeClassicalLow =
        GetNakedPerfectConeClassical<T>(rec_shv_low.SHV);
    os << "PERF: RankMat(ConeClassicalLow)=" << RankMat(ConeClassicalLow)
       << "\n";
    MyMatrix<T> ConeClassicalUpp =
        GetNakedPerfectConeClassical<T>(rec_shv_upp.SHV);
    os << "PERF: RankMat(ConeClassicalUpp)=" << RankMat(ConeClassicalUpp)
       << "\n";
#endif
    bool test1 = rec_shv_upp.min == rec_shv_in.min;
#ifdef SANITY_CHECK_FLIP
    if (rec_shv_upp.min > rec_shv_in.min) {
      std::cerr << "PERF: The shortest vectors should be of the same norm or "
                   "below the input\n";
      throw TerminalException{1};
    }
#endif
    bool test2 = TestInclusionSHV(rec_shv_in.SHV, rec_shv_low.SHV);
#ifdef DEBUG_FLIP_DISABLE
    os << "PERF: rec_shv_upp.min = " << rec_shv_upp.min << "\n";
    os << "PERF: rec_shv_in.min  = " << rec_shv_in.min << "\n";
    os << "PERF: test1=" << test1 << " test2=" << test2 << "\n";
#endif
    if (test1) {
#ifdef DEBUG_FLIP
      os << "PERF: Q_upp=\n";
      WriteMatrix(os, Q_upp);
#endif

#ifdef DEBUG_FLIP
      double coeff_d = UniversalScalarConversion<double, T>(bound_upp);
      os << "PERF: Exit flip (return Q_upp) iterLoop_second=" << iterLoop_second
         << " coeff=" << bound_upp << " coeff_d=" << coeff_d << "\n";
#endif
      return {std::move(Q_upp), std::move(rec_shv_upp)};
    }
    if (!test2 && test1) {
#ifdef DEBUG_FLIP_DISABLE
      os << "PERF: Qperf=\n";
      WriteMatrix(os, eMatIn);
      os << "PERF: Q_low=\n";
      WriteMatrix(os, Q_low);
      //
      os << "PERF: rec_shv_in.SHV=\n";
      WriteMatrix(os, rec_shv_in.SHV);
      os << "PERF: rec_shv_low.SHV=\n";
      WriteMatrix(os, rec_shv_low.SHV);
      os << "PERF: Return Q_low\n";
#endif
#ifdef DEBUG_FLIP
      os << "PERF: Exit flip (return Q_low) iterLoop_second=" << iterLoop_second
         << " coeff=" << bound_low << "\n";
#endif
      return {std::move(Q_low), std::move(rec_shv_low)};
    }
    T TheGamma = get_mid_val(bound_low, bound_upp);
    MyMatrix<T> Q_gamma = eMatIn + TheGamma * eMatDir;
#ifdef DEBUG_FLIP_DISABLE
    os << "PERF: Q_gamma=\n";
    WriteMatrix(os, Q_gamma);
#endif
#ifdef DEBUG_FLIP
    double TheGamma_d = UniversalScalarConversion<double, T>(TheGamma);
    os << "PERF: Kernel_Flipping_Perfect, gamma=" << TheGamma
       << " gamma_d=" << TheGamma_d << "\n";
    os << "PERF: Kernel_Flipping_Perfect, gamma=" << TheGamma
       << " bound_low=" << bound_low << " bound_upp=" << bound_upp << "\n";
#endif
    Tshortest<T, Tint> const &rec_shv_gamma = memo_shortest.comp(Q_gamma);
#ifdef DEBUG_FLIP_DISABLE
    os << "|rec_shv_gamma.SHV|=" << rec_shv_gamma.SHV.rows() << "\n";
    WriteMatrix(os, rec_shv_gamma.SHV);
    MyMatrix<T> ConeClassicalGamma =
        GetNakedPerfectConeClassical<T>(rec_shv_gamma.SHV);
    os << "PERF: RankMat(ConeClassicalGamma)=" << RankMat(ConeClassicalGamma)
       << "\n";
#endif
    if (rec_shv_gamma.min >= rec_shv_in.min) {
      bound_low = TheGamma;
    } else {
#ifdef DEBUG_FLIP
      os << "PERF: Assigning bound_upp to TheGamma=" << TheGamma << "\n";
#endif
      bound_upp = TheGamma;
      int nbRow = rec_shv_gamma.SHV.rows();
      for (int iRow = 0; iRow < nbRow; iRow++) {
        MyVector<Tint> V = GetMatrixRow(rec_shv_gamma.SHV, iRow);
        T rVal = EvaluationQuadForm<T, Tint>(eMatDir, V);
#ifdef DEBUG_FLIP
        os << "PERF: iRow=" << iRow << " / " << nbRow
           << " V=" << StringVectorGAP(V) << "\n";
#endif
        if (rVal < 0) {
          T qVal = EvaluationQuadForm<T, Tint>(eMatIn, V);
          T TheVal = (rec_shv_in.min - qVal) / rVal;
          if (TheVal < bound_upp) {
#ifdef DEBUG_FLIP
            os << "PERF: iRow=" << iRow
               << " Assigning bound_upp to TheVal=" << TheVal << "\n";
#endif
            bound_upp = TheVal;
          }
        }
      }
    }
  }
}

// clang-format off
#endif  // SRC_PERFECT_SHORTEST_FLIPPING_H_
// clang-format on
