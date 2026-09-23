// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "InvariantVectorFamily.h"
#include "LatticePleskenSouvignier.h"
// clang-format on

template <typename T>
void process(std::string choice, std::string MatFile,
             std::string const &OutFormat, std::string const &OutFile) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  MyMatrix<T> GramMat = ReadMatrixFile<T>(MatFile);
  if (!IsSymmetricMatrix(GramMat) ||
      !IsPositiveDefinite(GramMat, std::cerr)) {
    std::cerr << "LATT_GenerateCharacteristicVectorSet: The input Gram matrix in "
              << MatFile << " is not symmetric positive definite\n";
    throw TerminalException{1};
  }
  // Whether the returned family is closed under v -> -v. The "half"
  // families hold one vector per antipodal pair and so are not, which the
  // antipodality sanity check below must skip. Only read under SANITY_CHECK.
  [[maybe_unused]] bool is_antipodal = true;
  auto f = [&]() -> MyMatrix<Tint> {
    if (choice == "shortest") {
      Tshortest<T, Tint> rec = T_ShortestVector<T, Tint>(GramMat, std::cerr);
      return rec.SHV;
    }
    if (choice == "iterated_shortest") {
      return IteratedShortestVectorFamily<T, Tint>(GramMat, std::cerr);
    }
    if (choice == "iterated_shortest_half") {
      is_antipodal = false;
      return IteratedShortestVectorFamilyHalf<T, Tint>(GramMat, std::cerr);
    }
    // The iterated-shortest family completed to a Z-spanning one by the
    // generic saturation completion ivf_z_spanning.
    if (choice == "span_iterated_shortest") {
      MyMatrix<Tint> V =
          IteratedShortestVectorFamily<T, Tint>(GramMat, std::cerr);
      return ivf_z_spanning<T, Tint>(GramMat, V, std::cerr);
    }
    // The Z-spanning characteristic vector family built through the
    // orthogonal recursion (Hecke's _characteristic_vectors), keeping every
    // closest vector problem on the shortest-vector sublattice.
    if (choice == "inner_span_iterated_shortest") {
      return inner_span_iterated_shortest<T, Tint>(GramMat, std::cerr);
    }
    // The root variants: the tower starts from the roots (norm-2 vectors)
    // rather than the shortest vectors. Full-rank and its Z-spanning form.
    if (choice == "root_iterated_shortest") {
      return root_iterated_shortest<T, Tint>(GramMat, std::cerr);
    }
    if (choice == "inner_span_root_iterated_shortest") {
      return inner_span_root_iterated_shortest<T, Tint>(GramMat, std::cerr);
    }
    if (choice == "relevant_voronoi") {
      return ComputeVoronoiRelevantVector<T, Tint>(GramMat, std::cerr);
    }
    if (choice == "filtered_relevant_voronoi") {
      MyMatrix<Tint> M =
          ComputeVoronoiRelevantVector<T, Tint>(GramMat, std::cerr);
      return FilterByNorm(GramMat, M, std::cerr);
    }
    if (choice == "fullrank") {
      return ExtractInvariantVectorFamilyFullRank<T, Tint>(GramMat, std::cerr);
    }
    if (choice == "spanning") {
      return ExtractInvariantVectorFamilyZbasis<T, Tint>(GramMat, std::cerr);
    }
    // No limit here: the command line asks for the set itself, not for the
    // cheaper of two families, so the construction always runs to the end.
    CharVectSet_Budget budget;
    if (choice == "wr_cv") {
      return *CharacteristicVectorSetWellRoundedCV<T, Tint>(GramMat, true,
                                                            false, budget,
                                                            std::cerr);
    }
    if (choice == "cv") {
      return *CharacteristicVectorSetCV<T, Tint>(GramMat, true, false, budget,
                                                 std::cerr);
    }
    if (choice == "cv_fullrank") {
      return *CharacteristicVectorSetCV<T, Tint>(GramMat, false, false, budget,
                                                 std::cerr);
    }
    // The canonical vector family used by the graph-based canonicalization
    // (small even on the root-rich classes), whole and halved.
    if (choice == "canonic") {
      return GetCanonicVectorFamily<T, Tint>(GramMat, std::cerr).get_full();
    }
    if (choice == "canonic_half") {
      is_antipodal = false;
      return GetCanonicVectorFamily<T, Tint>(GramMat, std::cerr).SHVhalf;
    }
    // The Plesken-Souvignier norm-ball family: the vectors of norm at most
    // the largest diagonal entry of the LLL-reduced Gram, expressed in the
    // original basis. This is the family the PS engine searches, and it
    // explodes on the classes whose reduced diagonal has one large entry.
    if (choice == "plesken_souvignier" || choice == "ps_normball") {
      LLLreduction<T, Tint> rec = LLLreducedBasis<T, Tint>(GramMat, std::cerr);
      T bound = PleskenSouvignierBound(rec.GramMatRed);
      MyMatrix<Tint> half =
          PleskenSouvignierVectorFamily<T, Tint>(rec.GramMatRed, bound,
                                                 std::cerr);
      int n_pair = half.rows();
      int dim = half.cols();
      MyMatrix<Tint> res(2 * n_pair, dim);
      for (int i = 0; i < n_pair; i++) {
        MyVector<Tint> V = GetMatrixRow(half, i);
        MyVector<Tint> Vorig = rec.Pmat.transpose() * V;
        for (int k = 0; k < dim; k++) {
          res(i, k) = Vorig(k);
          res(n_pair + i, k) = -Vorig(k);
        }
      }
      return res;
    }
    // The half (one per antipodal pair) versions of the full-rank and
    // spanning invariant families.
    if (choice == "fullrank_half") {
      is_antipodal = false;
      return ExtractInvariantVectorFamilyFullRankHalf<T, Tint>(GramMat,
                                                               std::cerr);
    }
    if (choice == "spanning_half") {
      is_antipodal = false;
      return ExtractInvariantVectorFamilyZbasisHalf<T, Tint>(GramMat,
                                                             std::cerr);
    }
    std::cerr << "Failed to find a matching entry for choice\n";
    std::cerr << "Possible choices: shortest, iterated_shortest, "
                 "iterated_shortest_half, span_iterated_shortest, inner_span_iterated_shortest, root_iterated_shortest, inner_span_root_iterated_shortest, relevant_voronoi, "
                 "filtered_relevant_voronoi, fullrank, fullrank_half, "
                 "spanning, spanning_half, wr_cv, cv, cv_fullrank, canonic, "
                 "canonic_half, plesken_souvignier\n";
    throw TerminalException{1};
  };
  MyMatrix<Tint> M = f();
#ifdef SANITY_CHECK
  if (is_antipodal) {
    check_antipodality_mymatrix(M);
  }
#endif
  auto f_print = [&](std::ostream &osf) -> void {
    if (OutFormat == "structure") {
      int n_vect = M.rows();
      int n = GramMat.cols();
      int rank = RankMat(M);
      // Antipodality: whether the family is closed under v -> -v.
      bool antipodal = true;
      {
        std::unordered_set<MyVector<Tint>> Sset;
        for (int i_vect = 0; i_vect < n_vect; i_vect++) {
          Sset.insert(GetMatrixRow(M, i_vect));
        }
        for (int i_vect = 0; i_vect < n_vect; i_vect++) {
          MyVector<Tint> negv = -GetMatrixRow(M, i_vect);
          if (Sset.find(negv) == Sset.end()) {
            antipodal = false;
            break;
          }
        }
      }
      // Index of L = <rows of M> in its saturation (L (x) R) cap Z^n. With
      // B a Z-basis of L and Bsat one of the saturation, B = C Bsat for an
      // integer C, and the index is |det C| = |det(B Bsat^T)| /
      // det(Bsat Bsat^T). Works whatever the rank: for a full-rank family
      // the saturation is Z^n and this is the index of L in Z^n.
      Tint saturation_index(1);
      if (n_vect > 0) {
        MyMatrix<Tint> B = GetZbasis(M);
        MyMatrix<Tint> Bsat = IntegralSpaceSaturation(B);
        MyMatrix<Tint> Prod1 = B * Bsat.transpose();
        MyMatrix<Tint> Prod2 = Bsat * Bsat.transpose();
        Tint num = T_abs(DeterminantMatBareiss(Prod1));
        Tint den = DeterminantMatBareiss(Prod2);
        saturation_index = num / den;
      }
      std::map<T, size_t> map;
      for (int i_vect = 0; i_vect < n_vect; i_vect++) {
        MyVector<Tint> V = GetMatrixRow(M, i_vect);
        T norm = EvaluationQuadForm<T, Tint>(GramMat, V);
        map[norm] += 1;
      }
      osf << "|M| = " << n_vect << " / " << n << "\n";
      osf << "rank = " << rank << "\n";
      osf << "antipodal = " << (antipodal ? "true" : "false") << "\n";
      osf << "saturation_index = " << saturation_index << "\n";
      osf << "norms =";
      for (auto &[norm, multiplicity] : map) {
        osf << " [" << norm << " : " << multiplicity << " ]";
      }
      osf << "\n";
      return;
    }
    if (OutFormat == "GAP") {
      osf << "return ";
      WriteMatrixGAP(osf, M);
      osf << ";\n";
      return;
    }
    if (OutFormat == "CPP") {
      return WriteMatrix(osf, M);
    }
    std::cerr << "Failed to find a matching entry for OutFormat\n";
    std::cerr << "Allowed choices: structure, GAP, CPP\n";
    throw TerminalException{1};
  };
  FILE_PrintStderrStdoutFile(OutFile, f_print);
}

int main(int argc, char *argv[]) {
  maybe_install_gmp_pool();
  HumanTime time;
  try {
    if (argc != 6 && argc != 4) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "LATT_GenerateCharacteristicVectorSet [arith] choice "
                << "[MatFile] [OutFormat] [OutFile]\n";
      std::cerr << "       or\n";
      std::cerr << "LATT_GenerateCharacteristicVectorSet [arith] choice "
                << "[MatFile]\n";
      std::cerr << "allowed choices:\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "[arith]: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "[arith]: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      std::cerr << "choice: shortest, iterated_shortest, iterated_shortest_half, span_iterated_shortest, relevant_voronoi, "
                   "filtered_relevant_voronoi, fullrank, fullrank_half, "
                   "spanning, spanning_half, wr_cv, cv, cv_fullrank, canonic, "
                   "canonic_half, plesken_souvignier\n";
      std::cerr << "OutFormat: structure, GAP, CPP\n";
      std::cerr << "OutFile: stderr, stdout, my_file\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string choice = argv[2];
    std::string MatFile = argv[3];
    std::string OutFormat = "structure";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    auto f=[&]() -> void {
      if (arith == "gmp") {
        using T = mpq_class;
        return process<T>(choice, MatFile, OutFormat, OutFile);
      }
#ifdef ENABLE_BOOST_TYPES
      if (arith == "gmp_boost") {
        using T = boost::multiprecision::mpq_rational;
        return process<T>(choice, MatFile, OutFormat, OutFile);
      }
      if (arith == "multi_boost") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(choice, MatFile, OutFormat, OutFile);
      }
#endif
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint") {
        using T = fmpq_class;
        return process<T>(choice, MatFile, OutFormat, OutFile);
      }
#endif
      std::cerr << "process_A failure: No matching entry for arith\n";
#ifdef ENABLE_FLINT_SUPPORT
      std::cerr << "Allowed values: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << ", flint\n";
#else
      std::cerr << "Allowed values: gmp";
#ifdef ENABLE_BOOST_TYPES
      std::cerr << ", gmp_boost, multi_boost";
#endif
      std::cerr << " (build with ENABLE_FLINT_SUPPORT=1 for flint)\n";
#endif
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of LATT_GenerateCharacteristicVectorSet\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in LATT_GenerateCharacteristicVectorSet\n";
    exit(e.eVal);
  }
  runtime(time);
}
