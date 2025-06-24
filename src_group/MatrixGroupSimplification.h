// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GROUP_MATRIXGROUPSIMPLIFICATION_H_
#define SRC_GROUP_MATRIXGROUPSIMPLIFICATION_H_

#include "MAT_Matrix.h"


#ifdef DEBUG
#define DEBUG_MATRIX_GROUP_SIMPLIFICATION
#define DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
#endif

#ifdef TRACK_INFO
#define TRACK_INFO_MATRIX_GROUP_SIMPLIFICATION
#endif

//#define DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE


template<typename T>
struct ComplexityMeasure {
  T ell1;
  T ellinfinity;
};

template<typename T>
ComplexityMeasure<T> get_complexity_measure(MyMatrix<T> const& M) {
  int n = M.rows();
  T ell1(0);
  T ellinfinity(0);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      T val = M(i,j);
      T abs_val = T_abs(val);
      ell1 += abs_val;
      if (abs_val > ellinfinity) {
        ellinfinity = abs_val;
      }
    }
  }
  return {ell1, ellinfinity};
}

template<typename T>
T get_ell1_complexity_measure(MyMatrix<T> const& M) {
  int n = M.rows();
  T ell1(0);
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      T val = M(i,j);
      T abs_val = T_abs(val);
      ell1 += abs_val;
    }
  }
  return ell1;
}

// The tool for simplifying a list of generators.
// The transformations being applied are Tietze transformations.
// We could add some conjugacy operations like  U V U^{-1}.
// But those are more expensive
template<typename Tnorm, typename Ttype, typename Fcomplexity, typename Finvers, typename Fproduct>
std::vector<Ttype> ExhaustiveReductionComplexityKernel(std::vector<Ttype> const& ListM, Fcomplexity f_complexity, Finvers f_invers, Fproduct f_product, [[maybe_unused]] std::ostream& os) {
  size_t miss_val = std::numeric_limits<size_t>::max();
  using TtypePair = std::pair<Ttype, Ttype>;
  using TcombPair = std::pair<TtypePair, Tnorm>;
  using Tcomb = std::pair<Ttype, Tnorm>;
  auto f_comp=[](TcombPair const& a, TcombPair const& b) -> bool {
    if (a.second < b.second) {
      return true;
    }
    if (a.second > b.second) {
      return false;
    }
    return a.first.first < b.first.first;
  };
  auto get_comb=[&](Ttype const& eM) -> Tcomb {
    Tnorm comp = f_complexity(eM);
    return {eM, comp};
  };
  auto get_comb_pair=[&](Tcomb const& p) -> TcombPair {
    Ttype p_inv = f_invers(p.first);
    TtypePair p_pair{p.first, p_inv};
    return {p_pair, p.second};
  };
  std::map<TcombPair, size_t, decltype(f_comp)> map(f_comp);
  size_t nonce = 0;
  for (auto & eM: ListM) {
    Tcomb comb1 = get_comb(eM);
    TcombPair comb2 = get_comb_pair(comb1);
    map[comb2] = nonce;
    nonce += 1;
  }
  std::unordered_set<std::pair<size_t,size_t>> set_treated;
  // Generate the possible ways to simplify the pair of elements.
  auto f_generate_candidate=[&](TcombPair const& a, TcombPair const& b) -> std::vector<Tcomb> {
    Ttype const& a_dir = a.first.first;
    Ttype const& b_dir = b.first.first;
    Ttype const& a_inv = a.first.second;
    Ttype const& b_inv = b.first.second;
    //
    Ttype prod1 = f_product(a_dir, b_dir);
    Tcomb pair1 = get_comb(prod1);
    //
    Ttype prod2 = f_product(a_inv, b_dir);
    Tcomb pair2 = get_comb(prod2);
    //
    Ttype prod3 = f_product(a_dir, b_inv);
    Tcomb pair3 = get_comb(prod3);
    //
    Ttype prod4 = f_product(a_inv, b_inv);
    Tcomb pair4 = get_comb(prod4);
    return {pair1, pair2, pair3, pair4};
  };
  // Selects the best candidates in the 4 being generated.
  auto f_get_best_candidate=[&](TcombPair const& a, TcombPair const& b) -> Tcomb {
    std::vector<Tcomb> l_comb = f_generate_candidate(a, b);
    Tcomb best_comp = l_comb[0];
    for (size_t i=1; i<l_comb.size(); i++) {
      if (l_comb[i].second < best_comp.second) {
        best_comp = l_comb[i];
      }
    }
    return best_comp;
  };
  // Iterate the reduction algorithm over pairs of elements.
  // The result of the iteration might be a 0, 1 or 2 new elements.
  auto f_reduce=[&](TcombPair const& a, TcombPair const& b) -> std::pair<size_t, std::vector<TcombPair>> {
    TcombPair a_work = a;
    TcombPair b_work = b;
    size_t n_change = 0;
    while(true) {
      Tcomb cand = f_get_best_candidate(a_work, b_work);
      Tnorm a_norm = a_work.second;
      Tnorm b_norm = b_work.second;
      Tnorm cand_norm = cand.second;
      bool do_something = true;
      if (cand_norm < a_norm && cand_norm < b_norm) {
        if (a_norm < b_norm) {
          b_work = get_comb_pair(cand);
        } else {
          a_work = get_comb_pair(cand);
        }
      } else {
        if (cand_norm < b_norm) {
          b_work = get_comb_pair(cand);
        } else {
          if (cand_norm < a_norm) {
            a_work = get_comb_pair(cand);
          } else {
            do_something = false;
          }
        }
      }
      if (!do_something) {
        std::vector<TcombPair> npair{a_work, b_work};
        return {n_change, npair};
      } else {
        n_change += 1;
      }
    }
  };
  auto erase_entry=[&](TcombPair const& val) -> void {
    auto iter = map.find(val);
    if (iter == map.end()) {
      std::cerr << "SIMP: val shoulf be present in order to be removed\n";
      throw TerminalException{1};
    }
    map.erase(iter);
  };
  auto get_pos=[&](TcombPair const& val) -> size_t {
    auto iter = map.find(val);
    if (iter == map.end()) {
      std::cerr << "SIMP: val shoulf be present in order to get the position\n";
      throw TerminalException{1};
    }
    size_t distance = std::distance(map.begin(), iter);
    return distance;
  };
  Tnorm total_complexity(0);
  for (auto & kv: map) {
    total_complexity += kv.first.second;
  }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
  os << "SIMP: total_complexity=" << total_complexity << "\n";
#endif
  auto look_for_simplification=[&]() -> void {
    // Iterating over the elements and looking for simplifications.
    //
    // The iterators are unstable, so everytime we change the state,
    // the iterators are invalidated and need to be recomputed.
    //
    // As the iteration is being run, the number of elements may decrease.
    // All of this forces a single loop and to handle all the special
    // cases separately.
    size_t u = 0;
    size_t v = 1;
    auto increment_uv=[&]() -> bool {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      os << "SIMP: Starting increment_uv with u=" << u << " v=" << v << "\n";
#endif
      if (v < map.size() - 1) {
        v += 1;
      } else {
        u += 1;
        v = u + 1;
        if (u == map.size() - 1) {
          return false;
        }
      }
      return true;
    };
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
    size_t pos = 0;
    os << "SIMP: starting with the following matrices\n";
    for (auto & kv: map) {
      os << "SIMP: pos=" << pos << " nonce=" << kv.second << " norm=" << kv.first.second << " eM=\n";
      WriteMatrix(os, kv.first.first);
      pos += 1;
    }
#endif
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    size_t n_operation = 0;
#endif
    while(true) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      os << "SIMP: Addressing u=" << u << " v=" << v << " n_operation=" << n_operation << " |set|=" << map.size() << "\n";
#endif
      auto iter1 = map.begin();
      std::advance(iter1, u);
      TcombPair a = iter1->first;
      size_t nonce_a = iter1->second;
      //
      auto iter2 = map.begin();
      std::advance(iter2, v);
      TcombPair b = iter2->first;
      size_t nonce_b = iter2->second;
      std::pair<size_t, size_t> nonce_pair{nonce_a, nonce_b};
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      os << "SIMP:   Complexities, a.second=" << a.second << " b.second=" << b.second << "\n";
#endif
      //
      bool already_treated = false;
      std::pair<size_t, std::vector<TcombPair>> pair{0,{}};
      if (set_treated.find(nonce_pair) != set_treated.end()) {
        already_treated = true;
      } else {
        pair = f_reduce(a, b);
      }
      size_t n_changes = pair.first;
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
      os << "SIMP:   n_changes=" << n_changes << " already_treated=" << already_treated << " |set_treated|=" << set_treated.size() << "\n";
#endif
      if (n_changes > 0) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
        n_operation += 1;
#endif
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
        os << "SIMP:   Working with a=\n";
        WriteMatrix(os, a.first);
        os << "SIMP:   and b=\n";
        WriteMatrix(os, b.first);
#endif
        std::vector<size_t> att;
        std::vector<TcombPair> new_elt;
        size_t min_distance = miss_val;
        bool a_attained = false;
        bool b_attained = false;
        for (auto & ent: pair.second) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
          os << "SIMP:  ent.comp=" << ent.second << " elt.eM=\n";
          WriteMatrix(os, ent.first);
#endif
          auto iter = map.find(ent);
          if (iter == map.end()) {
            new_elt.push_back(ent);
          } else {
            size_t distance = std::distance(map.begin(), iter);
            if (distance < min_distance) {
              min_distance = distance;
            }
            if (distance == u) {
              a_attained = true;
            }
            if (distance == v) {
              b_attained = true;
            }
            if (distance != u && distance != v) {
              att.push_back(distance);
            }
          }
        }
        if (!a_attained) {
          erase_entry(a);
        }
        if (!b_attained) {
          erase_entry(b);
        }
        if (a_attained && b_attained) {
          std::cerr << "SIMP: if n_changes > 0 then either a or b is removed\n";
          throw TerminalException{1};
        }
        size_t min_distance_bis = miss_val;
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
        os << "SIMP:  |new_elt|=" << new_elt.size() << "\n";
#endif
        for (auto & elt: new_elt) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
          os << "SIMP:  elt.comp=" << elt.second << " elt.eM=\n";
          WriteMatrix(os, elt.first);
#endif
          map[elt] = nonce;
          nonce += 1;
          size_t distance = get_pos(elt);
          if (distance < min_distance_bis) {
            min_distance_bis = distance;
          }
        }
        if (min_distance_bis == miss_val) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
          os << "SIMP:  Scenario A\n";
#endif
          // This scenario occurs if the new found generators are already present
          // Two scenarios
          if (!a_attained) {
            // a was removed,
            if (u == map.size()) { // This can happen if a and b are removed and there is nothing after.
              break;
            }
            if (u == map.size() - 1) {
              // This can happen if a is removed but either b or something is after and that is it.
              // So, no more operation possible
              break;
            }
            // Left u to the same value as if u remains the same, we are in the next element.
            v = u + 1; // This is a valid position, continuing
          } else {
            // a still exists
            if (v == map.size()) {
              // v is not invalid, going to the next if possible.
              if (u < map.size() - 2) { // Enough space to make something work
                u += 1;
                v = u + 1;
              } else {
                break;
              }
            } else {
              // No change in value of v, but since b is dropped, we go to the next one.
            }
          }
        } else {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
          os << "SIMP:  Scenario B, |new_elt|=" << new_elt.size() << " min_distance_bis=" << min_distance_bis << "\n";
#endif
          // We have a new generator, adjusting accordingly.
          if (!a_attained) {
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION_EXTENSIVE
            os << "SIMP:  Scenario B, A\n";
#endif
            if (min_distance_bis == map.size() - 1) {
              // We reach the end of what we can do.
              break;
            }
            u = min_distance_bis;
            v = min_distance_bis + 1;
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
            os << "SIMP:   setup A\n";
#endif
          } else {
            size_t pos_a = get_pos(a);
            if (pos_a < min_distance_bis) {
              v = min_distance_bis; // Setting up to where we are and incrementing.
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
              os << "SIMP:   setup B(1)\n";
#endif
              bool test = increment_uv();
              if (!test) {
                break;
              }
            } else {
              u = min_distance_bis;
              v = min_distance_bis + 1;
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
              os << "SIMP:   setup B(2)\n";
#endif
            }
          }
        }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
        os << "SIMP:   a/b_attained=" << a_attained << "/" << b_attained << " min_distance_bis=" << min_distance_bis << "\n";
#endif
      } else {
        if (!already_treated) {
          set_treated.insert(nonce_pair);
        }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
        os << "SIMP:   no change, incrementing u / v\n";
#endif
        bool test = increment_uv();
        if (!test) {
          break;
        }
      }
    }
  };


  while(true) {
    look_for_simplification();
    Tnorm new_complexity(0);
    for (auto & kv: map) {
      new_complexity += kv.first.second;
    }
#ifdef DEBUG_MATRIX_GROUP_SIMPLIFICATION
    os << "SIMP: total_complexity=" << total_complexity << " new_complexity=" << new_complexity << "\n";
#endif
    if (total_complexity == new_complexity) {
      break;
    }
    total_complexity = new_complexity;
  }
  std::vector<Ttype> new_list_gens;
  for (auto & kv: map) {
    new_list_gens.push_back(kv.first.first.first);
  }
  return new_list_gens;
}

template<typename Tnorm, typename Ttype, typename Fcomplexity, typename Finvers, typename Fproduct>
std::vector<Ttype> ExhaustiveReductionComplexity(std::vector<Ttype> const& ListM, Fcomplexity f_complexity, Finvers f_invers, Fproduct f_product, [[maybe_unused]] std::ostream& os) {
  std::unordered_set<Ttype> SetMred;
  for (auto & eM : ListM) {
    Ttype eM_inv = Inverse(eM);
    if (eM < eM_inv) {
      SetMred.insert(eM);
    } else {
      SetMred.insert(eM_inv);
    }
  }
  std::vector<Ttype> ListMred;
  for (auto & eM: SetMred) {
    ListMred.push_back(eM);
  }
  if (ListMred.size() <= 1) {
    return ListMred;
  }
  return ExhaustiveReductionComplexityKernel<Tnorm,Ttype,Fcomplexity,Finvers,Fproduct>(ListMred, f_complexity, f_invers, f_product, os);
}

template<typename T>
std::vector<MyMatrix<T>> ExhaustiveReductionComplexityGroupMatrix(std::vector<MyMatrix<T>> const& ListM, std::ostream& os) {
#ifdef TRACK_INFO_MATRIX_GROUP_SIMPLIFICATION
  write_matrix_group(ListM, "Call_to_ExhaustiveReductionComplexityGroupMatrix");
#endif
  auto f_complexity=[&](MyMatrix<T> const& M) -> T {
    return get_ell1_complexity_measure(M);
  };
  auto f_invers=[](MyMatrix<T> const& M) -> MyMatrix<T> {
    return Inverse(M);
  };
  auto f_product=[](MyMatrix<T> const& A, MyMatrix<T> const& B) -> MyMatrix<T> {
    return A * B;
  };
  return ExhaustiveReductionComplexity<T,MyMatrix<T>,decltype(f_complexity),decltype(f_invers),decltype(f_product)>(ListM, f_complexity, f_invers, f_product, os);
}

template<typename T, typename Telt>
std::pair<std::vector<MyMatrix<T>>, std::vector<Telt>> ExhaustiveReductionComplexityGroupMatrixPerm(std::vector<MyMatrix<T>> const& ListM, std::vector<Telt> const& ListPerm, [[maybe_unused]] std::ostream& os) {
  using Ttype = std::pair<MyMatrix<T>, Telt>;
  auto f_complexity=[&](Ttype const& pair) -> T {
    return get_ell1_complexity_measure(pair.first);
  };
  auto f_invers=[](Ttype const& pair) -> Ttype {
    return {Inverse(pair.first), Inverse(pair.second)};
  };
  auto f_product=[](Ttype const& p1, Ttype const& p2) -> Ttype {
    return {p1.first * p2.first, p1.second * p2.second};
  };
  std::vector<Ttype> ListPair;
  size_t n_gen = ListM.size();
  for (size_t i_gen=0; i_gen<n_gen; i_gen++) {
    Ttype pair{ListM[i_gen], ListPerm[i_gen]};
    ListPair.push_back(pair);
  }
  std::vector<Ttype> RetPair = ExhaustiveReductionComplexity<T,Ttype,decltype(f_complexity),decltype(f_invers),decltype(f_product)>(ListPair, f_complexity, f_invers, f_product, os);
  std::vector<MyMatrix<T>> RetListM;
  std::vector<Telt> RetListPerm;
  size_t n_gen_ret = RetPair.size();
  for (size_t i_gr=0; i_gr<n_gen_ret; i_gr++) {
    RetListM.push_back(RetPair[i_gr].first);
    RetListPerm.push_back(RetPair[i_gr].second);
  }
  return {RetListM, RetListPerm};
}



template<typename T>
std::vector<MyMatrix<T>> Exhaust_get_total_generators(std::vector<MyMatrix<T>> const& list_mat) {
  std::vector<MyMatrix<T>> list_tot;
  for (auto & eMat: list_mat) {
    MyMatrix<T> eMat_inv = Inverse(eMat);
    list_tot.push_back(eMat);
    list_tot.push_back(eMat_inv);
  }
  return list_tot;
}

/*
  Apply a number of exhaustive tricks in order to reduce the size of the vector
 */
template<typename T>
MyVector<T> ExhaustiveVectorSimplificationKernel(MyVector<T> const& V, std::vector<MyMatrix<T>> const& list_mat) {
  int n = V.size();
  auto f_norm=[&](MyVector<T> const& v) -> T {
    T norm(0);
    for (int i=0; i<n; i++) {
      norm += T_abs(v(i));
    }
    return norm;
  };
  MyVector<T> V_work = V;
  T norm_work = f_norm(V);
  while(true) {
    int n_succ = 0;
    for (auto & mat : list_mat) {
      MyVector<T> V_cand = mat.transpose() * V_work;
      T norm_cand = f_norm(V_cand);
      if (norm_cand < norm_work) {
        V_work = V_cand;
        norm_work = norm_cand;
        n_succ += 1;
      }
    }
    if (n_succ == 0) {
      return V_work;
    }
  }
}

template<typename T>
MyVector<T> ExhaustiveVectorSimplification(MyVector<T> const& V, std::vector<MyMatrix<T>> const& list_mat) {
  std::vector<MyMatrix<T>> list_mat_tot = Exhaust_get_total_generators(list_mat);
  return ExhaustiveVectorSimplificationKernel(V, list_mat_tot);
}


template<typename T>
std::vector<MyVector<T>> ExhaustiveVectorSimplifications(std::vector<MyVector<T>> const& list_V, std::vector<MyMatrix<T>> const& list_mat) {
  std::vector<MyMatrix<T>> list_mat_tot = Exhaust_get_total_generators(list_mat);
  std::vector<MyVector<T>> list_V_red;
  for (auto& eV: list_V) {
    MyVector<T> eV_red = ExhaustiveVectorSimplificationKernel(eV, list_mat_tot);
    list_V_red.push_back(eV_red);
  }
  return list_V_red;
}

/*
  Take a matrix M in argument.
  Look for ways to multiply by elements of list_mat
  so that the L1 norm of M G decreases
 */
template<typename T>
MyMatrix<T> ExhaustiveMatrixCosetSimplificationKernel(MyMatrix<T> const& M, std::vector<MyMatrix<T>> const& list_mat) {
  int n = M.rows();
  auto f_norm=[&](MyMatrix<T> const& Hin) -> T {
    T norm(0);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        norm += T_abs(Hin(i,j));
      }
    }
    return norm;
  };
  MyMatrix<T> M_work = M;
  T norm_work = f_norm(M);
  while(true) {
    int n_succ = 0;
    for (auto & mat : list_mat) {
      MyMatrix<T> M_cand = M_work * mat;
      T norm_cand = f_norm(M_cand);
      if (norm_cand < norm_work) {
        M_work = M_cand;
        norm_work = norm_cand;
        n_succ += 1;
      }
    }
    if (n_succ == 0) {
      return M_work;
    }
  }
}


template<typename T>
std::vector<MyMatrix<T>> ExhaustiveMatrixCosetSimplifications(std::vector<MyMatrix<T>> const& list_cos, std::vector<MyMatrix<T>> const& list_mat) {
  std::vector<MyMatrix<T>> list_mat_tot = Exhaust_get_total_generators(list_mat);
  std::vector<MyMatrix<T>> list_cos_red;
  for (auto& eCos: list_cos) {
    MyMatrix<T> eCos_red = ExhaustiveMatrixCosetSimplificationKernel(eCos, list_mat_tot);
    list_cos_red.push_back(eCos_red);
  }
  return list_cos_red;
}

/*
  The result of the simplification algorithm for double cosets
 */
template<typename T>
struct DoubleCosetSimplification {
  MyMatrix<T> u_red; // The u reduction element
  MyMatrix<T> d_cos_red; // The reduced coset
  MyMatrix<T> v_red; // The v_reduction element
};


// The double coset is U x V
template<typename T>
DoubleCosetSimplification<T> ExhaustiveMatrixDoubleCosetSimplifications(MyMatrix<T> const& d_cos, std::vector<MyMatrix<T>> const& u_gens, std::vector<MyMatrix<T>> const& v_gens, size_t const& max_iter, [[maybe_unused]] std::ostream& os) {
  std::vector<MyMatrix<T>> u_gens_tot = Exhaust_get_total_generators(u_gens);
  std::vector<MyMatrix<T>> v_gens_tot = Exhaust_get_total_generators(v_gens);
  int n_gens_u = u_gens_tot.size();
  int n_gens_v = v_gens_tot.size();
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
  os << "DCOS_SIMP: ExhaustiveMatrixDoubleCosetSimplifications |u_gens_tot|=" << u_gens_tot.size() << " |v_gens_tot|=" << v_gens_tot.size() << "\n";
#endif
  int n = d_cos.rows();
  auto f_norm=[&](MyMatrix<T> const& Hin) -> T {
    T norm(0);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        norm += T_abs(Hin(i,j));
      }
    }
    return norm;
  };
  MyMatrix<T> u_red = IdentityMat<T>(n);
  MyMatrix<T> v_red = IdentityMat<T>(n);
  MyMatrix<T> d_cos_work = d_cos;
  T norm_work = f_norm(d_cos);
  std::vector<int> indices_u;
  for (int i=0; i<n_gens_u; i++) {
    indices_u.push_back(i);
  }
  std::vector<int> indices_v;
  for (int i=0; i<n_gens_v; i++) {
    indices_v.push_back(i);
  }


  auto f_search_uv=[&]() -> bool {
    for (int u=0; u<n_gens_u; u++) {
      for (int v=0; v<n_gens_v; v++) {
        int u2 = indices_u[u];
        int v2 = indices_v[v];
        MyMatrix<T> const& u_gen = u_gens_tot[u2];
        MyMatrix<T> const& v_gen = v_gens_tot[v2];
        MyMatrix<T> d_cos_cand = u_gen * d_cos_work * v_gen;
        T norm_cand = f_norm(d_cos_cand);
        if (norm_cand < norm_work) {
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
          os << "DCOS_SIMP: Improving with u2=" << u2 << " v2=" << v2 << " XXX u=" << u << " v=" << v << "\n";
#endif
          d_cos_work = d_cos_cand;
          norm_work = norm_cand;
          u_red = u_gen * u_red;
          v_red = v_red * v_gen;
          return true;
        }
      }
    }
    return false;
  };
  auto f_search_u=[&]() -> bool {
    for (int u=0; u<n_gens_u; u++) {
      int u2 = indices_u[u];
      MyMatrix<T> const& u_gen = u_gens_tot[u2];
      MyMatrix<T> d_cos_cand = u_gen * d_cos_work;
      T norm_cand = f_norm(d_cos_cand);
      if (norm_cand < norm_work) {
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
        os << "DCOS_SIMP: Improving with u2=" << u2 << "\n";
#endif
        d_cos_work = d_cos_cand;
        norm_work = norm_cand;
        u_red = u_gen * u_red;
        return true;
      }
    }
    return false;
  };
  auto f_search_v=[&]() -> bool {
    for (int v=0; v<n_gens_v; v++) {
      int v2 = indices_v[v];
      MyMatrix<T> const& v_gen = v_gens_tot[v2];
      MyMatrix<T> d_cos_cand = d_cos_work * v_gen;
      T norm_cand = f_norm(d_cos_cand);
      if (norm_cand < norm_work) {
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
        os << "DCOS_SIMP: Improving with v2=" << v2 << "\n";
#endif
        d_cos_work = d_cos_cand;
        norm_work = norm_cand;
        v_red = v_red * v_gen;
        return true;
      }
    }
    return false;
  };
  auto f_search=[&]() -> bool {
    int chosen_method = 3;
    f_random_transpose(indices_u);
    f_random_transpose(indices_v);
    if (chosen_method == 1) {
      return f_search_uv();
    }
    if (chosen_method == 2) {
      bool test_u = f_search_u();
      if (test_u) {
        return true;
      }
      return f_search_v();
    }
    if (chosen_method == 3) {
      bool test_u = f_search_u();
      if (test_u) {
        return true;
      }
      bool test_v = f_search_v();
      if (test_v) {
        return true;
      }
      return f_search_uv();
    }
    return false;
  };
  size_t n_iter = 0;
  while(true) {
    bool test = f_search();
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
    os << "DCOS_SIMP: ExhaustiveMatrixDoubleCosetSimplifications n_iter=" << n_iter << " norm_work=" << norm_work << "\n";
#endif
    if (!test) {
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
      os << "DCOS_SIMP: ExhaustiveMatrixDoubleCosetSimplifications n_final_iter(A)=" << n_iter << " norm_work=" << norm_work << "\n";
#endif
      return {u_red, d_cos_work, v_red};
    }
    if (n_iter == max_iter) {
#ifdef DEBUG_MATRIX_DOUBLE_COSET_SIMPLIFICATION
      os << "DCOS_SIMP: ExhaustiveMatrixDoubleCosetSimplifications n_final_iter(B)=" << n_iter << " norm_work=" << norm_work << "\n";
#endif
      return {u_red, d_cos_work, v_red};
    }
    n_iter += 1;
  }
}

// clang-format off
#endif  // SRC_GROUP_MATRIXGROUP_H_
