// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Group.h"
#include "Permutation.h"
#include "PolyNorm_PackingContacts.h"
#include <fstream>
#include <iomanip>
#include <sstream>
// clang-format on

// The classes of lattice cross-polytopes as written by
// enumerate_lattice_crosspolytopes (Paper_poly_norm/LatticeCrossPolytopes):
// "index i free|non-free denominator E elements: (a,b,...) ...".
std::vector<PCClass> ReadClasses(std::string const &file, int n) {
  std::ifstream in(file);
  std::vector<PCClass> L;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::string w;
    int64_t idx, E;
    ss >> w >> idx >> w >> w >> E >> w;
    PCClass C;
    C.E = E;
    std::string tok;
    while (ss >> tok) {
      std::vector<int64_t> v;
      std::string t = tok.substr(1, tok.size() - 2);
      std::stringstream ts(t);
      std::string x;
      while (std::getline(ts, x, ','))
        v.push_back(std::stoll(x));
      if (static_cast<int>(v.size()) != n) {
        std::cerr << "Wrong dimension in the classes file\n";
        throw TerminalException{1};
      }
      C.elts.push_back(v);
    }
    L.push_back(C);
  }
  return L;
}

// The contact sets of a class (stage 1 output): one set per line, the
// vectors (denominator E) concatenated.
std::vector<std::vector<std::vector<int64_t>>>
ReadSets(std::string const &file, int n, int64_t E) {
  std::vector<std::vector<std::vector<int64_t>>> L;
  std::ifstream in(file);
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ss(line);
    std::vector<int64_t> flat;
    int64_t x;
    while (ss >> x)
      flat.push_back(x);
    std::vector<std::vector<int64_t>> S;
    for (size_t i = 0; i + n <= flat.size(); i += n)
      S.push_back(std::vector<int64_t>(flat.begin() + i, flat.begin() + i + n));
    // The unit vectors first (the shortcut of SearchSet needs it).
    std::stable_partition(S.begin(), S.end(), [&](auto const &v) {
      int nz = 0;
      for (auto &c : v)
        if (c != 0)
          nz++;
      return nz == 1 && std::llabs(v[0] + v[1] + v[2] + v[3]) == E;
    });
    L.push_back(S);
  }
  return L;
}

int main(int argc, char *argv[]) {
  maybe_install_tls_gmp_pool();
  HumanTime time;
  try {
    if (argc != 8) {
      std::cerr << "POLYNORM_PackingContacts [FileEXT] [FileClasses] "
                   "[DirSets] [classes] [density] [n_thread] [OutFile]\n";
      std::cerr << "  classes : a comma separated list of class indices\n";
      std::cerr << "  density : only the lattices of density above it are "
                   "kept (the best known value)\n";
      return -1;
    }
    using T = mpq_class;
    using Tidx = uint32_t;
    using Telt = permutalib::SingleSidedPerm<Tidx>;
    using Tgroup = permutalib::Group<Telt, mpz_class>;
    std::string FileEXT = argv[1], FileClasses = argv[2], DirSets = argv[3];
    std::vector<int> ListIc;
    {
      std::stringstream ts(argv[4]);
      std::string x;
      while (std::getline(ts, x, ','))
        ListIc.push_back(std::stoi(x));
    }
    double density = std::stod(argv[5]);
    int n_thread = ParseScalar<int>(argv[6]);
    std::string OutFile = argv[7];
    MyMatrix<T> EXT = ReadMatrixFile<T>(FileEXT);
    PCBody<T> body = PC_BuildBody<T, Tgroup>(EXT, std::cerr);
    int n = body.n;
    std::vector<PCClass> ListClass = ReadClasses(FileClasses, n);
    double vol = UniversalScalarConversion<double, T>(body.vol_P);
    double value_bound = vol / density * (1 + 1e-9);
    auto prt = [&](std::ostream &os) -> void {
      os << std::setprecision(12);
      os << "n=" << n << " facets of D: " << body.m << " (" << body.n_orbit
         << " orbits), vol(P)=" << vol << ", density threshold " << density
         << "\n";
      for (int ic : ListIc) {
        HumanTime tc;
        PCSearch<T> search(body, ListClass, std::cerr);
        auto Sets = ReadSets(DirSets + "/class_" + std::to_string(ic) + ".sets",
                             n, ListClass[ic].E);
        std::vector<std::pair<size_t, int>> tasks;
        for (size_t is = 0; is < Sets.size(); is++)
          for (int f = 0; f < body.m; f++)
            if (body.is_rep[f])
              tasks.push_back({is, f});
        std::atomic<size_t> next(0);
        std::atomic<bool> done(false);
        std::thread monitor([&]() {
          int sec = 0;
          while (!done) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            if (++sec % 20 == 0)
              std::cerr << "PROGRESS class " << ic << " t=" << sec
                        << "s tasks started " << next << "/" << tasks.size()
                        << " nodes " << search.n_node << " lp " << search.n_lp
                        << " triple_lp " << search.n_triple_lp
                        << " assignments " << search.n_leaf_raw
                        << " spaces " << search.ListLeaf.size() << std::endl;
          }
        });
        std::vector<std::thread> th;
        for (int t = 0; t < n_thread; t++)
          th.emplace_back([&]() {
            while (true) {
              size_t i = next.fetch_add(1);
              if (i >= tasks.size())
                break;
              search.SearchSet(ic, Sets[tasks[i].first], tasks[i].second);
            }
          });
        for (auto &x : th)
          x.join();
        done = true;
        monitor.join();
        size_t n_leaf = search.ListLeaf.size();
        os << "class " << ic << " (index " << ListClass[ic].elts.size()
           << "): " << Sets.size() << " contact sets, " << search.n_node
           << " nodes, " << search.n_lp << " linear programs, "
           << search.n_leaf_raw << " assignments, " << n_leaf
           << " distinct affine spaces, stage 2 time " << tc << std::endl;
        // Stage 3.
        std::vector<PCCandidate<T>> cand;
        std::mutex mtx;
        std::atomic<size_t> next3(0), stalled(0);
        std::vector<std::thread> th3;
        for (int t = 0; t < n_thread; t++)
          th3.emplace_back([&]() {
            size_t st = 0;
            while (true) {
              size_t i = next3.fetch_add(1);
              if (i >= n_leaf)
                break;
              auto R = PC_Stage3(body, ListClass[ic], search.ListLeafEq[i],
                                 static_cast<int>(i), value_bound, st);
              std::lock_guard<std::mutex> lock(mtx);
              for (auto &c : R)
                cand.push_back(c);
            }
            stalled += st;
          });
        for (auto &x : th3)
          x.join();
        size_t n_adm = 0, n_unknown = 0;
        double best = std::numeric_limits<double>::max();
        for (auto &c : cand) {
          if (c.min_gauge < 0)
            n_unknown++;
          else if (c.min_gauge >= 1 - 1e-7) {
            n_adm++;
            best = std::min(best, c.value);
          }
        }
        os << "  stage 3: critical points below the threshold value "
           << cand.size() << ", admissible " << n_adm
           << ", admissibility not computed " << n_unknown
           << ", homotopy paths stopped in the middle " << stalled;
        if (n_adm > 0)
          os << ", best density " << vol / best;
        os << ", total time " << tc << std::endl;
      }
    };
    FILE_PrintStderrStdoutFile(OutFile, prt);
    std::cerr << "Normal termination of POLYNORM_PackingContacts\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in POLYNORM_PackingContacts\n";
    exit(e.eVal);
  }
  runtime(time);
}
