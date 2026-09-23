// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "JointCoveringDouble.h"
#include <iostream>
#include <csignal>
#include <sys/wait.h>
#include <unistd.h>
// clang-format on

/*
  Search tool for periodic sphere coverings in double precision: covering
  density evaluation, joint local descent on (Q, C), and multistart, for
  any dimension and any number of cosets. See JointCoveringDouble.h for
  the method; candidates found here are to be re-verified by the exact
  machinery.

  Configuration files: a header "n m", then the n rows of Q, then the m
  rows of the cosets (the first being zero).
 */

using namespace joint_covering_double;

int main(int argc, char *argv[]) {
  try {
    if (argc < 2) {
      std::cerr << "PERIODIC_JointCoveringSearch [mode] ...\n";
      std::cerr << "\n";
      std::cerr << "modes:\n";
      std::cerr << "  evaluate [config]\n";
      std::cerr << "      covering density, covering radius and cell-class "
                   "count of the configuration\n";
      std::cerr << "  descend [config] [out_config] [rounds]\n";
      std::cerr << "      joint local descent from the configuration; the "
                   "best found is written to out_config\n";
      std::cerr << "  multistart [n] [m] [count] [out_config] [seed] "
                << "[rounds]\n";
      std::cerr << "  multistart-alt [n] [m] [count] [out_config] [seed]\n";
      std::cerr << "  certify [config]\n";
      std::cerr << "      LP rigidity certificate at the configuration: is it "
                   "jammed (no first-order descent direction)?\n";
      std::cerr << "  descend-lp [config] [out_config] [rounds]\n";
      std::cerr << "      joint (Q,C) descent driven by the LP direction step; "
                   "reports whether the limit point is rigid (jammed)\n";
      std::cerr << "  polish [config] [out_config] [rounds]\n";
      std::cerr << "      L-BFGS then LP finish: the recommended optimizer; "
                   "reports rigidity of the limit\n";
      std::cerr << "  multistart-lp [n] [m] [count] [out_config] [seed] "
                   "[max_seconds]\n";
      std::cerr << "  basinhop [config] [out_config] [max_seconds] [seed] "
                   "[collect_below]\n";
      std::cerr << "  evaluate-pc [config]\n";
      std::cerr << "      packing-covering ratio gamma = 2 sqrt(mu^2 / "
                   "lambda^2) of the configuration\n";
      std::cerr << "  descend-pc [config] [out_config] [rounds]\n";
      std::cerr << "  basinhop-pc [config] [out_config] [max_seconds] "
                   "[seed]\n";
      std::cerr << "      the same searches for the packing-covering ratio; "
                   "records: gamma_4 = 1.3625 (Ho_4), gamma_5 = 1.44946 "
                   "(Ho_5), lattice-optimal, non-lattice open\n";
      std::cerr << "  addcoset [config] [out_config]\n";
      std::cerr << "      grow m -> m+1 by adding a coset at the deepest hole "
                   "(structured seeding for a richer coset set)\n";
      std::cerr << "      basin hopping from the seed configuration: cheap "
                   "deep-hole / jitter kicks, forked relaxation, Metropolis "
                   "with cooling; writes the best (and any record) found\n";
      std::cerr << "      LP-direction descents from random seeds, reporting "
                   "the rigid configurations found\n";
      std::cerr << "      alternating SDP/minimax descent from well-rounded "
                   "seeds (fast)\n";
      std::cerr << "      repeated descents from random configurations, "
                   "keeping the best in out_config\n";
      return -1;
    }
    std::string mode = argv[1];
    if (mode == "evaluate" && argc == 3) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      DensityResult dr = CoveringDensity(conf);
      std::cout << "n=" << conf.n << " m=" << conf.m
                << " classes=" << dr.cells.size() << "\n";
      std::cout << "mu2=" << dr.mu2 << "\n";
      printf("theta=%.13f\n", dr.theta);
      return 0;
    }
    if (mode == "qstepcheck" && argc == 3) {
      // gradient check of the Q-step barrier, then a full Q-step: verify it
      // improves (or holds) the covering density at fixed cosets.
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      DensityResult dr = CoveringDensity(conf);
      QStepData qd = BuildQStepData(dr.cells, conf.C);
      VectorXd x(qd.dim);
      for (int u = 0; u < qd.dim; u++) {
        x(u) = conf.Q(qd.basis[u].first, qd.basis[u].second);
      }
      double val;
      VectorXd g;
      MatrixXd H;
      QStepBarrier(qd, x, 100.0, val, &g, &H);
      double eps = 1e-6, worst = 0;
      for (int u = 0; u < qd.dim; u++) {
        VectorXd xp = x, xm = x;
        xp(u) += eps;
        xm(u) -= eps;
        double vp, vm;
        QStepBarrier(qd, xp, 100.0, vp, nullptr, nullptr);
        QStepBarrier(qd, xm, 100.0, vm, nullptr, nullptr);
        double fd = (vp - vm) / (2 * eps);
        worst = std::max(worst, std::abs(fd - g(u)) /
                                    std::max(1.0, std::abs(fd)));
      }
      printf("worst relative gradient error = %.3e\n", worst);
      MatrixXd Qopt = QStep(dr.cells, conf.C, conf.Q);
      PeriodicConfig c2 = conf;
      c2.Q = Qopt;
      DensityResult d2 = CoveringDensity(c2);
      printf("density before Q-step = %.12f\n", dr.theta);
      printf("density after  Q-step = %.12f\n", d2.theta);
      return 0;
    }
    if (mode == "gradcheck" && argc == 3) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      DensityResult dr = CoveringDensity(conf);
      Packing pk{conf.n, conf.m, conf.n * (conf.n + 1) / 2};
      Eigen::LLT<MatrixXd> llt(conf.Q);
      MatrixXd Lo = llt.matrixL();
      VectorXd x = Pack(pk, Lo, conf.C);
      VectorXd g;
      double beta = 300.0;
      double F0 = ObjectiveGradient(pk, dr.cells, beta, x, g);
      double eps = 1e-6;
      double worst = 0;
      for (int k = 0; k < x.size(); k++) {
        VectorXd xp = x, xm = x, gd;
        xp(k) += eps;
        xm(k) -= eps;
        double fd = (ObjectiveGradient(pk, dr.cells, beta, xp, gd) -
                     ObjectiveGradient(pk, dr.cells, beta, xm, gd)) /
                    (2 * eps);
        double err = std::abs(fd - g(k)) / std::max(1.0, std::abs(fd));
        worst = std::max(worst, err);
      }
      printf("F=%.12f dim=%d worst relative gradient error=%.3e\n", F0,
             int(x.size()), worst);
      return 0;
    }
    if (mode == "joint" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 60;
      DescendResult res = JointDescent(conf, rounds, std::cerr, true);
      if (!res.success) { std::cerr << "no config evaluated\n"; return 1; }
      printf("theta=%.13f\n", res.theta);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "descend-alt2" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 30;
      DescendResult res = DescendAlt2(conf, rounds, std::cerr, true);
      if (!res.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      printf("theta=%.13f\n", res.theta);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "descend-alt" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 20;
      DescendResult res = DescendAlt(conf, rounds, std::cerr, true);
      if (!res.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      printf("theta=%.13f\n", res.theta);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "descend" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 15;
      DescendResult res = Descend(conf, rounds, std::cerr, true);
      if (!res.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      printf("theta=%.13f\n", res.theta);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "multistart-alt" && (argc == 6 || argc == 7)) {
      int n = atoi(argv[2]);
      int m = atoi(argv[3]);
      int count = atoi(argv[4]);
      unsigned seed = argc == 7 ? unsigned(atol(argv[6])) : 1u;
      std::mt19937_64 gen(seed);
      std::normal_distribution<double> gauss(0.0, 1.0);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      double best = 1e30;
      for (int st = 0; st < count; st++) {
        PeriodicConfig conf;
        conf.n = n;
        conf.m = m;
        // well-rounded seed: R = I + small gaussian, Q = R^T R near the
        // identity, so the tessellation stays cheap and the descent starts
        // in the region where covering-optimal (well-rounded) forms live
        Eigen::MatrixXd R = Eigen::MatrixXd::Identity(n, n);
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
            R(i, j) += 0.25 * gauss(gen);
          }
        }
        conf.Q = R.transpose() * R;
        conf.C.resize(m, n);
        conf.C.row(0).setZero();
        for (int t = 1; t < m; t++) {
          for (int j = 0; j < n; j++) {
            conf.C(t, j) = unif(gen);
          }
        }
        DescendResult res = DescendAlt(conf, 25, std::cerr, false);
        if (res.success) {
          printf("start %d: theta=%.13f%s\n", st, res.theta,
                 res.theta < best ? "  *" : "");
          fflush(stdout);
          if (res.theta < best) {
            best = res.theta;
            WriteConfigFile(argv[5], res.conf);
          }
        } else {
          printf("start %d: failed\n", st);
          fflush(stdout);
        }
      }
      printf("best=%.13f\n", best);
      return 0;
    }
    if (mode == "multistart" && (argc >= 6 && argc <= 8)) {
      int n = atoi(argv[2]);
      int m = atoi(argv[3]);
      int count = atoi(argv[4]);
      unsigned seed = argc >= 7 ? unsigned(atol(argv[6])) : 1u;
      int rounds = argc == 8 ? atoi(argv[7]) : 15;
      std::mt19937_64 gen(seed);
      std::normal_distribution<double> gauss(0.0, 1.0);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      double best = 1e30;
      for (int st = 0; st < count; st++) {
        PeriodicConfig conf;
        conf.n = n;
        conf.m = m;
        Eigen::MatrixXd A(n, n);
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
            A(i, j) = gauss(gen);
          }
        }
        conf.Q = A * A.transpose() +
                 0.05 * Eigen::MatrixXd::Identity(n, n);
        conf.C.resize(m, n);
        conf.C.row(0).setZero();
        for (int t = 1; t < m; t++) {
          for (int j = 0; j < n; j++) {
            conf.C(t, j) = unif(gen);
          }
        }
        DescendResult res = Descend(conf, rounds, std::cerr, false);
        if (res.success) {
          printf("start %d: theta=%.13f%s\n", st, res.theta,
                 res.theta < best ? "  *" : "");
          fflush(stdout);
          if (res.theta < best) {
            best = res.theta;
            WriteConfigFile(argv[5], res.conf);
          }
        } else {
          printf("start %d: failed\n", st);
          fflush(stdout);
        }
      }
      printf("best=%.13f\n", best);
      return 0;
    }
    if (mode == "certify" && argc == 3) {
      // The pure LP rigidity certificate at a fixed configuration: no descent.
      // Reports the scale-free steepest-descent rate of log Theta over the
      // joint (Q, C) directions. rate ~ 0 means no direction lowers the
      // density to first order -- the configuration is rigid (jammed).
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      DensityResult dr = CoveringDensity(conf);
      auto basis = SymBasis(conf.n);
      int dQ = basis.size();
      int dim = dQ + conf.n * (conf.m - 1);
      std::vector<double> R2;
      std::vector<Eigen::VectorXd> g;
      JointGrad(dr.cells, conf.Q, conf.C, basis, R2, g);
      double mu2 = *std::max_element(R2.begin(), R2.end());
      Eigen::VectorXd gdet = Eigen::VectorXd::Zero(dim);
      Eigen::MatrixXd Qinv = conf.Q.inverse();
      for (int u = 0; u < dQ; u++) {
        int k = basis[u].first, l = basis[u].second;
        gdet(u) = (k == l) ? Qinv(k, k) : 2.0 * Qinv(k, l);
      }
      // A small box so the linearization stays valid (a large box lets the
      // linear model "see" descent that would in fact activate other radii):
      // the reported rate is then the genuine one-sided directional derivative.
      double rho_cert = 1e-3;
      DirLPResult lp = DirectionLP(R2, g, gdet, mu2, conf.n, rho_cert);
      double rate = -lp.predicted / rho_cert;
      printf("theta=%.13f\n", dr.theta);
      printf("orbits=%d active=%d rate=%.6e rigid=%d\n",
             (int)dr.cells.size(), lp.n_active, rate,
             rate < 1e-7 ? 1 : 0);
      return 0;
    }
    if (mode == "descend-lp" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 200;
      DescendResult res = JointDescentLP(conf, rounds, std::cerr, true);
      if (!res.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      printf("theta=%.13f\n", res.theta);
      printf("rigid=%d stationarity=%.3e active=%d\n", res.rigid ? 1 : 0,
             res.stationarity, res.n_active);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "polish" && (argc == 4 || argc == 5)) {
      // The recommended optimizer: L-BFGS descends fast to the neighbourhood
      // of a local minimum, then the LP direction step finishes the descent
      // along the non-smooth ridge and certifies whether the limit is rigid.
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 20;
      DescendResult a = Descend(conf, rounds, std::cerr, true);
      if (!a.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      std::cerr << "L-BFGS theta=" << a.theta << "; LP finishing...\n";
      DescendResult b = JointDescentLP(a.conf, 300, std::cerr, true);
      DescendResult r = (b.success && b.theta <= a.theta) ? b : a;
      printf("theta=%.13f\n", r.theta);
      printf("rigid=%d stationarity=%.3e active=%d\n", b.rigid ? 1 : 0,
             b.stationarity, b.n_active);
      WriteConfigFile(argv[3], r.conf);
      return 0;
    }
    if (mode == "multistart-lp" && (argc >= 6 && argc <= 8)) {
      int n = atoi(argv[2]);
      int m = atoi(argv[3]);
      int count = atoi(argv[4]);
      unsigned seed = argc >= 7 ? unsigned(atol(argv[6])) : 1u;
      // Optional wall-clock budget (seconds): when given, run starts until the
      // budget is exhausted rather than a fixed count.
      double max_seconds = argc == 8 ? atof(argv[7]) : 0.0;
      // Hard per-start wall-clock cap. Each start runs in a forked child so a
      // start whose tessellation sends qhull into a multi-hour convex-hull
      // grind (a single uninterruptible qh_new_qhull call, seen on degenerate
      // near-cospherical clouds) can be SIGKILLed and the CPU reclaimed.
      double start_budget = 300.0;
      // The lattice record to beat in this dimension (best lattice covering).
      double record = 0.0;
      if (n == 3) record = 1.4635030689668180;
      if (n == 4) record = 1.7655285081493524;
      if (n == 5) record = 2.1242859089916246;
      std::mt19937_64 gen(seed);
      std::normal_distribution<double> gauss(0.0, 1.0);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      double best = 1e30;
      int n_rigid = 0;
      int n_record = 0;
      auto t_start = std::chrono::steady_clock::now();
      for (int st = 0;; st++) {
        double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(
                             std::chrono::steady_clock::now() - t_start).count();
        if (max_seconds > 0.0) {
          if (elapsed >= max_seconds) break;
        } else {
          if (st >= count) break;
        }
        PeriodicConfig conf;
        conf.n = n;
        conf.m = m;
        Eigen::MatrixXd R = Eigen::MatrixXd::Identity(n, n);
        for (int i = 0; i < n; i++)
          for (int j = 0; j < n; j++)
            R(i, j) += 0.25 * gauss(gen);
        conf.Q = R.transpose() * R;
        conf.C.resize(m, n);
        conf.C.row(0).setZero();
        for (int t = 1; t < m; t++)
          for (int j = 0; j < n; j++)
            conf.C(t, j) = unif(gen);
        // Run this start in a child process with a hard time cap. The child
        // does the whole descent, writes its result and configuration to temp
        // files, and _exit()s; the parent waits at most start_budget seconds
        // and SIGKILLs a child that overruns. Only one child runs at a time.
        std::string tmpstat = std::string(argv[5]) + ".child.stat";
        std::string tmpconf = std::string(argv[5]) + ".child.conf";
        ::remove(tmpstat.c_str());
        pid_t pid = fork();
        if (pid == 0) {
          DescendResult a = Descend(conf, 15, std::cerr, false);
          DescendResult res = a;
          if (a.success) {
            DescendResult b = JointDescentLP(a.conf, 200, std::cerr, false, 90);
            if (b.success && b.theta <= a.theta) res = b;
          }
          FILE *cf = fopen(tmpstat.c_str(), "w");
          if (cf) {
            if (res.success) {
              fprintf(cf, "OK %.15g %d %d %.6e\n", res.theta,
                      res.rigid ? 1 : 0, res.n_active, res.stationarity);
              fclose(cf);
              WriteConfigFile(tmpconf, res.conf);
            } else {
              fprintf(cf, "FAIL\n");
              fclose(cf);
            }
          }
          _exit(0);
        }
        auto ct0 = std::chrono::steady_clock::now();
        bool finished = false;
        while (true) {
          int status;
          if (waitpid(pid, &status, WNOHANG) == pid) { finished = true; break; }
          double el = std::chrono::duration_cast<std::chrono::duration<double>>(
                          std::chrono::steady_clock::now() - ct0).count();
          if (el > start_budget) break;
          usleep(200000);
        }
        if (!finished) {
          kill(pid, SIGKILL);
          int status;
          waitpid(pid, &status, 0);
          printf("start %d: TIMEOUT killed after %.0fs\n", st, start_budget);
          fflush(stdout);
          continue;
        }
        FILE *rf = fopen(tmpstat.c_str(), "r");
        char tag[16] = {0};
        if (!rf || fscanf(rf, "%15s", tag) != 1 ||
            std::string(tag) != "OK") {
          if (rf) fclose(rf);
          printf("start %d: failed\n", st);
          fflush(stdout);
          continue;
        }
        double th = 0, stt = 0;
        int rg = 0, ac = 0;
        if (fscanf(rf, "%lf %d %d %lf", &th, &rg, &ac, &stt) != 4) {
          fclose(rf);
          printf("start %d: parse error\n", st);
          fflush(stdout);
          continue;
        }
        fclose(rf);
        PeriodicConfig rescfg = ReadConfigFile(tmpconf);
        if (rg) {
          n_rigid++;
          char fn[4096];
          snprintf(fn, sizeof(fn), "%s.rigid.%d", argv[5], n_rigid);
          WriteConfigFile(fn, rescfg);
        }
        bool beats = record > 0.0 && th < record - 1e-9;
        if (beats) {
          n_record++;
          char fn[4096];
          snprintf(fn, sizeof(fn), "%s.record.%d", argv[5], n_record);
          WriteConfigFile(fn, rescfg);
        }
        printf("start %d: theta=%.13f rigid=%d active=%d rate=%.2e%s%s\n", st,
               th, rg, ac, stt, th < best ? "  *" : "",
               beats ? "  <<< BEATS RECORD" : "");
        fflush(stdout);
        if (th < best) {
          best = th;
          WriteConfigFile(argv[5], rescfg);
        }
      }
      printf("best=%.13f rigid_count=%d record_beats=%d (record=%.13f)\n", best,
             n_rigid, n_record, record);
      return 0;
    }
    if (mode == "basinhop" && (argc >= 5 && argc <= 7)) {
      // Basin hopping: a random walk in the space of local minima. Each hop
      // applies a cheap geometric kick (a coset moved to the deepest hole, or
      // a single-coset jitter), re-relaxes, and accepts by Metropolis with a
      // cooling temperature. The relaxation runs in a forked child with a hard
      // time cap, so a kick that lands on a degenerate (near-cospherical)
      // configuration -- which the deep-hole move does by construction -- can
      // never wedge the walk in an uninterruptible qhull call.
      PeriodicConfig seed_conf = ReadConfigFile(argv[2]);
      double max_seconds = atof(argv[4]);
      unsigned seed = argc >= 6 ? unsigned(atol(argv[5])) : 1u;
      // Collect every distinct local optimum below this threshold into a
      // catalog; after the walk the best entries are re-polished by the full
      // L-BFGS descent and certified, so the run returns the remarkable
      // configurations it met, not only the single best.
      double collect_below = argc == 7 ? atof(argv[6]) : 0.0;
      std::vector<std::pair<double, PeriodicConfig>> catalog;
      auto collect = [&](double th, PeriodicConfig const &c) {
        if (collect_below <= 0.0 || th >= collect_below) return;
        for (auto &e : catalog)
          if (std::abs(e.first - th) < 1e-6) return;
        catalog.emplace_back(th, c);
        std::sort(catalog.begin(), catalog.end(),
                  [](auto const &a, auto const &b) { return a.first < b.first; });
        if (catalog.size() > 40) catalog.pop_back();
      };
      int n = seed_conf.n, m = seed_conf.m;
      // cap the search-time point cloud: good configurations of m cosets need
      // ~1500*m points, so 4000*m rejects the ball-growth blow-ups (a bad kick
      // fails fast) while leaving healthy configurations untouched
      SearchMaxPts() = 4000 * m;
      double record = 0.0;
      if (n == 3) record = 1.4635030689668180;
      if (n == 4) record = 1.7655285081493524;
      if (n == 5) record = 2.1242859089916246;
      std::mt19937_64 rng(seed);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      std::normal_distribution<double> gauss(0.0, 1.0);
      double start_budget = 120.0;
      std::string tmpstat = std::string(argv[3]) + ".child.stat";
      std::string tmpconf = std::string(argv[3]) + ".child.conf";
      auto relax_forked = [&](PeriodicConfig const &c, PeriodicConfig &out,
                              double &theta) -> bool {
        ::remove(tmpstat.c_str());
        pid_t pid = fork();
        if (pid == 0) {
          DescendResult r = RelaxCheap(c, 4);
          FILE *f = fopen(tmpstat.c_str(), "w");
          if (f) {
            if (r.success) {
              fprintf(f, "OK %.15g\n", r.theta);
              fclose(f);
              WriteConfigFile(tmpconf, r.conf);
            } else {
              fprintf(f, "FAIL\n");
              fclose(f);
            }
          }
          _exit(0);
        }
        auto ct0 = std::chrono::steady_clock::now();
        bool fin = false;
        while (true) {
          int status;
          if (waitpid(pid, &status, WNOHANG) == pid) { fin = true; break; }
          double el = std::chrono::duration_cast<std::chrono::duration<double>>(
                          std::chrono::steady_clock::now() - ct0).count();
          if (el > start_budget) break;
          usleep(200000);
        }
        if (!fin) {
          kill(pid, SIGKILL);
          int status;
          waitpid(pid, &status, 0);
          return false;
        }
        FILE *f = fopen(tmpstat.c_str(), "r");
        char tag[16] = {0};
        if (!f || fscanf(f, "%15s", tag) != 1 || std::string(tag) != "OK") {
          if (f) fclose(f);
          return false;
        }
        if (fscanf(f, "%lf", &theta) != 1) { fclose(f); return false; }
        fclose(f);
        out = ReadConfigFile(tmpconf);
        return true;
      };
      PeriodicConfig curc;
      double curth;
      bool init_ok = false;
      for (int attempt = 0; attempt < 8 && !init_ok; attempt++) {
        // a random seed may need a larger enumeration ball than a healthy
        // configuration just to be evaluated once: raise the cloud cap
        // progressively during the init attempts (the fork watchdog still
        // bounds the cost), then restore the tight cap for the walk
        SearchMaxPts() = 4000 * m * (1 + attempt);
        PeriodicConfig s2 = seed_conf;
        if (attempt > 0) {
          // the seed relaxed badly (an expensive/degenerate first tessellation);
          // jitter the cosets to escape it before giving up
          for (int t = 1; t < m; t++)
            for (int j = 0; j < n; j++)
              s2.C(t, j) += 0.10 * gauss(rng);
        }
        init_ok = relax_forked(s2, curc, curth);
        if (!init_ok)
          printf("init attempt %d failed, retrying with a jittered seed\n",
                 attempt);
        fflush(stdout);
      }
      if (!init_ok) {
        std::cerr << "basinhop: initial relaxation failed after retries\n";
        return 1;
      }
      SearchMaxPts() = 4000 * m;
      PeriodicConfig bestc = curc;
      double bestth = curth;
      collect(curth, curc);
      WriteConfigFile(argv[3], bestc);
      printf("init: theta=%.13f\n", curth);
      fflush(stdout);
      double T = 0.04 * curth, Tmin = 1e-4;
      auto t0 = std::chrono::steady_clock::now();
      int n_acc = 0, n_hop = 0;
      for (int hop = 0;; hop++) {
        double el = std::chrono::duration_cast<std::chrono::duration<double>>(
                        std::chrono::steady_clock::now() - t0).count();
        if (el > max_seconds) break;
        PeriodicConfig trial = curc;
        std::vector<CellClass> cells;
        try {
          cells = DelaunayCellClasses(trial);
        } catch (std::runtime_error &e) {
          continue;
        }
        int t = (m > 2) ? 1 + int(rng() % (m - 1)) : 1;
        double u = unif(rng);
        if (u < 0.15) {
          // occasional deep-hole placement -- a big, combinatorics-changing
          // jump; a sizeable jitter keeps it off the exactly-cospherical
          // locus (which would make the child tessellation very slow)
          Eigen::VectorXd z = DeepestHole(cells, trial.Q, trial.C);
          for (int j = 0; j < n; j++)
            trial.C(t, j) = (z(j) - std::floor(z(j))) + 0.05 * gauss(rng);
        } else {
          // small local single-coset jitter -- keeps the configuration near a
          // good one so its tessellation stays cheap and hops stay frequent
          double sigma = (u < 0.6) ? 0.06 : 0.12;
          for (int j = 0; j < n; j++) trial.C(t, j) += sigma * gauss(rng);
        }
        PeriodicConfig candc;
        double candth;
        n_hop++;
        if (!relax_forked(trial, candc, candth)) {
          printf("hop %d: relax timeout/fail\n", hop);
          fflush(stdout);
          continue;
        }
        collect(candth, candc);
        double dth = candth - curth;
        bool acc = (dth < 0.0) || (unif(rng) < std::exp(-dth / T));
        if (acc) { curc = candc; curth = candth; n_acc++; }
        bool rec = record > 0.0 && candth < record - 1e-9;
        if (candth < bestth - 1e-12) {
          bestth = candth;
          bestc = candc;
          WriteConfigFile(argv[3], bestc);
          if (rec) {
            char fn[4096];
            snprintf(fn, sizeof(fn), "%s.record", argv[3]);
            WriteConfigFile(fn, bestc);
          }
        }
        printf("hop %d: theta=%.13f %s cur=%.13f best=%.13f T=%.4f%s\n", hop,
               candth, acc ? "acc" : "rej", curth, bestth, T,
               rec ? "  <<< BEATS RECORD" : "");
        fflush(stdout);
        T = std::max(Tmin, T * 0.995);
      }
      printf("best=%.13f accept_rate=%.2f (record=%.13f)\n", bestth,
             n_hop ? double(n_acc) / n_hop : 0.0, record);
      WriteConfigFile(argv[3], bestc);
      // Post-polish: the cheap-relax catalog values are ~0.2% off the true
      // optima, so re-descend the best entries with the full L-BFGS (forked,
      // capped) and certify each limit with the LP. These are the remarkable
      // configurations the run returns.
      int n_polish = std::min<int>(catalog.size(), 8);
      if (n_polish > 0) printf("--- catalog (%d entries, polishing %d) ---\n",
                               (int)catalog.size(), n_polish);
      for (int kx = 0; kx < n_polish; kx++) {
        ::remove(tmpstat.c_str());
        pid_t pid = fork();
        if (pid == 0) {
          DescendResult r = Descend(catalog[kx].second, 12, std::cerr, false);
          FILE *f = fopen(tmpstat.c_str(), "w");
          if (f) {
            if (r.success) {
              fprintf(f, "OK %.15g\n", r.theta);
              fclose(f);
              WriteConfigFile(tmpconf, r.conf);
            } else {
              fprintf(f, "FAIL\n");
              fclose(f);
            }
          }
          _exit(0);
        }
        auto ct0 = std::chrono::steady_clock::now();
        bool fin = false;
        while (true) {
          int status;
          if (waitpid(pid, &status, WNOHANG) == pid) { fin = true; break; }
          double el = std::chrono::duration_cast<std::chrono::duration<double>>(
                          std::chrono::steady_clock::now() - ct0).count();
          if (el > 420.0) break;
          usleep(200000);
        }
        if (!fin) {
          kill(pid, SIGKILL);
          int status;
          waitpid(pid, &status, 0);
          printf("catalog %d: theta~%.7f polish TIMEOUT\n", kx,
                 catalog[kx].first);
          fflush(stdout);
          continue;
        }
        FILE *f = fopen(tmpstat.c_str(), "r");
        char tag[16] = {0};
        double th = 0;
        if (!f || fscanf(f, "%15s %lf", tag, &th) != 2 ||
            std::string(tag) != "OK") {
          if (f) fclose(f);
          printf("catalog %d: theta~%.7f polish failed\n", kx,
                 catalog[kx].first);
          fflush(stdout);
          continue;
        }
        fclose(f);
        PeriodicConfig pc2 = ReadConfigFile(tmpconf);
        // certify (inline: one tessellation, JointGrad, the small-box LP)
        int rg = -1, act = -1;
        double rate = -1;
        try {
          DensityResult drx = CoveringDensity(pc2);
          auto basis = SymBasis(pc2.n);
          int dQx = basis.size();
          int dimx = dQx + pc2.n * (pc2.m - 1);
          std::vector<double> R2x;
          std::vector<Eigen::VectorXd> gx;
          JointGrad(drx.cells, pc2.Q, pc2.C, basis, R2x, gx);
          double mu2x = *std::max_element(R2x.begin(), R2x.end());
          Eigen::VectorXd gdet = Eigen::VectorXd::Zero(dimx);
          Eigen::MatrixXd Qinv = pc2.Q.inverse();
          for (int u2 = 0; u2 < dQx; u2++) {
            int kk = basis[u2].first, ll = basis[u2].second;
            gdet(u2) = (kk == ll) ? Qinv(kk, kk) : 2.0 * Qinv(kk, ll);
          }
          double rho_cert = 1e-3;
          DirLPResult lp = DirectionLP(R2x, gx, gdet, mu2x, pc2.n, rho_cert);
          rate = -lp.predicted / rho_cert;
          act = lp.n_active;
          rg = (lp.status == 0 && rate < 1e-7) ? 1 : 0;
        } catch (std::runtime_error &e) {
        }
        char fn[4096];
        snprintf(fn, sizeof(fn), "%s.opt.%d", argv[3], kx);
        WriteConfigFile(fn, pc2);
        printf("catalog %d: theta=%.13f rigid=%d active=%d rate=%.2e -> %s\n",
               kx, th, rg, act, rate, fn);
        fflush(stdout);
      }
      return 0;
    }
    if (mode == "addcoset" && argc == 4) {
      // Structured growth m -> m+1: add one coset at the deepest hole (the
      // worst-covered point) of the input configuration, with a small jitter
      // to stay off the exactly-cospherical locus. Placing a point where the
      // covering radius is attained is the natural way to seed a richer coset
      // set from a good one, instead of a random start.
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int n = conf.n, m = conf.m;
      DensityResult dr = CoveringDensity(conf);
      Eigen::VectorXd z = DeepestHole(dr.cells, conf.Q, conf.C);
      std::mt19937_64 rng(1234567);
      std::normal_distribution<double> gauss(0.0, 1.0);
      PeriodicConfig nc;
      nc.n = n;
      nc.m = m + 1;
      nc.Q = conf.Q;
      nc.C.resize(m + 1, n);
      nc.C.topRows(m) = conf.C;
      for (int j = 0; j < n; j++)
        nc.C(m, j) = (z(j) - std::floor(z(j))) + 0.02 * gauss(rng);
      WriteConfigFile(argv[3], nc);
      printf("added coset at deepest hole: m=%d -> m=%d\n", m, m + 1);
      return 0;
    }
    if (mode == "evaluate-pc" && argc == 3) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      PCResult pc = PackingCovering(conf);
      std::cout << "n=" << conf.n << " m=" << conf.m
                << " classes=" << pc.cells.size() << "\n";
      printf("mu2=%.13f lambda2=%.13f\n", pc.mu2, pc.lambda2);
      printf("gamma=%.13f\n", pc.gamma);
      return 0;
    }
    if (mode == "descend-pc" && (argc == 4 || argc == 5)) {
      PeriodicConfig conf = ReadConfigFile(argv[2]);
      int rounds = argc == 5 ? atoi(argv[4]) : 15;
      DescendResult res = DescendPC(conf, rounds, std::cerr, true);
      if (!res.success) {
        std::cerr << "the descent could not evaluate any configuration\n";
        return 1;
      }
      printf("gamma=%.13f\n", res.theta);
      WriteConfigFile(argv[3], res.conf);
      return 0;
    }
    if (mode == "basinhop-pc" && (argc >= 5 && argc <= 7)) {
      // Basin hopping on the packing-covering ratio gamma. Same walk as
      // basinhop, with DescendPC as the (forked, time-capped) relaxation.
      PeriodicConfig seed_conf = ReadConfigFile(argv[2]);
      double max_seconds = atof(argv[4]);
      unsigned seed = argc >= 6 ? unsigned(atol(argv[5])) : 1u;
      // collect every distinct local optimum below this threshold; after the
      // walk the best entries are re-polished by a deeper DescendPC and
      // reported with their active counts -- the remarkable configurations
      double collect_below = argc == 7 ? atof(argv[6]) : 0.0;
      std::vector<std::pair<double, PeriodicConfig>> catalog;
      auto collect = [&](double gm, PeriodicConfig const &c) {
        if (collect_below <= 0.0 || gm >= collect_below) return;
        for (auto &e : catalog)
          if (std::abs(e.first - gm) < 1e-7) return;
        catalog.emplace_back(gm, c);
        std::sort(catalog.begin(), catalog.end(),
                  [](auto const &a, auto const &b) { return a.first < b.first; });
        if (catalog.size() > 40) catalog.pop_back();
      };
      int n = seed_conf.n, m = seed_conf.m;
      SearchMaxPts() = 4000 * m;
      // best LATTICE packing-covering values: gamma_3 (A_3^*, optimal even
      // among non-lattices by Boeroeczky), gamma_4 (Ho_4), gamma_5 (Ho_5).
      // For n = 4, 5 the non-lattice problem is open: below the record is a
      // genuine discovery.
      double record = 0.0;
      if (n == 3) record = 1.2909944487358056;   // sqrt(5/3), A_3^*
      if (n == 4) record = 1.3625000772664266;   // sqrt(8 sqrt3 - 12), Ho_4
      if (n == 5) record = 1.4494568681327882;   // sqrt(3/2 + sqrt(13)/6), Ho_5
      std::mt19937_64 rng(seed);
      std::uniform_real_distribution<double> unif(0.0, 1.0);
      std::normal_distribution<double> gauss(0.0, 1.0);
      double start_budget = 150.0;
      std::string tmpstat = std::string(argv[3]) + ".child.stat";
      std::string tmpconf = std::string(argv[3]) + ".child.conf";
      auto relax_forked = [&](PeriodicConfig const &c, PeriodicConfig &out,
                              double &gam) -> bool {
        ::remove(tmpstat.c_str());
        pid_t pid = fork();
        if (pid == 0) {
          DescendResult r = DescendPC(c, 6, std::cerr, false);
          FILE *f = fopen(tmpstat.c_str(), "w");
          if (f) {
            if (r.success) {
              fprintf(f, "OK %.15g\n", r.theta);
              fclose(f);
              WriteConfigFile(tmpconf, r.conf);
            } else {
              fprintf(f, "FAIL\n");
              fclose(f);
            }
          }
          _exit(0);
        }
        auto ct0 = std::chrono::steady_clock::now();
        bool fin = false;
        while (true) {
          int status;
          if (waitpid(pid, &status, WNOHANG) == pid) { fin = true; break; }
          double el = std::chrono::duration_cast<std::chrono::duration<double>>(
                          std::chrono::steady_clock::now() - ct0).count();
          if (el > start_budget) break;
          usleep(200000);
        }
        if (!fin) {
          kill(pid, SIGKILL);
          int status;
          waitpid(pid, &status, 0);
          return false;
        }
        FILE *f = fopen(tmpstat.c_str(), "r");
        char tag[16] = {0};
        if (!f || fscanf(f, "%15s", tag) != 1 || std::string(tag) != "OK") {
          if (f) fclose(f);
          return false;
        }
        if (fscanf(f, "%lf", &gam) != 1) { fclose(f); return false; }
        fclose(f);
        out = ReadConfigFile(tmpconf);
        return true;
      };
      PeriodicConfig curc;
      double curg;
      bool init_ok = false;
      for (int attempt = 0; attempt < 8 && !init_ok; attempt++) {
        // a random seed may need a larger enumeration ball than a healthy
        // configuration just to be evaluated once: raise the cloud cap
        // progressively during the init attempts (the fork watchdog still
        // bounds the cost), then restore the tight cap for the walk
        SearchMaxPts() = 4000 * m * (1 + attempt);
        PeriodicConfig s2 = seed_conf;
        if (attempt > 0) {
          for (int t = 1; t < m; t++)
            for (int j = 0; j < n; j++)
              s2.C(t, j) += 0.10 * gauss(rng);
        }
        init_ok = relax_forked(s2, curc, curg);
        if (!init_ok)
          printf("init attempt %d failed, retrying with a jittered seed\n",
                 attempt);
        fflush(stdout);
      }
      if (!init_ok) {
        std::cerr << "basinhop-pc: initial relaxation failed after retries\n";
        return 1;
      }
      SearchMaxPts() = 4000 * m;
      PeriodicConfig bestc = curc;
      double bestg = curg;
      collect(curg, curc);
      WriteConfigFile(argv[3], bestc);
      printf("init: gamma=%.13f\n", curg);
      fflush(stdout);
      double T = 0.04 * curg, Tmin = 1e-4;
      auto t0 = std::chrono::steady_clock::now();
      int n_acc = 0, n_hop = 0;
      for (int hop = 0;; hop++) {
        double el = std::chrono::duration_cast<std::chrono::duration<double>>(
                        std::chrono::steady_clock::now() - t0).count();
        if (el > max_seconds) break;
        PeriodicConfig trial = curc;
        std::vector<CellClass> cells;
        try {
          cells = DelaunayCellClasses(trial);
        } catch (std::runtime_error &e) {
          continue;
        }
        int t = (m > 2) ? 1 + int(rng() % (m - 1)) : 1;
        double u = unif(rng);
        if (u < 0.15) {
          Eigen::VectorXd z = DeepestHole(cells, trial.Q, trial.C);
          for (int j = 0; j < n; j++)
            trial.C(t, j) = (z(j) - std::floor(z(j))) + 0.05 * gauss(rng);
        } else {
          double sigma = (u < 0.6) ? 0.06 : 0.12;
          for (int j = 0; j < n; j++) trial.C(t, j) += sigma * gauss(rng);
        }
        PeriodicConfig candc;
        double candg;
        n_hop++;
        if (!relax_forked(trial, candc, candg)) {
          printf("hop %d: relax timeout/fail\n", hop);
          fflush(stdout);
          continue;
        }
        collect(candg, candc);
        double dg = candg - curg;
        bool acc = (dg < 0.0) || (unif(rng) < std::exp(-dg / T));
        if (acc) { curc = candc; curg = candg; n_acc++; }
        bool rec = record > 0.0 && candg < record - 1e-9;
        if (candg < bestg - 1e-12) {
          bestg = candg;
          bestc = candc;
          WriteConfigFile(argv[3], bestc);
          if (rec) {
            char fn[4096];
            snprintf(fn, sizeof(fn), "%s.record", argv[3]);
            WriteConfigFile(fn, bestc);
          }
        }
        printf("hop %d: gamma=%.13f %s cur=%.13f best=%.13f T=%.4f%s\n", hop,
               candg, acc ? "acc" : "rej", curg, bestg, T,
               rec ? "  <<< BEATS RECORD" : "");
        fflush(stdout);
        T = std::max(Tmin, T * 0.995);
      }
      printf("best=%.13f accept_rate=%.2f (record=%.13f)\n", bestg,
             n_hop ? double(n_acc) / n_hop : 0.0, record);
      WriteConfigFile(argv[3], bestc);
      // Post-polish: re-descend the best catalog entries with a deeper
      // DescendPC and report the active counts (Delaunay orbits attaining the
      // covering radius, pair vectors attaining the minimal distance), the
      // data that identifies a remarkable configuration.
      int n_polish = std::min<int>(catalog.size(), 8);
      if (n_polish > 0) printf("--- catalog (%d entries, polishing %d) ---\n",
                               (int)catalog.size(), n_polish);
      for (int kx = 0; kx < n_polish; kx++) {
        ::remove(tmpstat.c_str());
        pid_t pid = fork();
        if (pid == 0) {
          DescendResult r = DescendPC(catalog[kx].second, 15, std::cerr, false);
          FILE *f = fopen(tmpstat.c_str(), "w");
          if (f) {
            if (r.success) {
              fprintf(f, "OK %.15g\n", r.theta);
              fclose(f);
              WriteConfigFile(tmpconf, r.conf);
            } else {
              fprintf(f, "FAIL\n");
              fclose(f);
            }
          }
          _exit(0);
        }
        auto ct0 = std::chrono::steady_clock::now();
        bool fin = false;
        while (true) {
          int status;
          if (waitpid(pid, &status, WNOHANG) == pid) { fin = true; break; }
          double el = std::chrono::duration_cast<std::chrono::duration<double>>(
                          std::chrono::steady_clock::now() - ct0).count();
          if (el > 600.0) break;
          usleep(200000);
        }
        if (!fin) {
          kill(pid, SIGKILL);
          int status;
          waitpid(pid, &status, 0);
          printf("catalog %d: gamma~%.7f polish TIMEOUT\n", kx,
                 catalog[kx].first);
          fflush(stdout);
          continue;
        }
        FILE *f = fopen(tmpstat.c_str(), "r");
        char tag[16] = {0};
        double gm = 0;
        if (!f || fscanf(f, "%15s %lf", tag, &gm) != 2 ||
            std::string(tag) != "OK") {
          if (f) fclose(f);
          printf("catalog %d: gamma~%.7f polish failed\n", kx,
                 catalog[kx].first);
          fflush(stdout);
          continue;
        }
        fclose(f);
        PeriodicConfig pc2 = ReadConfigFile(tmpconf);
        int amu = -1, alam = -1;
        int n_orb = -1;
        try {
          PCResult fr = PackingCovering(pc2);
          n_orb = fr.cells.size();
          amu = 0;
          for (auto const &cl : fr.cells) {
            Eigen::MatrixXd P = CellPositions(cl, pc2.C);
            Eigen::VectorXd ctr;
            double R2;
            CircumcenterRadius2(P, pc2.Q, ctr, R2);
            if (R2 > fr.mu2 * (1 - 1e-7)) amu++;
          }
          alam = ShortPairVectors(pc2.Q, pc2.C, fr.lambda2 * (1 + 1e-7)).size();
        } catch (std::runtime_error &e) {
        }
        char fn[4096];
        snprintf(fn, sizeof(fn), "%s.opt.%d", argv[3], kx);
        WriteConfigFile(fn, pc2);
        printf("catalog %d: gamma=%.13f mu-active=%d/%d lambda-active=%d -> %s\n",
               kx, gm, amu, n_orb, alam, fn);
        fflush(stdout);
      }
      return 0;
    }
    std::cerr << "unrecognized arguments; run without arguments for usage\n";
    return -1;
  } catch (std::exception const &e) {
    std::cerr << "Error in PERIODIC_JointCoveringSearch: " << e.what() << "\n";
    return 1;
  }
}
