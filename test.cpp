#include <algorithm>
#include <vector>

#include "distance_algorithms/transposition3.hpp"
#include "distance_algorithms/reversal4.hpp"
#include "distance_algorithms/reversal_transposition45.hpp"
#include "cycle/cycles.hpp"
#include "heur/ga.hpp"
#include "misc/genome.hpp"
#include "misc/io.hpp"
#include "misc/permutation.hpp"
#include "misc/timer.hpp"
#include "quickcheck/quickcheck/quickcheck.hh"
using namespace quickcheck;

struct TestPerm {
  unique_ptr<Permutation> perm;
  TestPerm() {}
  TestPerm(const TestPerm &perm_t) {
    perm = unique_ptr<Permutation>(new Permutation(*perm_t.perm));
  }
};

ostream &operator<<(std::ostream &os, const TestPerm &t_perm) {
  os << *t_perm.perm;
  return os;
}

struct TestCycleGraph {
  unique_ptr<CycleGraph> cg;
  TestCycleGraph() {}
  TestCycleGraph(const TestCycleGraph &cg_t) {
    cg = unique_ptr<CycleGraph>(new CycleGraph(*cg_t.cg));
  }
};

ostream &operator<<(std::ostream &os, const TestCycleGraph &t_cg) {
  os << *t_cg.cg;
  return os;
}

struct Genomes {
  unique_ptr<Genome> g;
  unique_ptr<Genome> h;
};

Genomes generate_genomes(size_t n, bool add_sign, bool add_indels) {
  string sg, ig, sh, ih;
  int ir;
  int sum = 0;
  bool positive;
  Genomes gs;

  ir = generateInRange(50, 100);
  ig.append(to_string(ir));
  ig.append(" ");
  sum += ir;
  for (size_t i = 0; i < n; i++) {
    ir = generateInRange(0, 100);
    generate(n, positive);
    if (!add_sign) {
      sg.append("1");
    } else if (positive) {
      sg.append("+1");
    } else {
      sg.append("-1");
    }
    sg.append(" ");
    ig.append(to_string(ir));
    ig.append(" ");
    generate(n, positive);
    if (!add_sign) {
      sh.append("1");
    } else if (positive) {
      sh.append("+1");
    } else {
      sh.append("-1");
    }
    sh.append(" ");
    sum += ir;
  }

  vector<int> ir_list;
  ir_list.push_back(0);
  for (size_t i = 0; i < n; i++) {
    ir_list.push_back(generateInRange(0, sum));
  }
  ir_list.push_back(sum);
  sort(ir_list.begin(), ir_list.end());
  for (size_t i = 1; i <= n + 1; i++) {
    ir = ir_list[i] - ir_list[i - 1];
    ih.append(to_string(ir));
    ih.append(" ");
  }

  gs.g = unique_ptr<Genome>(new Genome(sg, ig, true));
  gs.h = unique_ptr<Genome>(new Genome(sh, ih, true));

  if (add_indels) {
      for (int i = 2; i < int(gs.g->size()); i++) {
        if(rand() % 100 < 10) {
            IR x = rand() % (gs.g->get_ir(i - 1) + gs.g->get_ir(i) + 1);
            gs.g->deletion(i,i+1,x);
            i--;
        }
      }

      for (int i = 2; i < int(gs.h->size()); i++) {
        if(rand() % 100 < 10) {
            IR x = rand() % (gs.h->get_ir(i - 1) + gs.h->get_ir(i) + 1);
            gs.h->deletion(i,i+1,x);
            i--;
        }
      }
  }
  /* cout << "g: " << *gs.g << endl; */
  /* cout << "h: " << *gs.h << endl; */

  return gs;
}

void generate(size_t n, TestPerm &t_perm) {
  bool add_sign;
  generate(n, add_sign);
  Genomes gs = generate_genomes(n, add_sign, false);
  t_perm.perm =
      unique_ptr<Permutation>(new Permutation(*gs.g, *gs.h, add_sign));
}

void generate(size_t n, TestCycleGraph &t_cg) {
  Genomes gs = generate_genomes(n, true, true);
  t_cg.cg = unique_ptr<CycleGraph>(new CycleGraph(*gs.g, *gs.h));
}

class PBounds3T : public Property<TestPerm> {
  bool holdsFor(const TestPerm &t_perm) {
    Transposition3 alg;

    int breaks_count = 0;
    for (size_t b = 1; b <= t_perm.perm->size() - 1; b++) {
      if (t_perm.perm->breakpoint((*t_perm.perm)[b], false))
        breaks_count++;
    }

    int dist = alg.estimate_distance(*t_perm.perm);

    bool sucess = true;
    sucess = sucess && dist >= (breaks_count / 3);
    sucess = sucess && dist <= (4 * breaks_count / 3);
    if (!sucess) {
      cout << "pi: " << *t_perm.perm << endl;
      cout << "dist: " << dist << endl;
      cout << "b(pi):" << breaks_count << endl;
    } else {
      cout << ".";
      cout.flush();
    }
    return sucess;
  }
};

class PBounds4R : public Property<TestPerm> {
  bool holdsFor(const TestPerm &t_perm) {
    Reversal4 alg;

    int breaks_count = 0;
    for (size_t b = 1; b <= t_perm.perm->size() - 1; b++) {
      if (t_perm.perm->breakpoint((*t_perm.perm)[b], true))
        breaks_count++;
    }

    int dist = alg.estimate_distance(*t_perm.perm);

    bool sucess = true;
    sucess = sucess && dist >= (breaks_count / 2);
    sucess = sucess && dist <= (4 * breaks_count / 2);
    if (!sucess) {
      cout << "pi: " << *t_perm.perm << endl;
      cout << "dist: " << dist << endl;
      cout << "b(pi):" << breaks_count << endl;
    } else {
      cout << ".";
      cout.flush();
    }
    return sucess;
  }
};

class PBounds45RT : public Property<TestPerm> {
  bool holdsFor(const TestPerm &t_perm) {
    ReversalTransposition45 alg;

    int breaks_count = 0;
    for (size_t b = 1; b <= t_perm.perm->size() - 1; b++) {
      if (t_perm.perm->breakpoint((*t_perm.perm)[b], true))
        breaks_count++;
    }

    int dist = alg.estimate_distance(*t_perm.perm);

    bool sucess = true;
    sucess = sucess && dist >= (breaks_count / 3);
    sucess = sucess && dist <= (4.5 * breaks_count / 3);
    if (!sucess) {
      cout << "pi: " << *t_perm.perm << endl;
      cout << "dist: " << dist << endl;
      cout << "b(pi):" << breaks_count << endl;
    } else {
      cout << ".";
      cout.flush();
    }
    return sucess;
  }
};

class GABounds : public Property<TestCycleGraph> {
  bool holdsFor(const TestCycleGraph &t_cg) {
    Timer timer;
    ostream *out = new ostream(0);

    GA ga = GA(new Chromossome(*t_cg.cg), 0.5, 0.5, 2, 2, 2, 3, out, timer);
    ga.solve(timer);
    vector<int> obj = ga.get_best_obj();

    bool sucess = true;
    sucess = sucess && obj[0] >= 0;
    sucess = sucess && obj[1] >= 1;
    sucess = sucess && obj[0] <= int(t_cg.cg->size()) / 4;
    sucess = sucess && obj[1] <= int(t_cg.cg->size()) / 4;
    if (!sucess) {
      cout << "cg: " << *t_cg.cg << endl;
      cout << "cg_size: " << t_cg.cg->size() << endl;
      cout << "obj: ";
      print_vec(cout, obj);
    } else {
      cout << ".";
      cout.flush();
    }
    return sucess;
  }
};

int main() {
  // set seed
  int seed = time(nullptr);
  cout << "seed: " << seed << endl;
  srand(seed);

  check<PBounds3T>(
      "upper and lower bounds hold for factor 4 algorithm for transposition");
  check<PBounds4R>(
      "upper and lower bounds hold for factor 4 algorithm for reversal");
  check<PBounds45RT>(
      "upper and lower bounds hold for factor 4.5 algorithm for reversal and transposition");
  check<GABounds>(
      "upper and lower bounds hold for cycle decomposition with ga");
  return 0;
}
