#include "distance_algorithms/reversal4.hpp"
#include "distance_algorithms/r_or_rt_noir.hpp"
#include "distance_algorithms/transposition3.hpp"
#include "distance_algorithms/reversal_transposition45.hpp"
#include "external/external.hpp"
#include "misc/genome.hpp"
#include "misc/io.hpp"
#include "misc/permutation.hpp"
#include "misc/timer.hpp"
#include <experimental/filesystem>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>
namespace fs = experimental::filesystem;
using namespace std;

#define N_POS_ARGS 1

struct Args {
  string input_file;
  string output_folder;
  int iterations = 1;
  bool extend = false;
  bool duplicate = false;
  bool fill_zero = false;
  string alg;
};

void help(char *name) {
  cout << "usage: Calculate distances for pairs of strings. If replicas are "
       << "presents multiple random mappings are generated." << endl
       << "\t" << name << " ALG [OPTIONS]" << endl
       << endl
       << "positional arguments:" << endl
       << "\tALG      the algorithm to use, the options are:" << endl
       << "\t\t - transposition3: Asymptotic 3-approximation for sorting by intergenic transpositions in unsigned permutations" << endl
       << "\t\t - reversal4: 4-approximation for sorting by intergenic reversals in unsigned permutations" << endl
       << "\t\t - reversal_transposition45: 4.5-approximation for sorting by intergenic reversals and transpositions in unsigned permutations" << endl
       << "\t\t - reversal_noir: 2-approximation for sorting by reversals and indels in signed permutations" << endl
       << "\t\t - reversal_transposition_noir: heuristic for sorting by reversals, transposition and indels in signed permutations" << endl
       << "\t\t - the name of an executable in the external folder"
       << endl
       << endl
       << "optional arguments:" << endl
       << "\t-h, --help              show this help message and exit" << endl
       << "\t-i, --input INPUT       input file (if not provided stdin is "
          "used). Each 4 lines of the input file correspond to a instance, "
          "each line has a list of space separated values, and represent in "
          "order the origin string, the origin intergenic region list, the "
          "target string, and the target intergenic region list."
       << endl
       << "\t-o, --output OUTPUT     output folder (if not provided stdout is used)"
       << endl
       << "\t-k, --iterations ITER   number of iterations (default 1)" << endl
       << endl
       << "\t-e, --extend            whether to extend the genomes before apply the algorithm"
       << endl
       << "\t-d, --duplicate         whether to duplicate sign genes to turn them into unsigned" << endl
       << "\t-z, --zeros             whether to fill intergenic regions with zeros instead of read it (necessary for reversal_noir and reversal_transposition_noir, which do not deal with intergenic regions)"
       << endl;

  exit(EXIT_SUCCESS);
}

void get_args(Args &args, int argc, char *argv[]) {
  extern char *optarg;
  extern int optind;
  int n_pos_args = 0;

  struct option longopts[] = {
      {"input", 1, NULL, 'i'},      {"output", 1, NULL, 'o'},
      {"iterations", 1, NULL, 'k'}, {"extend", 0, NULL, 'e'},
      {"duplicate", 0, NULL, 'd'},  {"zero", 0, NULL, 'z'},
      {"help", 0, NULL, 'h'}};

  char op;
  while ((op = getopt_long(argc, argv, "i:o:k:hedz", longopts, NULL)) != -1) {
    switch (op) {
    case 'i':
      args.input_file = optarg;
      break;
    case 'o':
      args.output_folder = optarg;
      break;
    case 'k':
      args.iterations = atoi(optarg);
      break;
    case 'e':
      args.extend = true;
      break;
    case 'd':
      args.duplicate = true;
      break;
    case 'z':
      args.fill_zero = true;
      break;
    default:
      help(argv[0]);
    }
  }
  for (int i = optind; i < argc; i++) {
    args.alg = argv[i];
    n_pos_args++;
  }

  if (n_pos_args != N_POS_ARGS) {
    help(argv[0]);
  }
}

int main(int argc, char *argv[]) {
  Args args;
  ifstream is;
  unique_ptr<vector<string>> input_lines;
  unique_ptr<DistAlg> alg;

  get_args(args, argc, argv);

  // set seed
  // srand(1);
  srand(time(0));

  // set algorithm
  if (args.alg == "transposition3") {
    alg.reset(new Transposition3());
  } else if (args.alg == "reversal4") {
    alg.reset(new Reversal4());
  } else if (args.alg == "reversal_noir") {
    alg.reset(new ReversalNOIR());
  } else if (args.alg == "reversal_transposition_noir") {
    alg.reset(new ReversalTranspositionNOIR());
  } else if (args.alg == "reversal_transposition45") {
    alg.reset(new ReversalTransposition45());
  } else {
    alg.reset(new ExternalDistAlg("external/" + args.alg));
  }

  unique_ptr<ReversalNOIR> alg_rnoir = unique_ptr<ReversalNOIR>(new ReversalNOIR());

  if (args.input_file != "") {
    is.open(args.input_file);
    input_lines.reset(read_lines(is));
    is.close();
  } else {
    input_lines.reset(read_lines(cin));
  }

  int div = (args.fill_zero) ? 2 : 4;
  try {
    if (input_lines->size() % div == 1) {
      throw invalid_argument("Number of lines is not multiple of 4 or 2 (if -z).");
    }
#pragma omp parallel for
    for (size_t i = 0; i < input_lines->size(); i += div) {
      Timer timer;
      ofstream os;
      unique_ptr<Permutation> pi, pi_best;
      int dist, dist_best = std::numeric_limits<int>::max();

      InputData data;
      if (args.fill_zero) {
          data = input((*input_lines)[i], (*input_lines)[i + 1], args.extend);
      } else {
          data = input((*input_lines)[i], (*input_lines)[i + 1], (*input_lines)[i + 2], (*input_lines)[i + 3], args.extend);
      }

      if (args.output_folder != "") {
        os.open((args.output_folder / fs::path(args.input_file).filename()).string() +
                string(5 - to_string(i / div).size(), '0') + to_string(i / div) +
                "-all");
      }

      for (int j = 1; j <= args.iterations; ++j) {
        timer.mark_time();
        pi.reset(new Permutation(*data.g, *data.h, args.duplicate));
        dist = alg->estimate_distance(*pi);
        /* InputData out_data = alg_rnoir->make_genomes(*pi); */
        /* output((args.output_folder != "") ? os : cout, out_data); */
        output((args.output_folder != "") ? os : cout, dist,
               timer.since_last_mark());
        if (dist < dist_best) {
          dist_best = dist;
          pi_best = move(pi);
        }
      }

      if (args.output_folder != "") {
        os.close();
        os.open((args.output_folder / fs::path(args.input_file).filename()).string() +
                string(5 - to_string(i / div).size(), '0') + to_string(i / div) +
                "-best");
      }

      if (args.output_folder != "") {
        os << *pi_best << endl;
        output(os, dist_best, timer.elapsed_time());
      } else {
        cout << *pi_best << endl;
        output(cout, dist_best, timer.elapsed_time());
      }
    }

  } catch (const invalid_argument &e) {
    cerr << "Something went wrong!!!" << endl;
    cerr << e.what() << endl;
  }

  return 0;
}
