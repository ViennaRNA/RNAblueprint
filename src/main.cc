/* RNAdesign
 * A program for designing RNA molecules.
 *
 * Created on: 25.03.2013
 * Author: Stefan Hammer <s.hammer@univie.ac.at>
 * License: GPLv3
 *
 */

// include header
#include "main.h"
#include "energy.h"

// boost components
#include <boost/graph/iteration_macros.hpp>

//declare global variables
bool debug = false;
std::mt19937 rand_gen;

//! main program starts here

int main (int ac, char* av[]) {

  // initialize command line options
  boost::program_options::variables_map vm = init_options(ac, av);
  int num_trees = 0;
  unsigned int number_of_designs = 4;
  if (vm.count("stat-trees")) {
    num_trees = vm["stat-trees"].as<int>();
  }
  if (vm.count("num")) {
    number_of_designs = vm["num"].as<unsigned int>();
  }
  bool verbose = vm["verbose"].as<bool>();

  // initialize mersenne twister with our seed
  unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();
  if (vm.count("seed")) {
    seed = vm["seed"].as<unsigned long>();
  }
  rand_gen.seed(seed);
  if (verbose) {
    std::cerr << "Using this seed: " << seed << std::endl;
  }

  // input handling ( we read from std::in per default and switch to a file if it is given in the --in option
  // std::in will provide a pseudo interface to enter structures directly!
  std::vector<std::string> structures;

  if (vm.count("in")) {
    std::ifstream* infile = new std::ifstream(vm["in"].as<std::string>(), std::ifstream::in);
    if (infile->is_open()) {
      structures = read_input(infile);
      infile->close();
      std::cerr << "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8" << std::endl;
    } else {
      std::cerr << "Unable to open file";
      return 1;
    }
  } else {
    std::cerr << "Input structures (one per line); @ to quit" << std::endl
        << "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8" << std::endl;
    structures = read_input(&std::cin); // read infile into array
  }

  std::ostream* out = &std::cout;
  if (vm.count("out")) {
    out = new std::ofstream(vm["out"].as<std::string>(), std::ofstream::out);
  }

  design::DependencyGraph dependency_graph(structures);

  if (!dependency_graph.is_bipartite()) {
    std::cerr << "Impossible to find a solution for this input: Dependency graph is not bipartite!" << std::endl;
    exit(1);
  }

  if (debug) {
    std::cerr << "Size of solution space: " << dependency_graph.number_of_sequences() << std::endl;
  }

  while (number_of_designs > 0) {
    dependency_graph.mutate(); // color the graph and get the sequence
    std::string sequence = dependency_graph.get_sequence();
    *out << sequence << std::endl;
    for (auto s : structures) {
      *out << s << "\t" << energy_of_structure(sequence, s) << std::endl; // calculate the energies
    }
    *out << "mfe: " << std::endl;
    std::string mfe_structure;
    float mfe = fold(sequence, mfe_structure);
    *out << mfe_structure << "\t" << mfe << std::endl << std::endl;
    number_of_designs--;
  }

  return 0;
}

std::vector<std::string> read_input (std::istream* in) {
  // read input file
  std::string line;
  std::vector<std::string> structures;
  while (!in->eof()) {
    getline(*in, line);
    if (line == "@") {
      std::fclose(stdin);
    } else if (line.length() != 0) {
      structures.push_back(line);
    }
  }

  // exit if there is no input
  if (structures.empty()) {
    exit(1);
  }

  if (debug) {
    std::cerr << "Read following structures:" << std::endl;
  }
  // check if structures have equeal length
  unsigned int length = 0;
  for (auto elem : structures) {
    if (debug) {
      std::cerr << elem << std::endl;
    }
    if ((length != elem.length()) && (length != 0)) {
      std::cerr << "Structures have unequal length." << std::endl;
      exit(1);
    }
    length = elem.length();
  }
  return structures;
}

boost::program_options::variables_map init_options (int ac, char* av[]) {
  // boost option parser
  // http://www.boost.org/doc/libs/1_53_0/doc/html/program_options/tutorial.html
  namespace po = boost::program_options;
  // Group of options that will be allowed only on command line
  po::options_description generic("Generic options");
  generic.add_options()
      ("help,h", "print help message")
      ("verbose,v", po::bool_switch()->default_value(false)->zero_tokens(), "be verbose")
      ("debug,d", po::value(&debug)->zero_tokens(), "be verbose for debugging")
      ;

  // Group of options that will be allowed on command line and in a config file
  po::options_description config("Program options");
  config.add_options()
      ("in,i", po::value<std::string>(), "input file which contains the structures [string]")
      ("out,o", po::value<std::string>(), "output file which contains the sequences [string]")
      ("seed,s", po::value<unsigned long>(), "random number generator seed [unsigned long]")
      ("num,n", po::value<unsigned int>(), "number of designs (default: 4) [unsigned int]")
      ("stat-trees,t", po::value<int>(), "only do ear-decomposition statistics: define amount of different spanning trees for every root to calculate [int]")
      ;

  po::positional_options_description p;
  p.add("in", 1).add("out", 2);

  po::options_description cmdline_options;
  cmdline_options.add(generic).add(config);

  po::variables_map vm;
  po::store(po::command_line_parser(ac, av).options(cmdline_options).positional(p).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << cmdline_options << "\n";
    exit(1);
  }

  return vm;
}

