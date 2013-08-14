/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "main.h"
#include "parsestruct.h"
#include "printgraph.h"
#include "decompose.h"
#include "pathcoloring.h"

//declare global variables
bool verbose = false;
std::string outfile;

//! main program starts here
int main(int ac, char* av[]) {

	// initialize command line options
	boost::program_options::variables_map vm = init_options(ac, av);
	int num_trees = 0;
	if (vm.count("stat-trees")) { num_trees = vm["stat-trees"].as<int>(); }
	bool ramachandran = vm["ramachandran"].as<bool>();
	bool no_bipartite_check = vm["noBipartiteCheck"].as<bool>();
	
	// initialize mersenne twister with our seed
	unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();
	if (vm.count("seed")) {
		seed = vm["seed"].as<unsigned long>();
	}
	rand_gen.seed(seed);
	std::cerr << "Using this seed: " << seed << std::endl;
	
	// input handling ( we read from std:in per default and switch to a file if it is given in the --in option
	// std::in will provide a pseudo interface to enter structures directly!
	std::vector<std::string> structures;
	
	if (vm.count("in")) {
		if (verbose) { std::cerr << "will read graphml file given in the options." << std::endl; }
			std::ifstream* infile = new std::ifstream(vm["in"].as<std::string>(), std::ifstream::in);
		if (infile->is_open()) {
			structures = read_input(infile);
			infile->close();
		} else {
			std::cerr << "Unable to open file";
			return 1;
		}
	} else {
		std::cerr << "Input structures (one per line); @ to quit" << std::endl
		<< "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8" << std::endl;
		structures = read_input(&std::cin);			// read infile into array
	}
	
	std::ostream* out = &std::cout;					// out stream
	
	// variables
	Graph graph = parse_structures(structures);			// generate graph from input vector
	*out << "dependency graph:";
	print_graph(graph, out, "root-graph");				// print the graph as GML to a ostream
	
	decompose_graph(graph, out, num_trees, 				// decompose the graph into its connected components, biconnected
		ramachandran, no_bipartite_check);			// components and decompose blocks via ear decomposition
	
	std::string sequence;
	*out << "Number of sequences: " << generate_path_seq (sequence, 'A', 'U', 10) << std::endl;
	*out << sequence << std::endl;
	
	return 0;
}

std::vector<std::string> read_input(std::istream* in) {
	// read input file
	std::string line;
	std::vector<std::string> structures;
	while (!in->eof()) {
		getline(*in,line);
		if (line == "@") {
			std::fclose(stdin);
		} else if (line.length() != 0) {
			structures.push_back(line);
		}
	}
	
	// exit if there is no input
	if (structures.empty()) { exit(1); }

	if (verbose) { std::cerr << "Read following structures:" << std::endl; }
	// check if structures have equeal length
	unsigned int length = 0;
	for (auto elem : structures) {
		if (verbose) { std::cerr << elem << std::endl; }
		if ((length != elem.length()) && (length != 0)){
			std::cerr << "Structures have unequal length." << std::endl;
			exit(1);
		}
		length = elem.length();
	}
	return structures;
}


boost::program_options::variables_map init_options(int ac, char* av[]) {
	// boost option parser
	// http://www.boost.org/doc/libs/1_53_0/doc/html/program_options/tutorial.html
	namespace po = boost::program_options;
	// Group of options that will be allowed only on command line
	po::options_description generic("Generic options");
	generic.add_options()
		("help,h", "print help message")
		("verbose,v", po::value(&verbose)->zero_tokens(), "be verbose")
	;
	
	// Group of options that will be allowed on command line and in a config file
	po::options_description config("Program options");
	config.add_options()
		("in,i", po::value<std::string>(), "input file which contains the structures [string]")
		("out,o", po::value<std::string>(&outfile), "write all (sub)graphs to gml files starting with given name [string]")
		("seed,s", po::value<unsigned long>(), "random number generator seed [unsigned long]")
		("ramachandran,r", po::bool_switch()->default_value(false)->zero_tokens(), "Use the Ramachandran ear decomposition algorithmus")
		("stat-trees,t", po::value<int>(), "only do ear-decomposition statistics: define amount of different spanning trees for every root to calculate [int]")
		("noBipartiteCheck,b", po::bool_switch()->default_value(false)->zero_tokens(), "Don't check if input dependency graph is bipartite")
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
	if (vm.count("out")) {
		if (verbose) { std::cerr << "graphml files will be written to file." << std::endl; }
	}
	
	return vm;
}

