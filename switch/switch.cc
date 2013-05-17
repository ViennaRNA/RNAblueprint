/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* Compile with: g++ -std=c++11 -g -lboost_program_options -o struct2graph struct2graph.cc
*/

// include header
#include "switch.h"

// overload << operator to print vectors with any content
template <typename T>
std::ostream& operator<< (std::ostream& os, std::vector<T>& vec) {
	int i = 0;
	for (auto elem : vec) {
		os << "(" << i++ << ") " << elem << std::endl;
	}
	return os;
}

bool verbose = false;
std::string outfile = "";


//! main program starts here
int main(int ac, char* av[]) {

	// variables
	std::ostream* out = &std::cout;			// out stream

	// initialize command line options
	boost::program_options::variables_map vm = init_options(ac, av);
	
	// input handling ( we read from std:in per default and switch to a file if it is given in the --in option
	// std::in will provide a pseudo interface to enter structures directly!
	std::vector<std::string> input;
	
	if (vm.count("in")) {
		if (verbose) { std::cerr << "will read graphml file given in the options." << std::endl; }
			std::ifstream* infile = new std::ifstream(vm["in"].as<std::string>(), std::ifstream::in);
		if (infile->is_open()) {
			read_input(infile, input);
			infile->close();
		} else {
			std::cerr << "Unable to open file";
			return 1;
		}
	} else {
		std::cerr << "Input structures (one per line); @ to quit" << std::endl
		<< "....,....1....,....2....,....3....,....4....,....5....,....6....,....7....,....8" << std::endl;
		read_input(&std::cin, input);		// read infile into array
	}
	
	for (auto line : input) {
		std::vector<std::string> structures;
		std::regex e ("[\(\)\.]");
		if (line.at(0) == '>') {
			structures.clear();
			std::string name = line.substr(1);
		} else if (std::regex_match(line, e)) {
			structures.push_back(line);
			if (structures.size() > 2) {
				std::cerr << "Too many structures in input!" << std::cerr;
				exit(1);
			} else if (structures.size() == 2) {
				// we got two structures -> start the engines!
				parse_structures(structures);
			}
		} // else if match nucleotid alphabet -> get constraints TODO
	}
	
	
	
	
	*out << input;
	*out << "done!";

	return 0;
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
		if (verbose) { std::cerr << "output will be written to file." << std::endl; }
	}
	
	return vm;
}

void read_input(std::istream* in, std::vector<std::string>& input) {
	// read input file
	std::string line;
	while (!in->eof()) {
		getline(*in,line);
		if (line == "@") {
			std::fclose(stdin);
		} else if (line.length() != 0) {
			input.push_back(line);
		}
	}
	
	// exit if there is no input
	if (input.empty()) { exit(1); }
	if (verbose) { std::cerr << "Read following input:" << std::endl; }
}

void parse_structures(std::vector<std::string> structures) {
	
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
	
	// count the number of positons
	int num_vertices = structures[0].length();
	if (verbose) { std::cerr << "Generating Graph with " << num_vertices << " vertices." << std::endl; }
	//Graph g(num_vertices);
	
	// give the vertices names
	//int vertex_name = 0;
	//Graph::vertex_iterator v, v_end;
	//for (boost::tie(v,v_end) = boost::vertices(g); v != v_end; ++v) {
	//	boost::put(boost::vertex_color_t(), g, *v, vertex_name++);
	//}

	// iterate over structures from input
	for (auto elem : structures) {
		std::vector<int> pair_table(structures[0].length(), 0);		// remember position of the open bracket
		unsigned int pos = 0;							// position of the character in the structure
		unsigned int open = 0;							// remember how many open brackets there are
		// iterate over characters from structure
		while (pos < elem.length()) {
			if (elem[pos] == '(') {
				pair_table[open] = pos;
				if (verbose) { std::cerr << elem[pos] << ", open count: "<< open; }
				open++;
			} else if (elem[pos] == ')') {
				open--;
				// check if edge already exists
				//bool exists_ab = boost::edge(boost::vertex(pair_table[open],g), boost::vertex(pos,g), g).second;
				//bool exists_ba = boost::edge(boost::vertex(pos,g), boost::vertex(pair_table[open],g), g).second;
				//if (!exists_ab && !exists_ba) {
					// add edge
				//	boost::add_edge(boost::vertex(pair_table[open],g), boost::vertex(pos,g), g);
				//}
				// reset value
				pair_table[open] = pos;
				if (verbose) { std::cerr << elem[pos] << ", open count: "<< open; }
			} else if (elem[pos] != '.') {
				std::cerr << std::endl << "Unknown character in dot bracked representation" << std::endl;
				exit(1);
			}
			// error handling: there can't be more closing brackets than opening ones
			if (open < 0) {
				std::cerr << std::endl << "Unbalanced brackets in make_pair_table" << std::endl;
				exit(1);
			}
			if (verbose) { std::cerr  << " pos count:" << pos << std::endl; }
			pos++;
		}
		// error handling: at the end all brackets must be closed again!
		if (open != 0) {
			std::cerr << std::endl << "too few closed brackets in make_pair_table" << std::endl;
			exit(1);
		}
	}
}
