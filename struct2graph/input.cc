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
#include "input.h"



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

