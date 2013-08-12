/* This program reads secundary RNA structures in dot-bracket and
* builds a graph for a latter ear-decomposition and bipartitness-check
*
* Created on: 25.03.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
*/

// include header
#include "parsestruct.h"

Graph parse_structures(std::vector<std::string> structures) {
	
	// count the number of positons
	int num_vertices = structures[0].length();
	if (verbose) { std::cerr << "Generating Graph with " << num_vertices << " vertices." << std::endl; }
	Graph g(num_vertices);
	
	// give the vertices names
	int vertex_name = 0;
	Graph::vertex_iterator v, v_end;
	for (boost::tie(v,v_end) = boost::vertices(g); v != v_end; ++v) {
		boost::put(boost::vertex_color_t(), g, *v, vertex_name++);
	}

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
				bool exists_ab = boost::edge(boost::vertex(pair_table[open],g), boost::vertex(pos,g), g).second;
				bool exists_ba = boost::edge(boost::vertex(pos,g), boost::vertex(pair_table[open],g), g).second;
				if (!exists_ab && !exists_ba) {
					// add edge
					boost::add_edge(boost::vertex(pair_table[open],g), boost::vertex(pos,g), g);
				}
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
	
	return g;
}

