/* This file tests the generate_path_seq funktion of the file
* pathcoloring.h
*
* Created on: 21.08.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* compile with g++ -Wall -std=c++11 -lstdc++ -g -o test_pathcoloring ../common.cc ../printgraph.cc ../pathcoloring.cc test_pathcoloring.cc
*
*/

// include header
#include "../printgraph.h"
#include "../pathcoloring.h"
#include "../common.h"
#include <random>
#include <chrono>

//declare global variables
bool verbose = false;
std::string outfile;

//! main program starts here
int main(int ac, char* av[]) {

	std::ostream* out = &std::cout;					// out stream
	
	// random generator
	unsigned long seed = std::chrono::system_clock::now().time_since_epoch().count();
	rand_gen.seed(seed);
	std::cerr << "Using this seed: " << seed << std::endl;
	
	Graph g(10);
	
	int vertex_name = 0;	
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		boost::put(boost::vertex_color_t(), g, v, vertex_name++);
	}
	
	for (unsigned int i = 0; i < boost::num_vertices(g)-1; i++) {
		boost::add_edge(boost::vertex(i,g), boost::vertex(i+1,g), g);
	}
	
	// set ear_integer
	BGL_FORALL_EDGES_T(e, g, Graph) {
		g[e].ear = 1;
	}
	
	print_graph(g, out, "path");
	
	std::string sequence;
	*out << "Number of sequences (G,G,16): " << generate_path_seq (sequence, G, G, 16) << std::endl;
	*out << sequence << std::endl;
	
	sequence.clear();
	*out << "Number of sequences (X,G,3): " << generate_path_seq (sequence, X, G, 3) << std::endl;
	*out << sequence << std::endl;
	
	sequence.clear();
	*out << "Number of sequences (C,X,4): " << generate_path_seq (sequence, C, X, 4) << std::endl;
	*out << sequence << std::endl;
	
	sequence.clear();
	*out << "Number of sequences (X,X,1): " << generate_path_seq (sequence, X, X, 1) << std::endl;
	*out << sequence << std::endl;
	
	sequence.clear();
	*out << "Number of sequences (X,X,100): " << generate_path_seq (sequence, X, X, 100) << std::endl;
	*out << sequence << std::endl;
	
	*out << "And now two cycles:" << std::endl;
	sequence.clear();
	*out << "Number of sequences (A,10): " << generate_cycle_seq (sequence, A, 10) << std::endl;
	*out << sequence << std::endl;
	
	sequence.clear();
	*out << "Number of sequences (X,4): " << generate_cycle_seq (sequence, X, 4) << std::endl;
	*out << sequence << std::endl;
	return 0;
}
