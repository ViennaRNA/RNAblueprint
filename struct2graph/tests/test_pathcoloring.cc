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

//declare global variables
bool verbose = true;
std::string outfile;

//! main program starts here
int main(int ac, char* av[]) {

	std::ostream* out = &std::cout;					// out stream
	
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
	
	Pairing pairing(6);						// init pairing matrix
	std::string sequence;
	*out << "Number of sequences: " << generate_path_seq (sequence, C, G, 14) << std::endl;
	*out << sequence << std::endl;
	
	return 0;
}
