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
	
	std::vector < std::vector < int > > testcase;
	testcase.push_back({G,G,2});
	testcase.push_back({X,G,3});
	testcase.push_back({C,X,4});
	testcase.push_back({X,X,1});
	testcase.push_back({X,X,5});
	testcase.push_back({-1,A,6});
	testcase.push_back({-1,G,6});
	testcase.push_back({-1,X,4});
	
	for (auto elem : testcase) {
		sequence.clear();
		*out << "-------------------------------------------------------------------" << std::endl;
		*out << "Sequence (" << enum_to_char(elem[0]) << "," << enum_to_char(elem[1]) << "," << elem[2] << "): " << std::endl;
		if (elem[0] >= 0) {
			*out << "Number of Sequenes: " << generate_path_seq (sequence, elem[0], elem[1], elem[2]) << std::endl;
		} else {
			*out << "Number of Sequenes: " << generate_cycle_seq (sequence, elem[1], elem[2]) << std::endl;
		}
		*out << sequence << std::endl;
	}

	return 0;
}
