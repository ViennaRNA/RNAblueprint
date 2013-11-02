/* This file is a boost.test unit test and provides tests for decompose.cc
*
*
* Created on: 23.10.2013
* Author: Stefan Hammer <s.hammer@univie.ac.at>
* License: GPLv3
*
* 
*
*/

// include header
#include "test_common.h"

// include headers containing functions to test
#include "../decompose.h"
#include "../printgraph.h"

// include std components

// include boost components
#include <boost/graph/iteration_macros.hpp>

// define heads

BOOST_AUTO_TEST_SUITE(Decompose)

enum type {CC, BC, ED};

Graph createGraph(int type) {
	Graph g(9);
	int vertex_name = 0;	
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		boost::put(boost::vertex_color_t(), g, v, vertex_name++);
	}
	
	// create a path between 0 and 10
	for (unsigned int i = 0; i < 7; i++) {
		boost::add_edge(boost::vertex(i,g), boost::vertex(i+1,g), g);
	}
	
	// create a block
	boost::add_edge(boost::vertex(0,g), boost::vertex(8,g), g);
	boost::add_edge(boost::vertex(0,g), boost::vertex(7,g), g);
	boost::add_edge(boost::vertex(2,g), boost::vertex(8,g), g);
	boost::add_edge(boost::vertex(6,g), boost::vertex(8,g), g);
	
	if ((type == BC) || (type == CC)) {
		// add a biconnected component
		for (int i = 9; i <=11; i++) {
			Vertex v = boost::add_vertex(g);
			boost::put(boost::vertex_color_t(), g, v, i);
		}
		boost::add_edge(boost::vertex(7,g), boost::vertex(9,g), g);
		boost::add_edge(boost::vertex(9,g), boost::vertex(10,g), g);
		boost::add_edge(boost::vertex(10,g), boost::vertex(11,g), g);
	}
	
	if (type == CC) {
		// add a connected component
		Vertex v = boost::add_vertex(g);
		boost::put(boost::vertex_color_t(), g, v, 12);
		v = boost::add_vertex(g);
		boost::put(boost::vertex_color_t(), g, v, 13);
		boost::add_edge(boost::vertex(12,g), boost::vertex(13,g), g);
	}
	
	/* full graph looks like this now (block plus path as biconnected component and connected component):
		5---6---7---9---10---11
		|   |   |
		4   8---0
		|   |   |
		3---2---1       12---13
	*/
	return g;
}

std::unordered_set<int> getVertexSet(Graph g) {
	std::unordered_set<int> result;
	BGL_FORALL_VERTICES(v, g, Graph) {
			result.insert(boost::get(boost::vertex_color_t(), g, v));
		}
	return result;
}


BOOST_AUTO_TEST_CASE(connectedComponents) {

	// create a graph
	Graph g = createGraph(CC);
	BOOST_TEST_MESSAGE("decompose connected components");
	connected_components_to_subgraphs(g);
	
	int number_of_children = 0;
	
	Graph::children_iterator child, child_end;
	for (boost::tie(child, child_end) = g.children(); child != child_end; ++child) {
		number_of_children++;
		// for the smaller connected component (12---13)
		if (boost::num_vertices(*child) == 2) {
			// check if both vertices exist and are labeled right
			BOOST_CHECK(boost::get(boost::vertex_color_t(), *child, 0) == 12);
			BOOST_CHECK(boost::get(boost::vertex_color_t(), *child, 1) == 13);
			// check if just one edge exists here
			BOOST_CHECK(boost::num_edges(*child) == 1);
		}
	}
	// check if it is just 2 connected components
	BOOST_CHECK(number_of_children == 2);
}

BOOST_AUTO_TEST_CASE(biconnectedComponents) {

	// get graph
	Graph g = createGraph(BC);
	// print_graph(g, &std::cout, "graph");
	BOOST_TEST_MESSAGE("decompose biconnected components 1");
	biconnected_components_to_subgraphs(g);
	
	int number_of_children = 0;
	
	Graph::children_iterator child, child_end;
	for (boost::tie(child, child_end) = g.children(); child != child_end; ++child) {
		number_of_children++;
		// print_graph(*child, &std::cout, "bc");
		// for the smaller connected component (12---13)
		if (boost::num_vertices(*child) == 4) {
			std::unordered_set<int> testcase {7, 9, 10, 11};			
			BOOST_CHECK(getVertexSet(*child) == testcase);
			// check if just one edge exists here
			BOOST_CHECK(boost::num_edges(*child) == 3);
		}
	}
	// check if it is just 2 connected components
	BOOST_CHECK(number_of_children == 2);
	
	// ---------------------------------------------------------------
	// do another one with 2 paths connected on the block
	Graph h = createGraph(BC);
	for (int i = 12; i <=13; i++) {
			Vertex v = boost::add_vertex(h);
			boost::put(boost::vertex_color_t(), h, v, i);
		}
	boost::add_edge(boost::vertex(0,h), boost::vertex(12,h), h);
	boost::add_edge(boost::vertex(12,h), boost::vertex(13,h), h);
	
	BOOST_TEST_MESSAGE("decompose biconnected components 2");
	biconnected_components_to_subgraphs(h);
	
	number_of_children = 0;
	
	for (boost::tie(child, child_end) = h.children(); child != child_end; ++child) {
		number_of_children++;
		// print_graph(*child, &std::cout, "bc");
		// for the smaller connected component (12---13)
		if (boost::num_vertices(*child) == 4) {
			std::unordered_set<int> testcase {7, 9, 10, 11};
			BOOST_CHECK(getVertexSet(*child) == testcase);
			// check if just one edge exists here
			BOOST_CHECK(boost::num_edges(*child) == 3);
		} else if (boost::num_vertices(*child) == 3) {
			std::unordered_set<int> testcase {0, 12, 13};
			BOOST_CHECK(getVertexSet(*child) == testcase);
			// check if just one edge exists here
			BOOST_CHECK(boost::num_edges(*child) == 2);
		}
	}
	// check if it is just 3 connected components
	BOOST_CHECK(number_of_children == 3);
}

BOOST_AUTO_TEST_CASE(schieberEarDecomposition) {

	// set random generator to a static seed;
	rand_gen.seed(1);
	
	// create a graph
	Graph g = createGraph(ED);
	BOOST_TEST_MESSAGE("decompose connected components");
	schieber_ear_decomposition(g);
	
	int number_of_children = 0;
	
	Graph::children_iterator child, child_end;
	for (boost::tie(child, child_end) = g.children(); child != child_end; ++child) {
		number_of_children++;
		// print_graph(*child, &std::cout, "ed");
		// for the smallest ear (6---8)
		if (boost::num_vertices(*child) == 2) {
			// check if both vertices exist and are labeled right
			std::unordered_set<int> testcase {6, 8};			
			BOOST_CHECK(getVertexSet(*child) == testcase);
			// check if just one edge exists here
			BOOST_CHECK(boost::num_edges(*child) == 1);
		// for the last cycle
		} else if (boost::num_vertices(*child) == 4) {
			// check if both vertices exist and are labeled right
			std::unordered_set<int> testcase1 {1, 0, 2, 8};			
			BOOST_CHECK(getVertexSet(*child) == testcase1);
			// check if just one edge exists here
			BOOST_CHECK(boost::num_edges(*child) == 4);
		}
	}
	// check if it is just 2 connected components
	BOOST_CHECK(number_of_children == 3);
}

BOOST_AUTO_TEST_CASE(partsBetweenArticulationPoints) {
	// TODO write a testfunction for partsbetween articulation points and fix this function!
}

BOOST_AUTO_TEST_SUITE_END()
