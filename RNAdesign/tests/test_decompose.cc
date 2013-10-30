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

// include boost components
#include <boost/graph/iteration_macros.hpp>

// define heads
<<<<<<< HEAD
=======
namespace Decompose {
	class TestCase {
		public:
			TestCase(int first, int last, int length, int nos, Sequence sequence);
			int first;
			int last;
			int length;
			int nos;
			Sequence sequence;
	};
	
	Sequence get_vertex_colors(Graph& g);
	void reset (Graph& g);
}
>>>>>>> 30105234d0a84fdf40003a593b66b603e2079a67

BOOST_AUTO_TEST_SUITE(Decompose)

Sequence get_vertex_colors(Graph& g) {
	Sequence sequence;
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		sequence.push_back(g[v].base);
	}
	return sequence;
}

Graph createGraph() {
	Graph g(14);
	int vertex_name = 0;	
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		boost::put(boost::vertex_color_t(), g, v, vertex_name++);
	}
	
	// create a path between 0 and 10
	for (unsigned int i = 0; i < 10; i++) {
		boost::add_edge(boost::vertex(i,g), boost::vertex(i+1,g), g);
	}
	
	// create a block
	boost::add_edge(boost::vertex(0,g), boost::vertex(11,g), g);
	boost::add_edge(boost::vertex(2,g), boost::vertex(11,g), g);
	boost::add_edge(boost::vertex(6,g), boost::vertex(11,g), g);
	// add a connected component
	boost::add_edge(boost::vertex(12,g), boost::vertex(13,g), g);
	
	/* graph looks like this now (block plus path as biconnected component and connected component):
		5---6---7---8---9---10
		|   |   |
		4  11---0
		|   |   |
		3---2---1      12---13
	*/
	return g;
}


BOOST_AUTO_TEST_CASE(connectedComponents) {

	// create a graph
	Graph g = createGraph();
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

/*
	BOOST_TEST_MESSAGE("decompose biconnected components");
	biconnected_components_to_subgraphs(g);
	
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
*/
}

BOOST_AUTO_TEST_SUITE_END()
