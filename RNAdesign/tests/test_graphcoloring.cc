/* This file is a boost.test unit test and provides tests for pathcoloring.cc
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
#include "../graphcoloring.h"

// define heads
BOOST_AUTO_TEST_SUITE(GraphColoring)

class compound_test {
public:
    void test_construction() {
        BOOST_CHECK_THROW( new class_under_test( -1 ) );

        v = new class_under_test( 1 );

        BOOST_CHECK( v is valid );
        ...
    }

    void test_access_methods() {
        BOOST_CHECK_EQUAL( v->get_value(), 1 );
        ...
    }
private:
    class_under_test* v;
};
...

boost::shared_ptr<compound_test> instance( new compound_test );
Ts>add( BOOST_CLASS_TEST_CASE( &compound_test::constructor, instance ) );
Ts>add( BOOST_CLASS_TEST_CASE( &compound_test::test_access_methods, instance ) );


Sequence get_vertex_colors(Graph& g) {
	Sequence sequence;
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		sequence.push_back(g[v].base);
	}
	return sequence;
}

void reset (Graph& g) {
	// reset color
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		g[v].base = X;
	}
}

BOOST_AUTO_TEST_CASE(colorPathGraph) {
	// set random generator to a static seed;
	rand_gen.seed(1);
	// create a graph
	Graph g(10);
	int vertex_name = 0;	
	BGL_FORALL_VERTICES_T(v, g, Graph) {
		boost::put(boost::vertex_color_t(), g, v, vertex_name++);
	}
	for (unsigned int i = 0; i < boost::num_vertices(g)-1; i++) {
		boost::add_edge(boost::vertex(i,g), boost::vertex(i+1,g), g);
	}
	
	// color this graph!
	BOOST_TEST_MESSAGE("color path-graph");
	color_path_cycle_graph (g);
	//print_graph(g, out, "path");
	Sequence sequence {U, G, C, G, U, A, U, G, U, A};
	BOOST_CHECK(sequence == get_vertex_colors(g));
	
	// color first base and try all over again
	BOOST_TEST_MESSAGE("path_starts_A");
	reset(g);
	g[boost::vertex(0,g)].base = A;
	color_path_cycle_graph (g);
	Sequence sequence1 {A, U, G, U, G, U, G, U, G, C};
	BOOST_CHECK(sequence1 == get_vertex_colors(g));
	
	// color first base and try all over again
	BOOST_TEST_MESSAGE("path_ends_U");
	reset(g);
	g[boost::vertex(boost::num_vertices(g),g)].base = U;
	color_path_cycle_graph (g);
	Sequence sequence2 {C, G, U, G, U, A, U, G, U, A};
	BOOST_CHECK(sequence2 == get_vertex_colors(g));
	
	// color both ends and try all over again
	BOOST_TEST_MESSAGE("path_ends_AG");
	reset(g);
	g[boost::vertex(0,g)].base = A;
	g[boost::vertex(boost::num_vertices(g),g)].base = G;
	color_path_cycle_graph (g);
	Sequence sequence3 {A, U, G, C, G, U, G, U, A, U};
	BOOST_CHECK(sequence3 == get_vertex_colors(g));
	
	// make a cycle and color this cycle!
	BOOST_TEST_MESSAGE("cycle");
	reset(g);
	boost::add_edge(boost::vertex(0,g), boost::vertex(boost::num_vertices(g)-1,g), g);
	color_path_cycle_graph (g);
	Sequence sequence4 {C, G, C, G, C, G, C, G, U, G};
	BOOST_CHECK(sequence4 == get_vertex_colors(g));
	
	// set one base and color again this cycle
	BOOST_TEST_MESSAGE("cycle_starts_G");
	reset(g);
	g[boost::vertex(5,g)].base = G;
	color_path_cycle_graph (g);
	Sequence sequence5 {C, G, C, G, U, G, U, G, C, G};
	BOOST_CHECK(sequence5 == get_vertex_colors(g));
}

BOOST_AUTO_TEST_SUITE_END()
